#!/usr/bin/env python3

# modified version of generate_multihetstep.py from Stephan Schiffels' msmc-tools

import sys
import gzip
import string
import copy
import argparse
import io

class MaskIterator:
	def __init__(self, filename, negative=False):
		if filename.endswith(".gz"):
			self.file = io.TextIOWrapper(gzip.open(filename, "r"))
		else:
			self.file = open(filename, "r")
		self.eof = False
		self.lastPos = 1
		self.negative = negative
		self.readLine()

	def readLine(self):
		try:
			line = next(self.file)
			fields = line.strip().split()
			if len(fields) == 2:
				self.start = int(fields[0])
				self.end = int(fields[1])
			else:
				self.start = int(fields[1]) + 1
				self.end = int(fields[2])
		except StopIteration:
			self.eof = True
	
	def getVal(self, pos):
		assert pos >= self.lastPos
		self.lastPos = pos
		while pos > self.end and not self.eof:
			self.readLine()
		if pos >= self.start and pos <= self.end:
			return True if not self.negative else False
		else:
			return False if not self.negative else True

class MergedMask:
	def __init__(self, mask_iterators):
		self.maskIterators = mask_iterators

	def getVal(self, pos):
		return all((m.getVal(pos) for m in self.maskIterators))

class VcfIterator:
	def __init__(self, filename):
		self.file = io.TextIOWrapper(gzip.open(filename, "r"))
	
	def __iter__(self):
		return self
	
	def __next__(self):
		line = next(self.file)
		while line.startswith("#"):
			line = next(self.file)
		fields = line.strip().split()
		chrom = fields[0]
		pos = int(fields[1])
		alleles = [fields[3]]
		for alt_a in fields[4].split(","):
			alleles.append(alt_a)
		geno = fields[9][:3]
		phased = geno[1] == "|"
		return (chrom, pos, tuple(alleles), (int(geno[0]), int(geno[2])), phased)

class OrderedAlleles:
	def __init__(self):
		self.ordered_alleles = []
	
	def addGenotype(self, a1, a2, phasing):
		if len(self.ordered_alleles) == 0:
			self.ordered_alleles = [[a1, a2]]
			if not phasing and a1 != a2:
				self.ordered_alleles.append([a2, a1])
		else:
			new = []
			for o in self.ordered_alleles:
				new.append(o + [a1, a2])
				if not phasing and a1 != a2:
					new.append(o + [a2, a1])
			self.ordered_alleles = new

	def getPrint(self):
		if len(self.ordered_alleles[0]) == 2 or len(self.ordered_alleles) == 1:
			return ''.join(self.ordered_alleles[0])
		else:
			return ','.join([''.join(o) for o in self.ordered_alleles])

class JoinedVcfIterator:
	def __init__(self, filenames):
		self.vcfIterators = [VcfIterator(f) for f in filenames]
		self.current_lines = [next(v) for v in self.vcfIterators]
	
	def __iter__(self):
		return self
	
	def __next__(self):
		minIndices = self.getMinIndices()
		chrom = self.current_lines[minIndices[0]][0]
		pos = self.current_lines[minIndices[0]][1]
		ref = self.current_lines[minIndices[0]][2][0]
	
		ordered_alleles = OrderedAlleles()
		
		for i, l in enumerate(self.current_lines):
			if i not in minIndices:
				ordered_alleles.addGenotype(ref, ref, True)
			else:
				alleles = self.current_lines[i][2]
				geno = self.current_lines[i][3]
				phased = self.current_lines[i][4]
				ordered_alleles.addGenotype(alleles[geno[0]], alleles[geno[1]], phased)
				try:
					self.current_lines[i] = next(self.vcfIterators[i])
				except StopIteration:
					self.current_lines[i] = None
		return (chrom, pos, ordered_alleles.getPrint())
	
	def getMinIndices(self):
		activeLines = [(i, l) for i, l in enumerate(self.current_lines) if l]
		if len(activeLines) == 0:
			raise StopIteration
		if len(activeLines) == 1:
			return [activeLines[0][0]]
		else:
			minIndices = [activeLines[0][0]]
			minPos = activeLines[0][1][1]
			for a in activeLines[1:]:
				if a[1][1] == minPos:
					minIndices.append(a[0])
				if a[1][1] < minPos:
					minPos = a[1][1]
					minIndices = [a[0]]
			return minIndices
		

parser = argparse.ArgumentParser()
parser.add_argument("files", nargs="+", help="Input VCF files")
parser.add_argument("--snpfile", help="Output file for SNPs (instead of stdout)", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument("--logfile", help="Log file (instead of stderr)", type=argparse.FileType('w'), default=sys.stderr)
parser.add_argument("--coverfile", help="Output file for coverage of windows")
parser.add_argument("--masks", nargs="+", help="Apply masks in bed format. Should have one calling mask for each individual. Can also have additional masks for e.g. mappability or admixture")
parser.add_argument("--negative_masks", nargs="*", help="same as mask, but interpreted as negative mask, so places where sites should be excluded")
parser.add_argument("--min_pos", help="Data starts at this chrom position",type=int,default=1)
parser.add_argument("--baselength", help="minimum window size",type=int,default=80)
parser.add_argument("--phase", help="list all possible phases?",action='store_true')
parser.set_defaults(phase=False)
args = parser.parse_args()

nrIndidividuals = len(args.files)
nrHaplotypes = 2 * nrIndidividuals

print("generating msmc input file with {} haplotypes".format(nrHaplotypes), file=args.logfile)

joinedVcfIterator = JoinedVcfIterator(args.files)
maskIterators = []
if args.masks:
	for f in args.masks:
		print("adding mask: {}".format(f), file=args.logfile)
		maskIterators.append(MaskIterator(f))
if args.negative_masks:
	for nm in args.negative_masks:
		print("adding negative mask: {}".format(nm), file=args.logfile)
		maskIterators.append(MaskIterator(nm, True))

mergedMask = MergedMask(maskIterators)

def is_segregating(alleles):
	if len(set(alleles.partition(",")[0])) > 1:
		return True
	return False

pos = args.min_pos-1
nr_called = 0
nr_called_window=0

if args.coverfile:
	coverage = open(args.coverfile,'w')
	
for chrom, snp_pos, alleles in joinedVcfIterator:
	while pos < snp_pos:
		pos += 1
		if mergedMask.getVal(pos):
			nr_called += 1
			nr_called_window +=1
		if pos % 1000000 == 0:
			print("processing pos {}".format(pos), file=args.logfile)
		if pos % args.baselength == 0:
			if args.coverfile:
				print(nr_called_window, file=coverage)
			nr_called_window = 0
	if mergedMask.getVal(snp_pos):
		if is_segregating(alleles):
			if args.phase:
				print(chrom, snp_pos, nr_called, alleles, sep="\t", file=args.snpfile)
			else:
				print(chrom, snp_pos, nr_called, alleles.split(",")[0], sep="\t", file=args.snpfile)
			nr_called = 0
			
# read to the end of the window containing the last SNP:
if args.coverfile:
	while pos % args.baselength != 0:
		pos += 1
		if mergedMask.getVal(pos):
			nr_called_window += 1
		if pos % args.baselength == 0:
			print(nr_called_window, end="\n", file=coverage)
	coverage.close()
	
