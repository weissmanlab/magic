#!/usr/bin/env python3

import math, numpy, sys, os, argparse, io, scipy, itertools
from Bio import Phylo

numpy.set_printoptions(threshold=numpy.nan, linewidth=numpy.nan)

def extractTree(x):
	#x is a string in the format "[<length>](1:<height>,2:<height>);\n", with <length> an int and <height> a decimal
	temp = x.partition(']')
	length = int(temp[0].lstrip('['))
	height = float(temp[-1].partition(',')[0].partition(':')[-1])
	return [length, height]

def removeGhosts(run):
	#run is a list of trees in form [length,height]. this combines trees separated by ghost recombinations
	i = 0
	while i +1 < len(run) :
		while  i +1< len(run) and run[i+1][1] == run[i][1]:
			run[i][0] += run.pop(i+1)[0]
		i += 1

parser = argparse.ArgumentParser()
parser.add_argument("path", help="path/prefix for input/output files")
parser.add_argument("--outpath",help="path/prefix for output files, if different than input")
parser.add_argument("--form",help="format of ms output", type=str,choices=['trees','snps','both'], default='both')
parser.add_argument("--chromL",help="length of simulated chromosome (not the full genome length)",type=int,default=pow(10,7))
parser.add_argument("--baselength",help="minimum window size to use",type=int,default=80)
parser.add_argument("--scales",help="number of different length scales to use",type=int,default=16)
parser.add_argument("--ratio",help="Ratio between successive length scales. Must be integer for this simple algorithm to work",type=int,default=2)
parser.add_argument("--deghost", help = "Remove ghost recombination events", type=int, default=1)
args = parser.parse_args()

baselength=args.baselength
scales=args.scales
ratio=args.ratio


if args.outpath:
	outname = args.outpath
else:
	outname = args.path
# make sure that there's a directory for the output:
if not os.path.exists(outname):
	os.makedirs(outname)

# Get the simulation parameters:
with open(args.path + '.out.txt', 'r') as infile:
	firstline = next(infile).split()

samples = [''] * int(firstline[1])

if args.form:
	form = args.form
elif '-T' in firstline:
	if '-t' in firstline:
		form = 'both'
	else:
		form = 'trees'
elif '-t' in firstline:
	form = 'snps'
else:
	sys.exit('Need to specify format for ms results')
if form not in ('snps', 'trees', 'both'):
	sys.exit('Unknown format for ms results')


if form in ('trees', 'both'):
	#clear the times file if it already exists:
	with open(outname + '/times.txt', 'w') as timefile:
		pass
	with open(args.path + '.out.txt', 'r') as infile:
		numrun = -1 #start at -1 so first run is run 0
		for line in infile:
			if line.startswith('//'): #we've come to the beginning of the data for a run
				numrun += 1
				#print 'Calculating stats for run ' + str(numrun)
				temp = next(infile)
				run = []
				while temp.startswith('['): #go through the whole run
					run.append(temp)
					try:
						temp = next(infile)
					except StopIteration:
						break
				# now we've got the whole run read in
				# convert run from list of strings to list of trees:
				run = [extractTree(tree) for tree in run] #each tree is in format [length,height]
				if args.deghost == True:
					# remove ghost recombination events:
					removeGhosts(run)
				with open(outname + '/times.txt', 'a') as timefile:
					print('\n'.join(' '.join('{0:.5g}'.format(x) for x in tree) for tree in run), file=timefile)
					
if form in ('snps','both'):
	SNPs = []
	#collect all SNPs into one big superchromosome:
	with open(args.path+'.out.txt', 'r') as infile:
		numrun = -1 #start at -1 so first run is run 0
		datastarted = False
		for line in infile:
			if datastarted:
				if line.startswith('positions:'):
					temp=line.split()
					numrun += 1
					SNPs.extend([int(args.chromL*(numrun+float(x))) for x in temp[1:]])
					samplenum = 0
				if line.rstrip().isdigit():
					samples[samplenum] += line.rstrip()
					samplenum += 1
			else:
				if line.startswith('//'):
					datastarted = True
	supchromL = args.chromL * (numrun + 1)
	pi0 = len(SNPs) / float(supchromL)
	segsums = []
	if len(samples) > 2:
		loci = list(zip(*samples))
		hetsums = []
	for fold in range(scales):
		L = baselength * pow(ratio, fold)
		SNPhist = numpy.histogram(SNPs, numpy.arange(0, supchromL, L))[0]
		segsums.append(numpy.bincount(SNPhist))
		if len(samples) > 2:
			# list of SNPs at bdries of windows:
			breakSNPs = numpy.pad( numpy.cumsum( SNPhist ), (1,0), 'constant')
			windowhets = [sum(pair[0][i] != pair[1][i] for i in range(j,k)) for pair in zip(samples[::2], samples[1::2]) for j,k in zip(breakSNPs, breakSNPs[1:])]
			hetsums.append(numpy.bincount(windowhets))
	with open(outname + '/counts.txt', 'w') as countfile:
		print('\n'.join(' '.join(str(x) for x in segsums[fold]) for fold in range(scales)), file=countfile)
	if len(samples) > 2:
		with open(outname + '/hets.txt', 'w') as hetfile:
			print('\n'.join(' '.join(str(x) for x in hetsums[fold]) for fold in range(scales)), file=hetfile)

	
						

	
	
	
	
	
	
	
	
	
	
	