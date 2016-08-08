#!/usr/bin/env python3

import sys

def convert(path, outpath=None, chromL=10**8, infsites=False):		
	if outpath:
		outname = outpath
	else:
		outname = path
	SNPs = []
	run = 0
	datastarted = False
	with open(path+'_ms.out.txt','r') as infile:
		#get number of runs from infile:
		params=next(infile).split()
		samples = [''] * int(params[1])
		nruns = int(params[2])
		#collect all SNPs into one big superchromosome:
		for line in infile:
			if datastarted:
				if line.startswith('positions:'):
					SNPs.extend([int(chromL*(run+float(x))/nruns) for x in line.split()[1:]])
					run += 1
					samplenum = 0
				if line.rstrip().isdigit():
					samples[samplenum] += line.rstrip()
					samplenum += 1
			elif line.startswith('//'):
				datastarted = True
	with open(outname+'_msmc.in.txt','w') as outfile:
		prev=0
		for j, snp in enumerate(SNPs):
			if snp <= prev:
				if prev >= chromL:
					print('Pileup at end of chromsome')
					break
				elif infsites: #need to move the snp to the next site
					if snp == prev:
						print('Collision at '+str(snp), file=sys.stderr)
					else:
						print('Out of order: '+str(snp)+' after '+str(prev), file=sys.stderr)
					snp = prev+1
				else: #throw away additional mutations to already-mutated sites
					continue
			print('\t'.join(['chr1', str(snp), str(snp-prev), ''.join(sample[j] for sample in samples)]), file=outfile)
			prev=snp

if __name__ == '__main__':	
	import argparse

	parser = argparse.ArgumentParser()
	parser.add_argument("path", help="path/prefix for input/output files")
	parser.add_argument("--outpath", help="path/prefix for output files, if different than input")
	parser.add_argument("--chromL", help="length of full genome to output",type=int,default=pow(10,8))
	parser.add_argument("--infsites", help="Use an infinite sites model where all mutations occur at unique loci", action='store_true')
	args = parser.parse_args()

	convert(args.path, args.outpath, args.chromL, args.infsites)

# 
# 	SNPs = []
# 	run = 0
# 	datastarted = False
# 	with open(args.path+'_ms.out.txt','r') as infile:
# 		#get number of runs from infile:
# 		params=next(infile).split()
# 		samples = [''] * int(params[1])
# 		nruns = int(params[2])
# 		#collect all SNPs into one big superchromosome:
# 		for line in infile:
# 			if datastarted:
# 				if line.startswith('positions:'):
# 					SNPs.extend([int(args.chromL*(run+float(x))/nruns) for x in line.split()[1:]])
# 					run += 1
# 					samplenum = 0
# 				if line.rstrip().isdigit():
# 					samples[samplenum] += line.rstrip()
# 					samplenum += 1
# 			elif line.startswith('//'):
# 				datastarted = True
# 	with open(outname+'_msmc.in.txt','w') as outfile:
# 		prev=0
# 		for j, snp in enumerate(SNPs):
# 			if snp <= prev:
# 				if prev >= args.chromL:
# 					print('Pileup at end of chromsome')
# 					break
# 				elif args.infsites: #need to move the snp to the next site
# 					if snp == prev:
# 						print('Collision at '+str(snp), file=sys.stderr)
# 					else:
# 						print('Out of order: '+str(snp)+' after '+str(prev), file=sys.stderr)
# 					snp = prev+1
# 				else: #throw away additional mutations to already-mutated sites
# 					continue
# 			print('\t'.join(['chr1', str(snp), str(snp-prev), ''.join(sample[j] for sample in samples)]), file=outfile)
# 			prev=snp
# 			