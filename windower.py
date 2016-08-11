#!/usr/bin/env python3

import sys, math, numpy, os.path, collections


def merge(windows, num):
	if num % 1:
		sys.exit('Ratio of succesive window sizes must be an integer')
	while len(windows) % num: #delete from the ends
	# to squeeze the most out of data, better to do chromosome arms separately, so that centromeric region also gets dropped
		if windows[0][1] < windows[-1][1]:
			windows.pop(0)
		else:
			windows.pop()
	return [numpy.sum(windows[i:i+num], axis=0) for i in range(0, len(windows), num)]

		
def snpdist(windows, L, args):
	dist = {}
	for win in windows:
		 if win[1] >= args.min_coverage * L :
		 # the window is sequenced at high enough coverage
		 	if args.sample_gaps == 'up' and win[1] < L and win[0]:
		 		# upsample the window
		 		count = win[0] + numpy.random.binomial(L-win[1], win[0]/L)
		 	elif args.sample_gaps == 'down' and win[1] > math.ceil(args.min_coverage * L)  and  win[0]:
		 		# downsample the window
		 		count = numpy.random.hypergeometric(win[0], win[1]-win[0], math.ceil(args.min_coverage*L)) 
		 	else:
		 		# record window as-is
		 		count = win[0]
		 	try:
		 		dist[count] += 1
		 	except KeyError:
		 		dist[count] = 1
	return dist
		
def initwindows(SNPs, basecoverage, args):
	L0 = args.baselength
	#initialize windows
	windowcounts = numpy.histogram(SNPs, bins=list(range(1, len(basecoverage)*L0+1, L0)))[0]
	return [numpy.array(x) for x in zip(windowcounts, basecoverage)]
		
if __name__ == '__main__':	
	import argparse

	parser = argparse.ArgumentParser()
	parser.add_argument("path", help="Path/prefix for SNP file, and maybe coverage and output files")
	parser.add_argument("--coverfile", help="Coverage file, if different from <path>_cover.txt")
	parser.add_argument("--outpath", help="path/prefix for output files, if different from input")
	parser.add_argument("--scales", help="number of length scales to use for windows", type=int, default=16)
	parser.add_argument("--ratio", help="Ratio between successive length scales. Must be integer", type=int, default=2)
	parser.add_argument("--baselength", help="Smallest length scale", type=int, default=80)
	parser.add_argument("--min_coverage", help="Minimum coverage for a window to be counted", type=float, default=0.8)
	parser.add_argument("--sample_gaps", help="For counts, up- or downsample windows to the same coverage", choices=['up','down',None], default='down')
	parser.add_argument("--stat", help="Statistic to calculate (tbl, individual tips, folded sfs component)", default='tbl')
	parser.add_argument("--input_form", help="Format of input", choices=['msmc','windows','reverse_windows'], default='msmc')
	parser.add_argument("--final_windows", help="Print out the longest windows", action='store_true')
	args = parser.parse_args()


	if args.path.endswith('.txt'):
		args.path = args.path[:-4]	
	if args.coverfile:
		covername = args.coverfile
	else:
		covername = args.path + '_cover.txt'
	if args.outpath:
		outname = args.outpath
	else:
		outname = args.path

	if args.input_form == 'msmc':
		with open(args.path+'.txt', 'r') as infile:
			samplesize = len(next(infile).split()[-1]) # might need this if we're calculating SFS or individuals' singletons
		with open(args.path+'.txt', 'r') as infile:
			if args.stat == 'tbl':
				# list of positions of every polymorphic site
				SNPs = [int(line.split()[1]) for line in infile]
			elif args.stat == 'tbl_alleles':
				if args.sample_gaps:
					sys.exit('Up/down-sampling only work with number of segregating sites, not number of alleles')
				# list of positions of every polymorphic site, with multiplicity = num alleles - 1
				SNPs = [int(line.split()[1]) for line in infile for _ in range(1, len(set(line.split()[-1])))]
			elif args.stat == 'indiv_tips':
				if args.final_windows:
					print("No final_windows option yet for indiv_singletons", file=sys.stderr)
				# break out singletons in each individual separately
				SNPs = {i:[] for i in range(0, samplesize//2)}
				for line in infile:
					alleles = line.split()[-1]
					for allele, num in collections.Counter(alleles).items():
						if num == 1:
							SNPs[alleles.find(allele)//2].append(int(line.split()[1]))
			else: #guess that we're trying to use sites with allele in n copies for some n
				try:
					copies = int(args.stat)
				except ValueError:
					sys.exit("Not sure what statistic to calculate")
				if copies not in range(1,samplesize):
					sys.exit("SFS component out of range")
				# list of positions of every site (without multiplicity) where an allele has n copies
				SNPs = [int(line.split()[1]) for line in infile if copies in collections.Counter(line.split()[-1]).values()]
		with open(covername, 'r') as coverfile:
			basecoverage = [int(x) for x in coverfile]
	elif args.input_form == 'windows':
		with open(args.path+'.txt', 'r') as infile:
			windows = [numpy.array([int(x) for x in line.split()[:2]]) for line in infile]
	elif args.input_form == 'reverse_windows':
		with open(args.path+'.txt', 'r') as infile:
			windows = [numpy.array([int(x) for x in reversed(line.split()[:2])]) for line in infile]



	#count diversity in windows:
	if args.stat == 'indiv_tips':
		windowsums = {indiv: [] for indiv in SNPs.keys()}
	else:
		windowsums = []
	with open(outname+'.log','w') as outfile:	
		for fold in range(args.scales):
			print('Calculating stats for lengthscale '+str(fold), file=outfile)
			L = args.baselength * pow(args.ratio,fold)
			if args.stat == 'indiv_tips':
				try:
					windows = {indiv: merge(indivwindows, args.ratio) for indiv, indivwindows in windows.items()}
				except NameError:
					windows = {indiv: initwindows(indivSNPs, basecoverage, args) for indiv, indivSNPs in SNPs.items()}
				for indiv, indivwindows in windows.items():
					if len(indivwindows) < 2:
						print('Reached chromosome length at length scale {} for individual {}'.format(fold, indiv), file=outfile)
						del windows[indiv]
					else:
						windowsums[indiv].append(snpdist(indivwindows, L, args))
				if not windows:
					break
			else:
				try:
					windows = merge(windows, args.ratio)
				except NameError: 
					windows = initwindows(SNPs, basecoverage, args)
				if len(windows) < 2:
					print('Reached chromosome length at length scale ' + str(fold), file=outfile)
					break
				windowsums.append(snpdist(windows, L, args)) 
			
	if args.stat == 'indiv_tips':
		for indiv, indivwindows in windowsums.items():
			with open(outname+'{}_counts.txt'.format(indiv),'w') as outfile:
				for scale in indivwindows:
					print(', '.join(' '.join(str(x) for x in pair) for pair in sorted(scale.items())), file=outfile)			
	else:
		with open(outname+'_counts.txt','w') as outfile:
			for scale in windowsums:
				print(', '.join(' '.join(str(x) for x in pair) for pair in sorted(scale.items())), file=outfile)
		if args.final_windows:
			with open(outname+'_windows'+str(fold)+'.txt','w') as outfile:
				print('\n'.join((' '.join(str(x) for x in window)) for window in windows), file=outfile)




	
	
	
	
	
	