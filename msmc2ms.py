#!/usr/bin/env python3

# takes a PSMC' (msmc) output file as input and outputs an input file for ms

def rho(prefix):
	'''finds the estimated mutation and recombination rates and returns the ratio'''
	with open(prefix + '_out.log', 'r') as infile:
		for line in infile:
			if line.startswith('mutationRate'):
				m0 = float(line.split()[-1])
				break
	with open(prefix + '_out.loop.txt', 'r') as infile:
		for line in infile:
			pass
	r = float(line.split()[0])
	return r/m0


# code to run as script:
if __name__ == "__main__":
	import sys, argparse, math

	parser = argparse.ArgumentParser()
	parser.add_argument("msmc", help="path/prefix for MSMC files to take as input")
	parser.add_argument("--form", help="desired format of ms output", type=str, choices=['trees','snps','both'], default='both') 
	parser.add_argument("--chromL", help="Length of chromosome to simulate. Set to 0 to do single-locus sims.", type=int, default=10**6)
	args = parser.parse_args()
	

	MSNeT = []

	with open(args.msmc + '_out.final.txt','r') as infile:
		for line in infile:
			# skip the header:
			if not line.startswith('t'):
				MSMCparams = [float(x) for x in line.split()]
				# MSMCparams should be [index, mu t_left, mu t_right, 1/(2Ne(t) mu)]
				# ms needs [t_left/(4Ne(0)), Ne(t)/Ne(0)]
				if MSMCparams[0] == 0: 
					#the first line of data; only use it to normalize the rest (and the mutation and recombination rates)
					norm = MSMCparams[-1] # = 1/(2Ne(0) mu)
				else:
					#check if there has been a change in Ne from the previous time step:
					if MSNeT == [] or norm / MSMCparams[-1] != MSNeT[-1][-1]: 
						MSNeT.append([MSMCparams[1] * norm / 2.0, norm / MSMCparams[-1]])
	
	outstring = ''		
	if args.form in ('trees','both'):
		outstring += '-T '			
	if args.form in ('snps','both'):
		outstring += '-t {} '.format(2 * args.chromL / norm)
	if args.chromL:
		rho = rho(args.msmc)
		outstring += '-r {} {} -p {} '.format(2*args.chromL/norm*rho, args.chromL, math.ceil(math.log10(args.chromL)))
	outstring += '-eN ' + ' -eN '.join(' '.join(str(x) for x in timestep) for timestep in MSNeT)
	
	print(outstring)

