#!/usr/bin/env python3

# takes a PSMC' (msmc) output file as input and prints a simple list of the minimal parameters

import sys, argparse

parser = argparse.ArgumentParser()
parser.add_argument("msmc", help="path/prefix for input MSMC files")
args = parser.parse_args()


TCs = []

with open(args.msmc + '_out.final.txt','r') as infile:
	# skip the header
	for line in infile:
		if not line.startswith('time'):
			MSMCparams = [float(x) for x in line.split()]
			# MSMCparams should be [index, mu t_left, mu t_right, 1/(2Ne(t) mu)]
			# we just want [mu t_left, 1/(2Ne(t) mu)]
			if TCs == [] or MSMCparams[-1] != TCs[-1][-1]: #there has been a pop size change from the previous time step
				TCs.append([MSMCparams[1], MSMCparams[-1]])
				
print('\n'.join(' '.join(str(x) for x in timestep) for timestep in TCs))


