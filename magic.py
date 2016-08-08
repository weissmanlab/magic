#!/usr/bin/env python3

import sys, math, scipy, scipy.stats, scipy.optimize, itertools
import numpy as np



# Plotting: this is just for Jupyter notebooks
if __name__ != "__main__":
	import matplotlib.pyplot as plt
	import matplotlib
	plt.style.use('seaborn-talk')
	matplotlib.rcParams['axes.labelsize'] = 20
	matplotlib.rcParams['axes.titlesize'] = 24
	matplotlib.rcParams['xtick.labelsize'] = 20
	matplotlib.rcParams['xtick.major.size'] = 10
	matplotlib.rcParams['xtick.minor.size'] = 5
	matplotlib.rcParams['ytick.labelsize'] = 20
	matplotlib.rcParams['ytick.major.size'] = 10
	cMAGIC = 'm'
	cMSMC = 'c'
	cTrue = 'k'
	cMLT = 'b'



class ProbPoint():
	def __init__(self, x=None, p=None, e=None):
		self.x = x
		self.p = p
		self.e = e
	def check(self, emax=None, pmin=0, pmax=1):
		if pmin < self.p < pmax:
			# require that the error bars be positive (don't allow points that claim to be known perfectly) but not too large:
			if math.isfinite(self.e) and 0 < self.e < 1: 
				if emax:
					if emax > self.e/self.p/(1-self.p):
						return True
				else:
					return True
		return False

			

# Building the LT curve

# estimate underlying coalescence times from a mutation histogram
# tries to get it right on a per-window basis, rather than matching overall distribution
def mhist2Tlist(hist, method='Ghoshetal5'):
	if method == 'Ghoshetal5':
		x1 = hist.nonzero()[0][0]
		if len(hist) > x1+1 and sum(hist[x1+1:]) > 0:
			Nm1oD = (sum(hist[x1+1:]) - 1) / (sum(j*k for j, k in enumerate(hist[x1:])) - 1)
			return np.array([(k, j - Nm1oD * max(j-x1-1, 0), j) for j, k in enumerate(hist) if k])
		else:
			return np.array([(hist[x1], x1, x1)])
	elif method == 'ClevensonZidek':
		tot = sum(j*k for j,k in enumerate(hist))
		if tot:
			corr = 1 + (np.count_nonzero(hist)-1) / tot
			return np.array([(k,j/corr,j) for j, k in enumerate(hist) if k])
		else:
			return np.array([(hist[0],0,0)])
	elif method == 'ML':
		return np.array([(k,j,j) for j, k in enumerate(hist) if k])
	else:
		raise ValueError('Error: invalid method')

# given a mutation histogram, return a generating function with error bars
def gfe(hist,method='Ghoshetal5'):
	'''given a mutation histogram, return a generating function with error bars'''
	tot = sum(hist)
	Tlist = mhist2Tlist(hist,method)
	#need to allow for fact that homozygous windows don't really have T=0:
	if Tlist[0,1] == 0:
		Tlist[0,1] = np.log(2)/Tlist[0,0] # could use something else
	def gfpme(z):
		zpows = np.power( z, range(len(hist)) )
		p0 = zpows @ hist/tot
		with np.errstate(over='raise'):
			try:
				err = np.sqrt((np.exp(-(1-z**2) * Tlist[:,1]) - np.exp(-2*(1-z) * Tlist[:,1])) @ Tlist[:,0]) / tot
			except FloatingPointError:
				err = np.inf
		return [p0, err]
	return gfpme

def probCheck(p, emax=None, pmin=0, pmax=1, var=False):
	'''Check that a putative probability with error makes sense'''
	# last two entries in p should be value and error
	# error is standard error, unless var=True, in which case it is variance
	if pmin < p[-2] < pmax:
		#require that the error bars be positive but not too large -- we don't allow points that claim to be known perfectly:
		if math.isfinite(p[-1]) and 0 < p[-1] < 1: 
			if emax:
				if var: #need to take sqrt to get std err
					err = math.sqrt(p[-1])
				else:
					err = p[-1]
				if emax > err/p[-2]/(1-p[-2]):
					return True
			else:
				return True
	return False


# only take points if they look like real probabilities with reasonable error bars
def probFilter(points, emax=None, pmin=0, pmax=1):
	return [p for p in points if (p and probCheck(p,emax,pmin,pmax))]
	# include the "if p" to filter out Nones

# sigmoid function
def sigmoid(x,yleft,yright,xmid,slope):
	'''Shifted, scaled sigmoid function'''
	return yleft + (yright-yleft) * scipy.special.expit( slope * ( xmid - x ) )


# fit a sigmoid curve to some points
def sigmoidFit(points):
	xvals, yvals, sigmavals = zip(*points)
	p0vals = [ min(yvals[0]*1.1,1.0), yvals[-1]*.9, .5*(xvals[0]+xvals[-1]), 0.1 ]
	if max(yvals) - min(yvals) > .1 * np.median(sigmavals) :
		try:
			return scipy.optimize.curve_fit( sigmoid,xvals, yvals, p0=p0vals, sigma=sigmavals, absolute_sigma=True, maxfev=10**5 )
		except Exception as error:
			# print(error)
			return None
	else:
		return None


def H0e(s, genFs, Lratio=2, baseL=80, coverage=1, emaxL=1, pmin=0, pmax=1, extrapolation=.5):
	'''Infer the pointwise LT with error at s given a list of window-averaged generating functions genf'''
	LTLpts = probFilter([[i] + f(1 - s / (coverage * baseL * Lratio**i)) for i, f in enumerate(genFs)], emax=emaxL, pmin=pmin, pmax=pmax)
	if len(LTLpts) < 4 : #not going to be able to fit sigmoid
		return None
	fit = sigmoidFit(LTLpts)
	if fit:
		# check that we have data close to left asymptote (ie, not extrapolating too much):
		if (fit[0][2] - LTLpts[0][0]) * np.abs(fit[0][3]) * extrapolation > 1:
			fit0 = sigmoid(0, *fit[0])
			#approximate the variance at 0:
			esx = np.exp(fit[0][2] * fit[0][3])
			var0 = (esx**2 * fit[1][0][0] + fit[1][1][1] + 2 * esx * fit[1][1][0]) / (1+esx)**2
			if probCheck([s, fit0, var0], var=True) :
				return [s, fit0, np.sqrt(var0)]
	return None
	
		
def inferSLT(counts, svals=None, sratio=np.sqrt(2), Lratio=2, baseL=80, coverage=1, maxHom=.99, emaxL=1, pmin=0, pmax=1, emax0=.1, extrapolation=.5, failtol=2):
	'''From diversity histograms across a range of window lengths, calculate the Laplace transform of the coalescence time distribution at a range of points'''
	# generating functions of histograms at every length scale:
	genFs = [gfe(hist) for hist in counts]
	# define function turn given s into [s, LT{p_T}(s), error]:
	def s2pt(s):
		return H0e(s, genFs, Lratio=Lratio, baseL=baseL, coverage=coverage, emaxL=emaxL, pmin=pmin, pmax=pmax, extrapolation=extrapolation)
	if svals:  # list of desired s values provided
		return np.array(probFilter([s2pt(s) for s in svals], emax=emax0))
	else:  # need to determine appropriate s values
	# start with a LT variable value that should have a nice curve
	s = 0.1 * coverage * baseL * sum(counts[0]) / (np.arange(len(counts[0])) @ counts[0])
	allSLT = [s2pt(s)]
	# go to lower values until estimated homozygosity exceeds maxHom:
	while allSLT[-1] and probCheck(allSLT[-1], pmax=maxHom, emax=emax0):
		s /= sratio
		allSLT.append(s2pt(s))
	allSLT = sorted(slt for slt in allSLT if slt)
	try:
		s = allSLT[-1][0]
	except: #no s values worked so far; just give up
		return None
		#s = sratio * coverage * baseL * sum(counts[0]) / ( np.arange( len(counts[0]) ) @ counts[0] )
	# go to higher values until we run out of data:
	failures = 0 
	while failures <= failtol: # we will allow for some gaps
		s *= sratio
		slt = s2pt(s)
		if slt and probCheck(slt, emax=emax0):
			allSLT.append(slt)
		else:
			failures += 1
	return np.array(probFilter(allSLT, emax=emax0))

	
# Inferring a mixture of gamma distributions:

def GammaObj(GParams, mLTobs, zeroPt=False):
	norm = np.sum(GParams[::3])
	if zeroPt:
		# last entry of GParams is weight at t=0
		return np.linalg.norm([((GParams[:-1:3] @ np.power(1 + GParams[2::3]*obs[0], -GParams[1::3]) + GParams[-1]) / norm - obs[1]) / obs[2] for obs in mLTobs])
	else:
		return np.linalg.norm([(GParams[::3] @ np.power(1 + GParams[2::3]*obs[0], -GParams[1::3]) / norm - obs[1]) / obs[2] for obs in mLTobs])

class GammaParamStep(object):
	def __init__(self, stepsize=0.5):
		self.stepsize = stepsize
	def __call__(self, gp):
		s = self.stepsize
		# component weights are drawn from a Dirichlet distribution to normalize them:
		gp[::3] = np.random.dirichlet(gp[::3]/sum(gp[::3])/s)
		# component means and scales are log-normal to keep them positive
		# we take steps such that the *median*, not mean, is equal to the starting value
		gp[1::3] = np.random.lognormal(np.log(gp[1::3]), s)
		gp[2::3] = np.random.lognormal(np.log(gp[2::3]), s)
		return gp
		

def InferGParams(mLTobs, method='basinhopping', zeroPt=False, fullout=False, guess=None, npieces=None, bndries=None, m=None, T=None, niter=100, factr=1e3, pgtol=1e-6, maxfun=1e4, maxiter=1e4):
	if guess is None:
		if zeroPt:
			if npieces is None:
				npieces = 1 + len(mLTobs) // 3
			guess = np.append( np.ravel( [ [ 1/npieces, npieces-.1-i, 1/mLTobs[-i,0]/(npieces-.1-i) ] for i in range(npieces-1) ] ) ,1/npieces )
		else:
			if npieces is None:
				npieces = (len(mLTobs) + 1) // 3
			guess = np.ravel([ [1/npieces, npieces-.1-i, 1/mLTobs[-i,0]/(npieces-.1-i)] for i in range(npieces)])
	if bndries is None:
		bndries = [ [ (0,1), (0,None), (0,None) ][x % 3] for x in range(len(guess)) ]
	if m is None:
		m = max( 10, len(mLTobs)**2 )
	if method=='basinhopping':
		step = GammaParamStep()
		if T is None:
			# we expect the differences among peak heights to scale with the number of points being fitted:
			T = len(mLTobs)/10
		ans = scipy.optimize.basinhopping(GammaObj, guess, niter=niter, T=T, take_step=step, minimizer_kwargs={"method":"L-BFGS-B", "args":(mLTobs,zeroPt), "bounds":bndries, "options": {"ftol":factr*1e-17, "gtol":pgtol, "maxfun":maxfun, "maxiter":maxiter, "maxcor":m}},)
		if fullout:
			return ans
		else:
			return [ans.x, ans.fun]
	else:
		return scipy.optimize.fmin_l_bfgs_b(GammaObj, guess, args=(mLTobs,zeroPt), approx_grad=True, bounds=bndries, factr=factr, pgtol=pgtol, maxfun=maxfun, maxiter=maxiter, m=m)



# Gamma mixture distributions:

class GammaMix(scipy.stats.rv_continuous):
	def __init__(self, params):
		scipy.stats.rv_continuous.__init__(self, a=0)
		self.params = np.copy(params)
		self.params[::3] /= np.sum(self.params[::3])
		self.parray = np.copy(self.params[:3*(len(self.params)//3)]) #leave off any trailing weight at t=0
		self.parray.shape = (len(self.parray)//3, 3)
	def _pdf(self, t):
		# Note that we just omit any possible point mass at t=0 (but it's included in normalization)
		return np.sum(self.params[i] * scipy.stats.gamma.pdf(t, self.params[i+1], scale=self.params[i+2]) for i in range(0, len(self.params)-1, 3))	
	def _cdf(self, t):
		#include possibility for delta functions at 0
		return np.sum(self.params[i] * scipy.stats.gamma.cdf(t, self.params[i+1], scale=self.params[i+2]) if (i<len(self.params)-2 and np.prod(self.params[i:i+3])) else self.params[i] for i in range(0,len(self.params),3))
	def lt(self, s):
		return np.sum(self.params[i] * np.power(1 + self.params[i+2]*s, -self.params[i+1]) if i<len(self.params)-2 else self.params[i] for i in range(0, len(self.params), 3))
	def blcdf(self, r):
		return (np.prod(self.parray, axis=1) / np.sum(np.prod(self.parray, axis=1))) @ np.power(1 + self.parray[:,2]*r, -self.parray[:,1] - 1)
	def ne(self, t):
		return self.sf(t)/self.pdf(t)
	def eg(self, trange):
		theta0 = (1-self.cdf(trange[0])) / self.pdf(trange[0])
		EG = [(t0 / theta0, np.log(self.pdf(t1)/self.pdf(t0)*(1-self.cdf(t0))/(1-self.cdf(t1))) / (t1-t0) * theta0) for t0, t1 in zip(trange, trange[1:])]
		EG.append((trange[-1] / theta0, 0))
		return EG
	def ms(self, trange, L=0, rho=0, trees=False):
		'''Produce parameter string for ms from gamma mixture parameters'''
		theta0 = (1-self.cdf(trange[0])) / self.pdf(trange[0])
		eGparams = self.eg(trange)
		msparams = ''
		if trees:
			msparams += '-T '
		if L:
			msparams += '-t {} -r {} {} -p {} '.format(L*theta0, L*theta0*rho, L, math.ceil(math.log10(L)))
		msparams += '-eG ' + ' -eG '.join(' '.join(str(param) for param in eg) for eg in eGparams)
		return msparams


# Processing windower output:

def combinecounts(counts, input="sparse"):
	'''Combine multiple histograms'''
	if input == "sparse":
		combokeys = set().union(*[hist.keys() for hist in counts])
		return {i:sum(hist[i] for hist in counts if i in hist.keys()) for i in combokeys}
	elif input == "full":
		total = np.zeros(max(len(hist) for hist in counts), dtype=np.int)
		for hist in counts:
			total[:len(hist)] += hist
		return total
	else:
		sys.exit("Unknown input format")


def dict2array(wdata):
	'''Given dictionary histogram of data {i:count_i}, returns an array with array[i]=count_i'''
	maxi = max(wdata.keys())
	a = np.zeros(maxi+1)
	for i,n in wdata.items():
		a[i] = n
	return a
	
	
def extractCounts(filenames, input="sparse"):
	'''from a list of files, returns a list of arrays with the j^th entry of the i^th array = number of windows at lengthscale i with j polymorphisms'''
	# if the input files are written as sparse dictionaries:
	if input == "sparse": 
		countdicts = []
		for file in filenames:
			with open(file, 'r') as infile:
				countdicts.append([{int(pair.split()[0]): int(pair.split()[1]) for pair in line.split(',')} for line in infile])
		combocountdicts = [combinecounts(hists, "sparse") for hists in itertools.zip_longest(*countdicts)]
		return [dict2array(hist) for hist in combocountdicts]
	#if the input files are written as full lists:
	elif input == "full":
		allcounts = []
		for file in filenames:
			with open(file, 'r') as infile:
				allcounts.append([np.array([int(x) for x in line.split()]) for line in infile])
		return [combinecounts(hists, "full") for hists in itertools.zip_longest(*allcounts)]
	else:
		sys.exit("Unknown input format")
	
	
	
# code to run as script:
if __name__ == "__main__":
	import argparse
	
	parser = argparse.ArgumentParser()
	parser.add_argument("countfiles", nargs="+", help="files with histograms of polymorphisms/window")
	parser.add_argument("--outpath", help="Output prefix (otherwise prints to stdout)")
	parser.add_argument("--baselength", help="Number of (called) bases in shortest windows", type=np.int, default=64)
	parser.add_argument("--maxLT", help="Max value of Laplace transform to fit", type=np.float, default=.99)
	parser.add_argument("--extrapolation", help="How far to extrapolate to small length scales. 1 is a lot, .1 is very little.", type=np.float, default=.5)
	parser.add_argument("--zero", help="Allow the coalescence time to be exactly 0 with some probability", action='store_true')
	parser.add_argument("--LT", nargs="?", help="Also return the Laplace transform values. Add `only' to return *only* the LT values.", default=False, const=True)
	parser.add_argument("--components", help="Number of components to fit in gamma mixture", type=int, default=None)
	parser.add_argument("--iterations", help="How many times to run optimization algorithm", type=int, default=100)
	parser.add_argument("--maxfun", help="Max number of function evaluations in each optimization run", type=int, default=1e4)
	parser.add_argument("--input", help="Format of input histograms (full or sparse)", choices=("full","sparse"), default="sparse")
	args = parser.parse_args()

	counts = extractCounts(args.countfiles, args.input)
# 	with open(args.outpath + '_totcounts.txt', 'w') as outfile:
# 		print('\n'.join(' '.join(str(x) for x in hist) for hist in counts), file=outfile)
	
	SLTpts = inferSLT(counts, baseL=args.baselength, maxHom=args.maxLT, extrapolation=args.extrapolation)
	
	if args.LT:
		LTstring = '\n'.join(' '.join(str(x) for x in pt) for pt in SLTpts)
		if args.outpath:
			with open(args.outpath + '_LT.txt', 'w') as outfile:
				print(LTstring, file=outfile)
		else:
			print(LTstring + '\n')
		if args.LT == 'only':
			sys.exit()

	GParams = InferGParams(SLTpts, zeroPt=args.zero, npieces=args.components, niter=args.iterations, maxfun=args.maxfun, fullout=True)
	
	if args.outpath:
		with open(args.outpath + '_full.txt', 'w') as outfile:
			print(GParams, file=outfile)
	else:
		print(GParams)
		print('\n')
	
	try:
		GP = GParams.x # syntax if we're getting full output from basinhopping
	except:
		GP = GParams[0] # syntax if we just did L-BFGS-B or summarized output
	norm = np.sum(GP[::3])
	niceparams = [ GP[i:i+3]/[norm,1,1] for i in range(0, len(GP)-1, 3) if all(GP[i:i+3])]
	#niceparams doesn't include any mass at zero
	#calculate it separately:
	wzero = 0
	for i in range(0, len(GP)-1, 3):
		if any(GP[i+1:i+3]==0):
			# this corresponds to a point mass at 0
			wzero += GP[i]/norm
	if args.zero:
		# assumed mass at 0
		wzero += GP[-1]/norm
	
	nicestring = '\n'.join(' '.join(str(x) for x in component) for component in niceparams)
	if wzero:
		nicestring += ' ' + str(wzero)
	
	if args.outpath:
		with open(args.outpath + '_final.txt', 'w') as outfile:
			print(nicestring, file=outfile)
	else:
		print(nicestring)

	
	
	
	
	
	
	
	
	
