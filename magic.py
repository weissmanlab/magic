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
	tableau20 = np.array([(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)])/255  
             
# Printing: utility function for optional output files		
def chooseprint(*objects, file=sys.stdout, method='w', **kwargs):
	'''Print to stdout, stderr, or a named file. If file, opens and closes the file.'''
	if file in (None, sys.stdout, sys.stderr):
		if 'end' in kwargs.keys():
			print(*objects, file=file, **kwargs)
		else:
			print(*objects, file=file, end = '\n\n', **kwargs)
	else:
		with open(file, method) as outfile:
			print(*objects, file=outfile, **kwargs)

# Inferring the Laplace transform curve:

class ProbPoint:
	'''A probability, optionally with standard error and ordinate.'''
	def __init__(self, p, e=0, x=None):
		self.p = p
		self.e = e
		if x is not None:
			self.x = x
	def pe(self):
		return self.p, self.e
	def xpe(self):
		return self.x, self.p, self.e
	def check(self, emax=None, pmin=0, pmax=1, var=False):
		'''Check that a putative probability with error makes sense'''
		# error is standard error, unless var=True, in which case it is variance
		if pmin < self.p < pmax:
			#require that the error bars be positive but not too large -- we don't allow points that claim to be known perfectly:
			if math.isfinite(self.e) and 0 < self.e < 1: 
				if emax:
					if var: #need to take sqrt to get std err
						err = math.sqrt(self.e)
					else:
						err = self.e
					if emax > err / self.p / (1-self.p):
						return True
				else:
					return True
		return False


class Histogram:
	def __init__(self, counts):
		self.counts = counts
		self.n_obs = sum(counts)
		self.tot_hits = np.arange(len(counts)) @ counts
		if self.n_obs:
			self.mean = self.tot_hits/self.n_obs

class SNPHistogram(Histogram):
	'''A histogram of SNP counts across windows'''
	def make_tarray(self, method='Ghoshetal5'):
		'''Quick estimate of window-averaged coalescence times, to use in calculating stochasticity in mutation accumulation'''
		hist = self.counts
		if method == 'Ghoshetal5':
			x1 = hist.nonzero()[0][0]
			if len(hist) > x1+1 and sum(hist[x1+1:]) > 0:
				Nm1oD = (sum(hist[x1+1:]) - 1) / (sum(j*k for j, k in enumerate(hist[x1:])) - 1)
				return np.array([(k, j - Nm1oD * max(j-x1-1, 0), j) for j, k in enumerate(hist) if k])
			else:
				return np.array([(hist[x1], x1, x1)])
		elif method == 'ClevensonZidek':
			if self.tot_hits:
				corr = 1 + (np.count_nonzero(hist)-1) / self.tot_hits
				return np.array([(k, j/corr, j) for j, k in enumerate(hist) if k])
			else:
				return np.array([(hist[0], 0, 0)])
		elif method == 'ML':
			return np.array([(k, j, j) for j, k in enumerate(hist) if k])
		else:
			raise ValueError('Error: invalid method')
	def __init__(self, counts, bases=None, coverage=1):
		Histogram.__init__(self, counts)
		self.bases = bases
		self.coverage = coverage
		self.tarray = self.make_tarray()
		# need to allow for fact that homozygous windows don't really have T=0:
		if self.tarray[0, 1] == 0:
			self.tarray[0, 1] = np.log(2)/self.tarray[0, 0] # could use something else
	def theta(self):
		'''Mean number of mutations per sequenced base.'''
		return self.mean / (self.bases * self.coverage)
	def gfe(self, z):
		'''Generating function with error bars.'''
		zpows = z**np.arange(len(self.counts))
		p0 = zpows @ self.counts / self.n_obs
		with np.errstate(over='raise'):
			try:
				err = np.sqrt((np.exp(-(1-z**2) * self.tarray[:,1]) - np.exp(-2 * (1-z) * self.tarray[:,1])) @ self.tarray[:,0]) / self.n_obs
			except FloatingPointError:
				err = np.inf
		return ProbPoint(p0, err, z)
	def lte(self, s):
		'''Estimated LaplaceTransform{p_T}(s), with error.'''
		return self.gfe(1 - s/(self.bases * self.coverage))
	def ltle(self, s):
		'''Estimated LaplaceTransform{p_T}(s), with error and log(#bases) as ordinate.'''
		lt = self.gfe(1 - s/(self.bases * self.coverage))
		lt.x = np.log(self.bases)
		return lt


def sigmoid(x, yleft, yright, xmid, slope, alpha=1):
	'''Sigmoid function (generalized logistic curve).'''
	return yleft + (yright - yleft) * scipy.special.expit(slope * (x - xmid))**alpha


def sigmoid_fit(points, anchor=None):
	''''Fit a sigmoid curve to points.'''
	xvals, yvals, sigmavals = zip(*points)
	# initial guesses:
	# bounded linear extrapolation for left asymptote:
	yleft0 = min(1, yvals[0] + (yvals[0] - yvals[1]) * xvals[0] / (xvals[1] - xvals[0]))
	p0vals = [yleft0, yvals[-1]*.9, np.mean(xvals), 1]
	if anchor is None:
		# make sure that there is some slope to detect above the noise:
		if max(yvals) - min(yvals) < .1 * np.median(sigmavals):
			return None
		else:
			bounds = (0, [1, 1, np.inf, np.inf]) 
			f = sigmoid
	else:
		# make sure that there is some slope to detect above the noise:
		if max(yvals) - anchor < .1 * np.median(sigmavals):
			return None
		else:
			# anchor the right asymptote of the sigmoid:
			del p0vals[1]
			p0vals.append(1)
			def f(x, yleft, xmid, slope, alpha):
				return sigmoid(x, yleft, anchor, xmid, slope, alpha)
			bounds = (0, [1, np.inf, np.inf, np.inf])
	try:
		return scipy.optimize.curve_fit(f, xvals, yvals, p0=p0vals, sigma=sigmavals, absolute_sigma=True, bounds=bounds, max_nfev=10**5)
	except Exception as error:
		# print(error)
		return None			



def h0e(LTLpts, extrapolation=.5, anchor=None):
	'''Infer the pointwise LT with error at s given an array of estimates based on different window sizes.'''
	if len(LTLpts) < 4 : #not going to be able to fit sigmoid
		return None
	fit = sigmoid_fit(LTLpts, anchor=anchor)
	if fit:
		if anchor:
			yleft, xmid, slope, alpha = fit[0]
		else:
			yleft, yright, xmid, slope = fit[0]
		# check that we have data close to left asymptote (ie, not extrapolating too much):
		if (xmid - LTLpts[0][0]) * slope * extrapolation > 1:
# 			# check that inferred asymptote is not that far from linear extrapolation:
# 			for pt1, pt2 in zip(LTLpts, LTLpts[1:]):
# 				if np.abs(pt1[1] - pt1[0] * (pt2[1]-pt1[1]) / (pt2[0]-pt1[0]) - yleft) < yleft * (1-yleft) * extrapolation / 1:
# 					break
# 			else:
# 				return None
			# check that we have a valid probability:
			try:
				pt = ProbPoint(yleft, np.sqrt(fit[1][0][0]))
			except:
				return None
			if pt.check():
				return pt
	return None
	
		
def infer_slt(counts, svals=None, sratio=np.sqrt(2), maxHom=.99, emaxL=1, pmin=0, pmax=1, emax0=.1, extrapolation=.5, failtol=2, anchor=True):
	'''From diversity histograms across a range of window lengths, calculate the Laplace transform of the coalescence time distribution at a range of points'''
	# define function h0e(s) = [s, LT{p_T}(s), error]:
	def s2pt(s):
		LTLpts = [pt.xpe() for pt in (hist.ltle(s) for hist in counts) if pt.check(emax=emaxL, pmin=pmin, pmax=pmax)]
		if anchor:
			pe = h0e(LTLpts, extrapolation=extrapolation, anchor=np.exp(-s*counts[-1].theta()))
		else:
			pe = h0e(LTLpts, extrapolation=extrapolation)
		if pe is None:
			return None
		else:
			pe.x = s
			return pe
	if svals is not None:  # list of desired s values provided
		allSLT = [s2pt(s) for s in svals]
		return np.array([slt for slt in allSLT if slt.check(emax=emax0)])
	else:  # need to determine appropriate s values
		# start with a LT variable value that should have a nice curve
		s = 0.1 / counts[0].theta()
		allSLT = [s2pt(s)]
		# go to lower values until estimated homozygosity exceeds maxHom:
		while allSLT[-1] and allSLT[-1].check(pmax=maxHom):
			s /= sratio
			allSLT.append(s2pt(s))
		allSLT = sorted([slt for slt in allSLT if slt], key=lambda pt: pt.x)
		try:
			s = allSLT[-1].x
		except: #no s values worked so far; just give up
			return None
		# go to higher values until we run out of data:
		failures = 0 
		while failures <= failtol: # we will allow for some gaps
			s *= sratio
			slt = s2pt(s)
			if slt and slt.check(emax=emax0):
				allSLT.append(slt)
			else:
				failures += 1
		return np.array([slt.xpe() for slt in allSLT if slt.check(emax=emax0)])

	
# Inferring a mixture of gamma distributions:

def gamma_obj(GParams, mLTobs, zeroPt=False):
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
		eps = 1e-8
		# component weights are drawn from a Dirichlet distribution to normalize them:
		gp[::3] = np.random.dirichlet(gp[::3]/sum(gp[::3])/s)
		# component means and scales are log-normal to keep them positive
		# we take steps such that the *median*, not mean, is equal to the starting value
		# add eps to make sure we don't try to try to take log of 0
		gp[1::3] = np.random.lognormal(np.log(gp[1::3] + eps), s)
		gp[2::3] = np.random.lognormal(np.log(gp[2::3] + eps), s)
		return gp
		

def infer_gamma_mix(mLTobs, method='basinhopping', zeroPt=False, fullout=False, guess=None, npieces=None, bndries=None, m=None, T=None, niter=100, factr=1e3, pgtol=1e-6, maxfun=1e4, maxiter=1e4):
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
		ans = scipy.optimize.basinhopping(gamma_obj, guess, niter=niter, T=T, take_step=step, minimizer_kwargs={"method":"L-BFGS-B", "args":(mLTobs,zeroPt), "bounds":bndries, "options": {"ftol":factr*1e-17, "gtol":pgtol, "maxfun":maxfun, "maxiter":maxiter, "maxcor":m}},)
		if fullout:
			return ans
		else:
			return [ans.x, ans.fun]
	else:
		return scipy.optimize.fmin_l_bfgs_b(gamma_obj, guess, args=(mLTobs,zeroPt), approx_grad=True, bounds=bndries, factr=factr, pgtol=pgtol, maxfun=maxfun, maxiter=maxiter, m=m)



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
		'''Produce parameter string for ms from gamma mixture parameters.'''
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

def combine_counts(counts, input="sparse"):
	'''Combine multiple histograms.'''
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
	'''Given dictionary histogram of data {i:count_i}, returns an array with array[i]=count_i.'''
	maxi = max(wdata.keys())
	a = np.zeros(maxi+1)
	for i,n in wdata.items():
		a[i] = n
	return a
	
	
def extract_counts(filenames, input="sparse"):
	'''from a list of files, returns a list of arrays with the j^th entry of the i^th array = number of windows at lengthscale i with j polymorphisms.'''
	# if the input files are written as sparse dictionaries:
	if input == "sparse": 
		countdicts = []
		for file in filenames:
			with open(file, 'r') as infile:
				countdicts.append([{int(pair.split()[0]): int(pair.split()[1]) for pair in line.split(',')} for line in infile])
		combocountdicts = [combine_counts(hists, "sparse") for hists in itertools.zip_longest(*countdicts)]
		return [SNPHistogram(dict2array(hist)) for hist in combocountdicts]
	#if the input files are written as full lists:
	elif input == "full":
		allcounts = []
		for file in filenames:
			with open(file, 'r') as infile:
				allcounts.append([np.array([int(x) for x in line.split()]) for line in infile])
		return [SNPHistogram(combine_counts(hists, "full")) for hists in itertools.zip_longest(*allcounts)]
	else:
		sys.exit("Unknown input format")
	
	
	
# code to run as script:
if __name__ == "__main__":
	import argparse
	
	parser = argparse.ArgumentParser()
	parser.add_argument("countfiles", nargs="+", help="files with histograms of polymorphisms/window")
	parser.add_argument("--outpath", help="Output prefix (otherwise prints to stdout)")
	parser.add_argument("--baselength", help="Number of bases in shortest windows", type=np.int, default=80)
	parser.add_argument("--coverage", help="Fraction of bases that are sequenced", type=np.float, default=0.8)
	parser.add_argument("--maxLT", help="Max value of Laplace transform to fit", type=np.float, default=.99)
	parser.add_argument("--extrapolation", help="How far to extrapolate to small length scales. 1 is a lot, .1 is very little.", type=np.float, default=.5)
	parser.add_argument("--zero", help="Allow the coalescence time to be exactly 0 with some probability", action='store_true')
	parser.add_argument("--LT", help="Set to False to hide Laplace transform values. Set to `only' to return *only* the LT values.", default=True)
	parser.add_argument("--components", help="Number of components to fit in gamma mixture", type=int, default=None)
	parser.add_argument("--iterations", help="How many times to run optimization algorithm", type=int, default=50)
	parser.add_argument("--maxfun", help="Max number of function evaluations in each optimization run", type=int, default=5e4)
	parser.add_argument("--input", help="Format of input histograms (full or sparse)", choices=("full","sparse"), default="sparse")
	args = parser.parse_args()

	outfiles = {}
	for key in ('LT', 'final', 'full'):
		if args.outpath:
			outfiles[key] = args.outpath + '_{}.txt'.format(key)
		else:
			outfiles[key] = sys.stdout
	
	counts = extract_counts(args.countfiles, args.input)
	for scale, count in enumerate(counts):
		count.bases = args.baselength * 2**scale
		count.coverage = args.coverage
			
	SLTpts = infer_slt(counts, maxHom=args.maxLT, extrapolation=args.extrapolation)
	
	if SLTpts is None:
		sys.exit("Unable to infer the Laplace transform. If you don't have any more data, you might want to try increasing the allowed extrapolation.")
	
	if args.LT:
		LTstring = '\n'.join(' '.join(str(x) for x in pt) for pt in SLTpts)
		chooseprint(LTstring, file=outfiles['LT'])
		if args.LT == 'only':
			sys.exit()

	GParams = infer_gamma_mix(SLTpts, zeroPt=args.zero, npieces=args.components, niter=args.iterations, maxfun=args.maxfun, fullout=True)
	
	chooseprint(GParams, file=outfiles['full'])
	
	try:
		GP = GParams.x # syntax if we're getting full output from basinhopping
	except:
		GP = GParams[0] # syntax if we just did L-BFGS-B or summarized output
	norm = np.sum(GP[::3])
	niceparams = [GP[i:i+3]/[norm,1,1] for i in range(0, len(GP)-1, 3) if all(GP[i:i+3])]
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
	
	chooseprint(nicestring, file=outfiles['final'])

	
	
	
	
	
	
	
	
	
