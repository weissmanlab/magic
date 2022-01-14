#!/usr/bin/env python3

import sys, math, scipy, scipy.stats, scipy.optimize, itertools, argparse
import numpy as np
             
             
# Parsing arguments (needed to run as script):

def parse_args(arglist):
	parser = argparse.ArgumentParser()
	parser.add_argument("countfiles", nargs="+", help="Files with histograms of polymorphisms per window (or, if LT=start, file with Laplace transform values)")
	parser.add_argument("--out", help="Output prefix (otherwise prints to stdout)")
	parser.add_argument("--baselength", help="Number of bases in shortest windows", type=int, default=80)
	parser.add_argument("--coverage", help="Fraction of bases that are sequenced", type=float, default=0.8)
	parser.add_argument("--maxLT", help="Max value of Laplace transform to fit", type=float, default=.99)
	parser.add_argument("--ltstep", help="Max spacing between inferred Laplace transform values", type=float, default=0.05)
	parser.add_argument("--extrapolation", help="How far to extrapolate to small length scales. 10 is a lot, .1 is very little.", type=float, default=.5)
	parser.add_argument("--family", help="Parametric form to use for coalescence time distribution", choices=("pieceexp", "gammamix"), default="pieceexp")
	parser.add_argument("--zero", help="Allow the coalescence time to be exactly 0 with some probability", action='store_true')
	parser.add_argument("--LT", help="Set to False to hide Laplace transform values. Set to 'only' to return *only* the LT values. Set to 'start' to fit a distribution starting from an existing Laplace transform", default=True)
	parser.add_argument("--components", help="Number of components to fit in probability distribution", type=int, default=None)
	parser.add_argument("--iterations", help="How many times to run optimization algorithm", type=int, default=50)
	parser.add_argument("--maxfun", help="Max number of function evaluations in each optimization run", type=int, default=5e4)
	parser.add_argument("--input", help="Format of input histograms ('full' or 'sparse'). You should leave this as 'sparse'; it's just included for backwards compatibility.", choices=("full", "sparse"), default="sparse")
	parser.add_argument("--smoothing", help="For piecewise-exponential distributions: how much of a penalty to assess for changes in coalescence rates", type=float, default=1)
	args = parser.parse_args(arglist)
	
	if args.family == 'pieceexp' and args.zero:
		sys.exit('Sorry, currently --zero can only be used with --family gammamix.')
	
	return args

             
# Printing:	

def chooseprint(*objects, file=sys.stdout, method='w', **kwargs):
	'''Print to stdout, stderr, or a named file. If file, opens and closes the file.'''
	if file in (None, sys.stdout, sys.stderr):
		if 'end' in kwargs.keys():
			print(*objects, file=file, **kwargs)
		else:
			# By default, add an extra newline when printing to stdout or stderr:
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
		self.counts = np.array(counts)
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
		'''Estimated LaplaceTransform{p_T}(s), with error and log_2(#bases) as ordinate.'''
		lt = self.gfe(1 - s/(self.bases * self.coverage))
		lt.x = np.log2(self.bases)
		return lt


def sigmoid(x, yleft=1, yright=0, xmid=0, slope=1, alpha=1):
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

def check_fit(fit, data, extrapolation=.5, anchor=None):
	'''Check that fit isn't extrapolating too much, i.e., we have data close to the left asymptote, and that we have a valid probability'''
	if fit is None:
		return False
	if anchor:
		yleft, xmid, slope, alpha = fit[0]
		yright = anchor
	else:
		yleft, yright, xmid, slope = fit[0]
		alpha = 1
	# check that we have data close to left asymptote (ie, not extrapolating too much):
	if (yleft - sigmoid(data[0][0], yleft, yright, xmid, slope, alpha))**2 < 0.001 * extrapolation * yleft * (1 - yleft):
		# check that we have a valid probability:
		try:
			pt = ProbPoint(fit[0][0], np.sqrt(fit[1][0][0]))
		except:
			return False
		if pt.check():
			return True
	return False

def h0e(LTLpts, extrapolation=.5, anchor=None, return_full=False):
	'''Infer the pointwise LT with error at s given a list of estimates based on different window sizes.'''
	if len(LTLpts) < 4 : #not going to be able to fit sigmoid
		return None
	# first fitting to find threshold where LT approaches small-scale value:
	fit = sigmoid_fit(LTLpts, anchor=anchor)
	if fit:
		if anchor:
			yleft, xmid, slope, alpha = fit[0]
		else:
			yleft, yright, xmid, slope = fit[0]
		# if we have enough points on both sides of the midpoint, try to do second fitting with just the short-scale points to focus in on left asymptote:
		shortLTLpts = [ltl for ltl in LTLpts if ltl[0] < xmid]
		longLTLpts = [ltl for ltl in LTLpts if ltl[0] > xmid]
		if len(shortLTLpts) >= 4 and len(longLTLpts) > 2:
			# if we have lots of short-scale points, we can restrict even further:
			while len(shortLTLpts) > 6 and (xmid - shortLTLpts[-1][0]) * slope < 1:
				shortLTLpts.pop()
			tmp_fit = sigmoid_fit(shortLTLpts, anchor=anchor)
			if check_fit(tmp_fit, shortLTLpts, extrapolation, anchor):
				fit = tmp_fit
		# check that we have a plausible fit:
		if check_fit(fit, LTLpts, extrapolation, anchor):
			pt = ProbPoint(fit[0][0], np.sqrt(fit[1][0][0]))
			if return_full:
				return fit
			else:
				return pt
	return None
	
		
def infer_slt(counts, svals=None, sratio=np.sqrt(2), maxHom=.99, emaxL=1, pmin=0, pmax=1, emax0=.1, extrapolation=.5, failtol=2, anchor=True, ltstep=0.05, min_s_pow=0.2):
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
		return np.array([slt.xpe() for slt in allSLT if slt and slt.check(emax=emax0)])
	else:  # need to determine appropriate s values
		if sratio <= 1:
			sys.exit("Error: sratio must be > 1.")
		# start with a wide range of s values:
		svals = np.geomspace(0.01 / counts[0].theta(), 100 / counts[0].theta())
		allSLT = [s2pt(s) for s in svals]
		# go to lower s values until estimated homozygosity exceeds maxHom or things get too noisy:
		while allSLT[0] and allSLT[0].check(pmax=maxHom):
			s = allSLT[0].x / sratio
			allSLT.insert(0, s2pt(s))
		allSLT = [slt for slt in allSLT if slt]
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
		# remove any points where the error bars exceed our tolerance:
		SLT = [slt for slt in allSLT if slt.check(emax=emax0)]
		# go back and fill in gaps where LT changed a lot between successive points
		# limit the number of points that can be added between any pair to avoid possible infinite loops caused by bad points
		min_s_ratio = sratio**min_s_pow
		i = 0
		while i < len(SLT) - 1:
			while SLT[i].p > SLT[i+1].p + ltstep and SLT[i+1].x/SLT[i].x > min_s_ratio:
				# add a point in between them
				# make a list of possible s-values in case the first doesn't work
				smids = np.logspace(np.log10(SLT[i].x), np.log10(SLT[i+1].x), num=4, endpoint=False)[1:]
				# start from the middle and work our way out
				order = np.argsort(np.abs(2 * np.log(smids) - np.log(SLT[i].x) - np.log(SLT[i+1].x)))
				for j in order:
					slt = s2pt(smids[j])
					if slt and slt.check(emax=emax0):
						SLT.insert(i+1, slt)
						break
				else: # give up on filling this gap, move on to the next one
					break
			i += 1
		return np.array([slt.xpe() for slt in SLT if slt.check(emax=emax0)])



# Inferring distributions from their Laplace transforms

## Helper functions for inferring a piecewise exponential distribution

def piece_exp_obj(params, breaks, mLTobs, smoothing=0, zeroPt=False):
	if zeroPt:
		# first entry of params is weight at t=0:
		p0 = params[0]
		rates = params[1:]
	else:
		p0 = 0
		rates = params
	breakPs = (1 - p0) * np.exp(np.cumsum(np.concatenate(((0,), -np.diff(breaks) * rates[:-1]))))
	prefactors = [breakPs * np.exp(-m * breaks) / (1 + m/rates) for m in mLTobs[:,0]]
	postfactors = [np.concatenate( (-np.expm1(-(rates[:-1] + m) * np.diff(breaks)), (1,)) ) for m in mLTobs[:,0]]
	return sum( ((prefactors[i] @ postfactors[i] + p0 - obs[1]) / obs[2])**2 for i, obs in enumerate(mLTobs) ) + smoothing * np.linalg.norm(np.diff(np.log(rates)))**2

	
# class PieceExpStep(object):
# 	def __init__(self, stepsize=0.5):
# 		self.stepsize = stepsize
# 	def __call__(self, rates):
# 		s = self.stepsize
# 		eps = 1e-8
# 		# rates are log-normal to keep them positive
# 		# we take steps such that the *median*, not mean, is equal to the starting value
# 		# add eps to make sure we don't try to try to take log of 0
# 		newrates = np.random.lognormal(np.log(rates + eps), s)
# 		return newrates

	
## Helper functions for inferring a mixture of gamma distributions:

def gamma_obj(GParams, mLTobs, zeroPt=False):
	norm = np.sum(GParams[::3])
	if zeroPt:
		# last entry of GParams is weight at t=0
		return np.linalg.norm([((GParams[:-1:3] @ np.power(1 + GParams[2::3]*obs[0], -GParams[1::3]) + GParams[-1]) / norm - obs[1]) / obs[2] for obs in mLTobs])
	else:
		return np.linalg.norm([(GParams[::3] @ np.power(1 + GParams[2::3]*obs[0], -GParams[1::3]) / norm - obs[1]) / obs[2] for obs in mLTobs])

class GammaMixStep(object):
	def __init__(self, stepsize=0.5):
		self.stepsize = stepsize
	def __call__(self, gp):
		s = self.stepsize
		eps = 1e-8
		gpnew = np.copy(gp)
		# component weights are drawn from a Dirichlet distribution to normalize them:
		gpnew[::3] = np.random.dirichlet(gp[::3]/sum(gp[::3])/s)
		# component means and scales are log-normal to keep them positive
		# we take steps such that the *median*, not mean, is equal to the starting value
		# add eps to make sure we don't try to try to take log of 0
		gpnew[1::3] = np.random.lognormal(np.log(gp[1::3] + eps), s)
		gpnew[2::3] = np.random.lognormal(np.log(gp[2::3] + eps), s)
		return gpnew

## Main inference function

def infer_distribution(mLTobs, method='basinhopping', family='gammamix', zeroPt=False, guess=None, npieces=None, bounds=None, 
	m=None, T=None, smoothing=0, eps=1e-9, niter=50, factr=1e3, pgtol=1e-6, maxfun=1e4, maxiter=1e4):
	'''Infer a probability distribution from its estimated Laplace transform'''
	if family not in ['gammamix', 'pieceexp']:
		sys.exit("Unknown functional family for distribution.")
	if len(mLTobs) < 2:
		sys.exit("Unable to infer enough of the Laplace transform to invert.")
	if m is None:
		m = max(10, len(mLTobs)**2)
	if family == 'gammamix':
		func = gamma_obj
		if guess is None:
			if zeroPt:
				if npieces is None:
					npieces = 1 + len(mLTobs) // 3
				guess = np.append(np.ravel([[ 1/npieces, npieces-.1-i, 1/mLTobs[-i,0]/(npieces-.1-i)] for i in range(npieces-1)]) , 1/npieces)
			else:
				if npieces is None:
					npieces = (len(mLTobs) + 1) // 3
				guess = np.ravel([[1 / npieces, npieces - .1 - i, 1 / mLTobs[-i,0] / (npieces-.1-i)] for i in range(npieces)])
		if bounds is None:
			bounds = [[(eps, 1), (eps, None), (eps, None)][x % 3] for x in range(len(guess))]
			# eps is a small number to keep quantities that should be positive from being set exactly to 0
		args = (mLTobs, zeroPt)
		step = GammaMixStep()
	elif family == 'pieceexp':
		func = piece_exp_obj
		if guess is None:
			if npieces is None:
				npieces = len(mLTobs)
			# try to infer out to time where ~95% of genome has coalesced:
			tmax = 1 / mLTobs[np.searchsorted(1 - mLTobs[:,1], 0.05), 0]
			breaks = np.concatenate( ((0,), np.geomspace(0.5/mLTobs[-1,0], tmax, num=npieces-1)) )
			guess = np.ones(npieces) * mLTobs[0,0]/(1 - mLTobs[0,1])
		if bounds is None:
			bounds = [(eps, None) for rate in guess]
			# eps is a small number to keep quantities that should be positive from being set exactly to 0
		args = (breaks, mLTobs, smoothing)
		step = None
	if method == 'basinhopping':
		if T is None: # need to set the "temperature" of the algorithm
			# we expect the differences among peak heights to scale with the number of points being fitted:
			T = len(mLTobs)/10
		ans = scipy.optimize.basinhopping(func, guess, niter=niter, T=T, take_step=step, minimizer_kwargs={"method":"L-BFGS-B", 
			"args":args, "bounds":bounds, "options": {"ftol":factr*1e-17, "gtol":pgtol, "maxfun":maxfun, "maxiter":maxiter, "maxcor":m}},)
		if family == 'pieceexp':
			ans.x = np.stack((breaks, ans.x)) # return the breakpoints along with the rates
	else:
		ans = list(scipy.optimize.fmin_l_bfgs_b(func, guess, args=args, approx_grad=True, bounds=bounds, factr=factr, pgtol=pgtol, maxfun=maxfun, maxiter=maxiter, m=m))
		if family == 'pieceexp':
			ans[0] = np.stack((breaks, ans[0])) # return the breakpoints along with the rates
	return ans
	
def clean_parameters(ans, args):
	'''Extract parameter values of a distribution from the output of infer_distribution'''
	try:
		params = ans.x # if we're getting full output from basinhopping
	except:
		params = ans[0] # if we just did L-BFGS-B
	if args.family == 'gammamix':
		norm = np.sum(params[::3])
		niceparams = [params[i:i+3]/[norm,1,1] for i in range(0, len(params)-1, 3) if all(params[i:i+3])]
		#niceparams doesn't include any mass at zero yet; calculate it separately:
		wzero = 0
		for i in range(0, len(params)-1, 3):
			if any(params[i+1:i+3]==0):
				# this corresponds to a point mass at 0
				wzero += params[i]/norm
		if args.zero:
			# assumed mass at 0
			wzero += params[-1]/norm
		if wzero:
			niceparams += [wzero,]
	elif args.family == 'pieceexp':
		niceparams = np.transpose(params)
	return niceparams
		
	


# Probability distributions:

## Gamma mixture distributions:

class GammaMix(scipy.stats.rv_continuous):
	'''Mixture of gamma distributions and optionally a discrete mass at 0.'''
	def __init__(self, params):
		scipy.stats.rv_continuous.__init__(self, a=0) # says that distribution is bounded below by 0
		self.params = np.ravel(np.copy(params))
		self.params[::3] /= np.sum(self.params[::3])
		self.parray = np.copy(self.params[:3*(len(self.params)//3)]) #leave off any trailing weight at t=0
		self.parray.shape = (len(self.parray)//3, 3)
	def _pdf(self, t):
		# Note that we just omit any possible point mass at t=0 (but it's included in normalization)
		return np.sum(self.params[i] * scipy.stats.gamma.pdf(t, self.params[i+1], scale=self.params[i+2]) for i in range(0, len(self.params)-1, 3))	
	def _cdf(self, t):
		#include possibility for delta functions at 0
		return np.sum(self.params[i] * scipy.stats.gamma.cdf(t, self.params[i+1], scale=self.params[i+2]) if (i<len(self.params)-2 and np.prod(self.params[i:i+3])) else self.params[i] for i in range(0,len(self.params),3))
	def _sf(self, t):
		return 1 - self._cdf(t)
	def lt(self, s):
		'''Laplace transform evaluated at s.'''
		return np.sum(self.params[i] * np.power(1 + self.params[i+2]*s, -self.params[i+1]) if i<len(self.params)-2 else self.params[i] for i in range(0, len(self.params), 3))
	def blcdf(self, r):
		'''Fraction of IBD blocks with map length less than r/(mutation rate).'''
		return np.sum(np.prod(component) * np.power(1 + r*component[2], -component[1] - 1) for component in self.parray) / np.sum(np.prod(self.parray, axis=1))
		# return (np.prod(self.parray, axis=1) / np.sum(np.prod(self.parray, axis=1))) @ np.power(1 + self.parray[:,2]*r, -self.parray[:,1] - 1)
	def ne(self, t):
		'''Inverse hazard rate ("effective population size" for pairwise coalescence time, but note that it is 4 * mu * N_e(2 * mu * t)).'''
		return self.sf(t)/self.pdf(t)
	def ms(self, trange=None, points=100, L=0, rho=0, trees=False):
		'''Produce parameter string for ms from gamma mixture parameters.'''
		if trange is None:
			trange = self.ppf(np.linspace(1/(points+1), 1, points, endpoint=False))
		theta0 = -trange[0] / np.log1p(-self.cdf(trange[0]))
		Nparams = [(trange[0] / theta0, self.ne(trange[0]) / theta0)]
		Nparams.extend([(t0 / theta0, np.log(self.ne(t0)/self.ne(t1)) / (t1-t0) * theta0) for t0, t1 in zip(trange, trange[1:])])
		Nparams.append((trange[-1] / theta0, 0))
		msparams = ''
		if trees:
			msparams += '-T '
		if L:
			msparams += '-t {} -r {} {} -p {} '.format(L*theta0, L*theta0*rho, L, math.ceil(math.log10(L)))
		msparams += '-eN ' + ' -eG '.join(' '.join(str(param) for param in eg) for eg in Nparams)
		return msparams

## Piecewise exponential distributions:
		
class PiecewiseExponential(scipy.stats.rv_continuous):
	"Piecewise-exponential probability distribution"
	def __init__(self, breaks, rates):
		scipy.stats.rv_continuous.__init__(self, a=0) # says that distribution is bounded below by 0
		self.breaks = np.copy(breaks)
		self.rates = np.copy(rates)
		if len(rates) != len(breaks):
			if len(rates) == len(breaks) + 1 and breaks[0] > 0:
				self.breaks = np.concatenate( ((0,), self.breaks) )
			else:
				raise Exception("Breaks and rates must match")
		self.breakPs = np.exp(np.cumsum(np.concatenate(((0,), -np.diff(self.breaks) * self.rates[:-1])))) # the survival function evaluated at the breakpoints
	def _pdf(self, t):
		i = np.searchsorted(self.breaks, t, side='right') - 1
		return self.rates[i] * self.breakPs[i] * np.exp(-self.rates[i] * (t - self.breaks[i]))
	def _sf(self, t):
		i = np.searchsorted(self.breaks, t, side='right') - 1
		return self.breakPs[i] * np.exp(-self.rates[i] * (t - self.breaks[i]))
	def _cdf(self, t):
		return 1 - self._sf(t)
	def mean(self):
		return -np.concatenate((np.diff(self.breakPs), [0])) @ (1/self.rates)
	def lt(self, s):
		'''Laplace transform evaluated at s.'''
		return self.breakPs[-1] * np.exp(-s * self.breaks[-1]) / (1 + s / self.rates[-1]) - (
			np.sum(self.breakPs[i] * np.exp(-s * self.breaks[i]) / (1 + s / self.rates[i]) * np.expm1(-(self.rates[i] + s) * gap) for i, gap in enumerate(np.diff(self.breaks))) )
	def blcdf(self, r):
		'''Fraction of IBD blocks with map length less than r/(mutation rate).'''
		rc = np.outer(r, np.ones(len(self.rates))) + self.rates
		prefactors = self.breakPs * np.exp(-np.outer(r, self.breaks)) / (1 + np.outer(r, 1/self.rates))
		postfactors = self.breaks + 1 / rc
		postfactors[:, :-1] -= (1 / rc[:, :-1] + self.breaks[1:]) * np.exp(-rc[:, :-1] * np.diff(self.breaks))
		return np.array([pre @ postfactors[i] for i, pre in enumerate(prefactors)]) / self.mean()
	def ne(self, t):
		'''Inverse hazard rate ("effective population size" for pairwise coalescence time, but note that it is 4 * mu * N_e(2 * mu * t)).'''
		return 1 / self.rates[np.searchsorted(self.breaks, t, side='right') - 1]
	def ms(self, L=0, rho=0, trees=False):
		'''Produce parameter string for ms.'''
		theta0 = 1 / self.rates[0]
		msparams = ''
		if trees:
			msparams += '-T '
		if L:
			msparams += '-t {} -r {} {} -p {} '.format(L * theta0, L * rho * theta0, L, math.ceil(math.log10(L)))
		msparams += '-eN ' + '-eN '.join('{} {} '.format(self.breaks[i] / theta0 / 2, 1 / self.rates[i] / theta0) for i in range(1, len(self.breaks)))
		return msparams
		


# Processing windower output:

def combine_counts(counts, input="sparse"):
	'''Combine multiple histograms.'''
	if input == "sparse":
		combokeys = set().union(*[hist.keys() for hist in counts if hist is not None])
		return {i:sum(hist[i] for hist in counts if hist is not None and i in hist.keys()) for i in combokeys}
	elif input == "full":
		total = np.zeros(max(len(hist) for hist in counts), dtype=int)
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
				countdicts.append([None if line=='\n' else {int(pair.split()[0]): int(pair.split()[1]) for pair in line.split(',')} for line in infile])
				# remove any trailing Nones arising from trailing whitespace in the input file:
				while countdicts[-1] and countdicts[-1][-1] is None:
					countdicts[-1].pop()
		combocountdicts = [combine_counts(hists, "sparse") for hists in itertools.zip_longest(*countdicts)]
		if {} in combocountdicts:
			sys.exit("There appears to be a lengthscale with no data. Try checking the *_counts.txt files for extraneous blank lines.")
		return [SNPHistogram(dict2array(hist)) for hist in combocountdicts]
	# if the input files are written as full lists:
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
	import warnings
	
	args = parse_args(sys.argv[1:])

	# Set up the output:
	if args.out:
		outfiles = {key: args.out + '_{}.txt'.format(key) for key in ('LT', 'final')}
		outfiles['log'] = args.out + '.log'
		# redirect warnings to the logfile:
		def warn2file(message, category, filename, lineno, file=None, line=None):
			'''Print warning to a file instead of sys.stderr'''
			with open(args.out+'.log', 'a') as warnfile:
				print(warnings.formatwarning(message, category, filename, lineno), file=warnfile)
		warnings.showwarning = warn2file
		# log the initial command:
		chooseprint(' '.join(sys.argv), file=outfiles['log'])
	else:
		outfiles = {key: sys.stdout for key in ('LT', 'final')}
		outfiles['log'] = sys.stderr

	if args.LT == 'start':
		# Import the Laplace transform:
		if len(args.countfiles) != 1:
			sys.exit('Error: please specify exactly one file for the Laplace transform')
		with open(args.countfiles[0], 'r') as infile:
			SLTpts = np.array([[float(x) for x in line.split()] for line in infile])
	else:
		# Import the diversity histograms:
		counts = extract_counts(args.countfiles, input=args.input)
		for scale, count in enumerate(counts):
			count.bases = args.baselength * 2**scale
			count.coverage = args.coverage
			
		# Infer the Laplace transform from the diversity histograms:
		SLTpts = infer_slt(counts, maxHom=args.maxLT, extrapolation=args.extrapolation, ltstep=args.ltstep)
	
		if SLTpts is None:
			sys.exit("Unable to infer the Laplace transform. If you don't have any more data, you might want to try increasing the allowed extrapolation.")
	
		# Output the Laplace transform:
		if args.LT:
			LTstring = '\n'.join(' '.join(str(x) for x in pt) for pt in SLTpts)
			chooseprint(LTstring, file=outfiles['LT'])
			if args.LT == 'only':
				sys.exit()
			
	# Infer the distribution from the Laplace transform:
	full_params = infer_distribution(SLTpts, zeroPt=args.zero, family=args.family, npieces=args.components, niter=args.iterations, maxfun=args.maxfun, smoothing=args.smoothing)
	
	# Save the full optimization result
	chooseprint(full_params, file=outfiles['log'], method='a')
	
	# Output just the parameters of inferred distribution:
	chooseprint('\n'.join(' '.join(str(x) for x in component) for component in clean_parameters(full_params, args)), file=outfiles['final'])

	
	
	
	
	
	
	
	
	
