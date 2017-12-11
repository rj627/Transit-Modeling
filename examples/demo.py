import numpy as np
import matplotlib.pyplot as plt
import batman
import emcee
import corner

#intialize a transit model
def initialize_model(t, t0, per, rp, a, inc, ecc, w, u, limb_dark):
	params = batman.TransitParams()
	params.t0 = t0 				
	params.per = per		
	params.rp = rp		
	params.a = a				
	params.inc = inc			
	params.ecc = ecc		
	params.w = w				
	params.u = u	 	      	       	
	params.limb_dark = limb_dark      	

	model = batman.TransitModel(params, t)	

	return params, batman.TransitModel(params, t)			#return parameters and model objects 

#prior
def lnprior(theta):
	return 0.							#assumes all priors have uniform probability	

#likelihood function 
def lnlike(theta, params, model, t, flux, err):
	params.rp, params.t0  = theta[0], theta[1]			#update parameters
	lc = model.light_curve(params)
	residuals = flux - lc
	ln_likelihood = -0.5*(np.sum((residuals/err)**2 + np.log(2.0*np.pi*(err)**2)))

	return ln_likelihood

#posterior probability
def lnprob(theta, params, model, t, flux, err):
	lp = lnprior(theta)
	if not np.isfinite(lp):
		return -np.inf
	return lp + lnlike(theta, params, model, t, flux, err)

########################################################################################################

### generate a light curve and add photon noise	###
n = 50									#number of exposures
t = np.linspace(-0.02, 0.02, n)  			  		#exposure times

#specify orbital parameters and limb darkening
t0 = 0. 								#time of inferior conjunction 
per = 1.								#orbital period	
rp = 0.1								#planet radius (in units of stellar radii)
a = 15.									#semi-major axis (in units of stellar radii)
inc = 87.								#orbital inclination (in degrees)	
ecc = 0.								#eccentricity	
w = 90.									#longitude of periastron (in degrees) 
u = [0.2] 		 	       			 		#limb darkening coefficients
limb_dark = "linear"      						#limb darkening model

#initialize model and parameters
params, m = initialize_model(t, t0, per, rp, a, inc, ecc, w, u, limb_dark)

#calculate light curve and add noise
flux = m.light_curve(params)					
err = 300.e-6
flux = flux + np.random.normal(0, err, n)		

########################################################################################################

### fit the light curve with MCMC ###

#initial guesses for MCMC fit parameters
guess_rp, guess_t0 = 0.1, 0.						
theta = [guess_rp, guess_t0]

#initialize sampler
ndim, nwalkers = len(theta), 50			
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = (params, m, t, flux, err))
pos = [theta + 1e-5*np.random.randn(ndim) for i in range(nwalkers)]

#run mcmc
sampler.run_mcmc(pos,500)

########################################################################################################

### look at output ###

#make a pairs plot from MCMC output
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
fig = corner.corner(samples, labels = ["rp", "t0"])
plt.show()

#plot a selection of fits from the chain
for rp, u in samples[np.random.randint(len(samples), size=100)]:
	params.rp = rp
	params.u = [u]
	lc = m.light_curve(params)
	plt.plot(t, lc, color = 'k', alpha = 0.05)
plt.errorbar(t, flux, yerr = err, fmt = 'ok')
plt.xlabel("Time")
plt.ylabel("Relative Flux")
plt.show()
