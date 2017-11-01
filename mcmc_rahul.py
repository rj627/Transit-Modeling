################################################################################
##########                        run_mcmc                 ####################
################################################################################  
    theta=[fit_params[0],fit_params[1]]
    ndim, nwalkers = len(theta), 50
    sampler=emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(time, lc, err, gw, nbr, orbparams))
    pos=[theta+1e-5*np.random.randn(ndim) for i in range(nwalkers)]
    sampler.run_mcmc(pos, 1000);

    samples=sampler.chain[:,50:,:].reshape((-1,ndim))
    fig=corner.corner(samples,labels=[ "t0", "Fp"])#, "A/R", "inc"])
    #fig.savefig(fpath+'analysis/pyplots/corners/'+str(bincen)+'.png')
    plt.show()



    return None


def lnprior(theta, t):

    return 0.


    
def lnlike_nnbr(freeparams,time, lc, err, gw, nbr, orbparams):
 
    res=nnbr_res(freeparams, time, lc, err, gw, nbr, orbparams)
    ln_likelihood=-0.5*(np.sum((res /err)**2+np.log(2.0*np.pi*(err)**2)))

    return ln_likelihood
    
#posterior probability
def lnprob(theta,time, lc, err, gw, nbr, orbparams):
    lp=lnprior(theta, time)
    
    if not np.isfinite(lp):
        return -np.inf
    return lp+lnlike_nnbr(theta, time, lc, err, gw, nbr, orbparams)

    def nnbr_res(coeffs, time, lc, err, gw, nbr, orbparams):
    ecl_model=make_model(coeffs, orbparams, np.squeeze(time), 'secondary')
    lc2=lc/ecl_model

    w1=lc2[nbr]

    w2=np.multiply(w1,gw)
    w3=np.sum(w2, 1)

    w4=np.divide(lc2,w3)
    resids=(w4-1.0)/err

    return resids



def make_model(freeparams, orbparams, ti, prisec):
 
    prisec='secondary'

    params = batman.TransitParams()                #object to store transit parameters
    
    params.per = orbparams[5]                      #orbital period
    params.rp = freeparams[1]                      #planet radius (in units of stellar radii)
    params.a = orbparams[0]                        #semi-major axis (in units of stellar radii)
    params.inc = orbparams[1]                      #orbital inclination (in degrees)
    params.ecc = 0.                                #eccentricity
    params.w = 90.                                 #longitude of periastron (in degrees)
    
    if prisec=='secondary':
        params.fp = freeparams[1]
        params.t_secondary = freeparams[0]   
        params.limb_dark = "uniform"   
        params.u=[]  
        md = batman.TransitModel(params, ti, transittype="secondary")
        flux = md.light_curve(params)    
    else: 
        params.limb_dark = "linear"                    #limb darkening model
        params.u = [0.5]                                 #limb darkening coefficients    
        params.t0 = freeparams[0]                      #time of inferior conjunction                                                                  
        md = batman.TransitModel(params, ti)            #initializes model
        flux = md.light_curve(params)
    model=flux
    #model = (model -  np.amax(model))/-(1.0 - np.amax(model))


    return model