
        ##################################################################
        #set initial values for MCMC
        ##################################################################
        t0=orbparams[7]
        per=orbparams[5]
        rp=hold_depth[b]
        a=orbparams[0]
        inc=orbparams[1]
        ecc=0.0
        w=90.0
        u=limb
        print(u)
        limb_dark="nonlinear"
        params, m= initialize_model(t,t0,per,rp,a,inc,ecc,w,u,limb_dark)
        ##################################################################
        #             initialize MCMC
        ##################################################################        
        flux=workinglc
        guess_rp, guess_t0, =hold_depth[b], orbparams[7]
        theta=[guess_rp,guess_t0]
        print (theta)
        ndim, nwalkers = len(theta), 50
        sampler=emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(params, m, t, flux, err, up, down))
        pos=[theta+1e-5*np.random.randn(ndim) for i in range(nwalkers)]
        sampler.run_mcmc(pos, 1000);
        
        samples=sampler.chain[:,50:,:].reshape((-1,ndim))
        fig=corner.corner(samples,labels=["rp", "t0"])#, "A/R", "inc"])
        fig.savefig(fpath+'analysis/pyplots/corners/'+str(bincen)+'.png')
        plt.draw()

        #Derive error bars
        rp_mcmc, t0_mcmc=map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))




def lnprior(theta):
    return 0.
    
def lnlike(theta,params,model,t,flux,err, up, down):
    freeparams = theta

    
    residuals =nnbr_res(freeparams, time, lc, err, gw, nbr, orbparams):
    ln_likelihood=-0.5*(np.sum((residuals /err)**2+np.log(2.0*np.pi*(err)**2)))
    
    return ln_likelihood
    
#posterior probability
def lnprob(theta,params,model,t,flux,err, up, down):
    lp=lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp+lnlike(theta,params,model,t,flux,err, up, down)
    