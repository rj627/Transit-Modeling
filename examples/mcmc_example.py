
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
        err=np.ones(len(t))*workinglc**0.5
        plt.figure()
        plt.errorbar(t, workinglc, yerr=err, linestyle='None')
        plt.title(str(bincen))
        plt.draw()
        flux=workinglc
        guess_rp, guess_t0, guess_c0, guess_c1, guess_baseline1, guess_baseline2=hold_depth[b], orbparams[7], 1.0, 0.0, np.max(flux[up]), np.max(flux[down])
        theta=[guess_rp,guess_t0, guess_c0,guess_c1, guess_baseline1, guess_baseline2]
        print(theta)
        ndim, nwalkers = len(theta), 50
        sampler=emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(params, m, t, flux, err, up, down))
        pos=[theta+1e-5*np.random.randn(ndim) for i in range(nwalkers)]
        sampler.run_mcmc(pos, 1000);
        
        samples=sampler.chain[:,50:,:].reshape((-1,ndim))
        fig=corner.corner(samples,labels=["rp", "t0", "c0", "c1", "baseline1", "baseline2"])#, "A/R", "inc"])
        fig.savefig(fpath+'analysis/pyplots/corners/'+str(bincen)+'.png')
        plt.draw()

        #Derive error bars
        rp_mcmc, t0_mcmc,c0_mcmc, c1_mcmc, baseline1_mcmc, baseline2_mcmc=map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))