import numpy as np
import sys
import batman
import matplotlib.pyplot as plt
import os
from find_nbr import *
import glob
import emcee
from scipy.optimize import leastsq
import corner
import limb_darkening as ld
from exoplanets_dot_org import *

def main():
    sys.settrace
################################################################################
##########              READ IN THE RAW PHOTOMETRY          ####################
#################################################################################
    numecl=0
    plnm='WASP_79'

    aorlist=os.listdir('/Users/rahuljayaraman/Documents/Miscellany/Research (Tucker Group)/Python (Transits)/'+plnm)

    for aor in aorlist:
        aor=aorlist[1]
        prisec='primary'
        
        dd=np.load('/Users/rahuljayaraman/Documents/Miscellany/Research (Tucker Group)/Python (Transits)/'+plnm+'/'+aor + '.npz')
        
        t=dd['time']
        all_lc=dd['lc']
        hp=dd['hp']
        cp=dd['cp']
        exptime=dd['exptime']
        orbparams=get_orb_pars('WASP-79 b')
        holdpos=dd['hold_pos']

        #plt.figure()
        #axes = plt.gca()
        #axes.set_ylim([0.95,1.05])
        #plt.scatter(t[:], all_lc[:,1]/np.mean(all_lc[:,1]), s=1)
        #plt.show()

    ################################################################################
        pred_ecl_time=get_pred_time(orbparams, t, 'transit')
        
        freeparams=[pred_ecl_time, orbparams[11]]
        #print(freeparams)
        #if prisec == 'secondary': freeparams[1]=0.0011
        print (all_lc.shape[1])
      
        for apr in range(0,all_lc.shape[1]):
            lc=np.squeeze(all_lc[:,apr])
            time=(t-t[0])
            time=np.squeeze(time)
            norm=np.nanmedian(lc)
            err=lc**0.5
            lc=lc/norm
            err=err/norm
            xpos=holdpos[:,0]
            ypos=holdpos[:,1]    
            npix=dd['beta_np']

################################################################################
##########              NORMALIZE THE PIXEL VALUES          ####################
################################################################################
            timelength=len(t)
            #cp1=cp[1:4, 1:4, :]   
            cp1=cp
            dep_ind=cp1.shape[0]*cp1.shape[1]
            cp2=np.reshape(cp1, (dep_ind,timelength))
            cp3=cp2#[:,start:end]      
            for p in range(0,len(time)):
                norm=np.sum(cp3[:, p])
                cp3[:, p]/=norm
    ################################################################################
    ##########                  FILTER THE DATA                 ####################
    ################################################################################        
            lc, cp3, time, xpos, ypos, npix, err = filter_data(lc, cp3, time, xpos, ypos, npix, dep_ind, err)
            print (len(t))
           
            
            # plt.figure()
            # plt.axvline(x=pred_ecl_time-t[0])
            # plt.axvline(x=pred_ecl_time-orbparams[4]*0.5-t[0],color = 'r', linestyle = 'dashed')
            # plt.axvline(x=pred_ecl_time+orbparams[4]*0.5-t[0],color = 'r', linestyle = 'dashed')
            # plt.scatter(time, lc, s=1)
            # #plt.ylim(0.95,1.05)
            # plt.xlim(time[0], np.amax(time))
            # plt.title(aor+str(apr))
            # plt.draw()
            # plt.pause(0.2)
            # plt.close('all')

            time2=np.multiply(time, time)
            time=time[np.newaxis]
            time2=time2[np.newaxis]
            t2hours=time2*24.0**2.0
            thours=time*24.0
            #print( 'Out of Filter')

    ################################################################################
    ##########                  TRIM THE DATA                 ####################
    ################################################################################     
            trim_time=10.  #in minutes
            trim_time=trim_time/(60.*24.0)  #convert to days
            start_index=int(trim_time/(exptime/86400.0))
            end_ind=np.squeeze(lc)
            end_ind=end_ind.size
            
            

            lc=lc[start_index:end_ind]
            time=np.squeeze(time[0,start_index:end_ind])
            xpos=xpos[start_index:end_ind]
            ypos=ypos[start_index:end_ind]
            npix=npix[start_index:end_ind]
            err=err[start_index:end_ind]

            # plt.figure()
            # plt.axvline(x=pred_ecl_time-t[0])
            # plt.axvline(x=pred_ecl_time-orbparams[4]*0.5-t[0],color = 'r', linestyle = 'dashed')
            # plt.axvline(x=pred_ecl_time+orbparams[4]*0.5-t[0],color = 'r', linestyle = 'dashed')
            # plt.scatter(time, lc, s=1)
            # plt.ylim(0.95,1.02)
            # plt.xlim(time[0], np.amax(time))
            # plt.title(aor+str(apr))
            # plt.draw()
            # plt.pause(1.0)
            # plt.close('all')



################################################################################
##########             FIND NEIGHBORS                ####################
################################################################################

            gw, nbr=find_nbr_qhull(xpos, ypos, npix, sm_num = 50, a = 1.0, b = 1.0, c = 1.0, print_space = 10000.)
################################################################################
##########                  FIT THE DATA                 ####################
################################################################################   
            params, model =initialize_model(freeparams, orbparams, np.squeeze(time), prisec, 'quadratic')
            fit_params, pcov, infodict,flag, sucess=leastsq(nnbr_res, freeparams, args=(time, lc, err, gw, nbr, orbparams, params, model), full_output=1)

################################################################################
##########                  PLOT THE FIT                ####################
################################################################################   
            
            fluxcurve = model.light_curve(params)
            lc2=lc/fluxcurve

            w1=lc2[nbr]
            w2=np.multiply(w1,gw)
            w3=np.sum(w2, 1)
            w4=np.divide(lc2,w3)
            w5=w4*fluxcurve
            resids=(w4-1.)#/err
            res2=(lc/fluxcurve-1.0)/err
            blc=bin_anything(w5, 64)
            btime=bin_anything(time, 64)
            
            plt.figure()
            plt.scatter(btime, blc,s=5)
            plt.scatter(time, lc, alpha=0.1, color='b', s=1)
            plt.scatter(time,fluxcurve, color='r', s=1)
            plt.title(aor+str(apr))
            plt.ylim(0.98, 1.02)

            plt.show()

################################################################################
##########                 Get Red Noise              ####################
################################################################################   
            bred, sdnr=est_rednoise(resids, exptime)
            if apr ==0: red_all=np.zeros(shape=(all_lc.shape[1], 4))
            red_all[apr,:]=[bred, sdnr, bred*sdnr, fit_params[1]*1.e6]

            np.savez("apr_" + str(apr) + "_info", nbr=nbr, gw=gw, lc=lc, err=err, ti=time, fp=fit_params, op=orbparams, t=t)


        best_index = np.argmin(red_all[:,2])
        best_arr = red_all[best_index]
        data = np.load('apr_' + str(best_index) + "_info.npz")

        orbparams  = data['op']
        fit_params = data['fp']
        time = data['ti']
        t = data['t']
        flux = data['lc']
        gw = data['gw']
        err = data['err']
        nbr = data['nbr']
        pars, mdl = initialize_model(fit_params, orbparams, np.squeeze(time), 'transit', 'quadratic')
    
        ##################################################################
        #             initialize MCMC
        ##################################################################        
        #flux=mdl.light_curve(pars) #????#
        # plt.figure()
        # axes = plt.gca()
        # axes.set_ylim([0.95,1.05])
        # plt.scatter(time[:], lc[:]/np.mean(lc[:]), s=1)
        # plt.scatter(time[:], em, color='r', s=1)
        # plt.show()
        # #fit_params
        #guess_rp, guess_t0, = orbparams[11], orbparams[12]
        guess_rp, guess_t0 = 0.1, 0.2                        
        theta = [guess_t0, guess_rp]
        print ('post theta')
        ndim, nwalkers = len(theta), 50
        sampler=emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(time, flux, err, gw, nbr, orbparams, pars, mdl))
        print ('created sampler')
        pos=[theta+1e-5*np.random.randn(ndim) for i in range(nwalkers)]
        sampler.run_mcmc(pos, 1000)
        print ('post mcmc')
        samples=sampler.chain[:,50:,:].reshape((-1,ndim))
        fig=corner.corner(samples,labels=["rp", "t0"])#, "A/R", "inc"])
        #fig.savefig(fpath+'analysis/pyplots/corners/'+str(bincen)+'.png')
        plt.draw()

        for rp, t0 in samples[np.random.randint(len(samples), size=100)]:
            pars.rp = rp
            pars.t0 = t0
            lc = mdl.light_curve(pars)
            plt.plot(time, lc, color = 'k', alpha = 0.05)
        #plt.errorbar(t, b, yerr = err, fmt = 'ok')
        plt.xlabel("Time")
        plt.ylabel("Relative Flux")
        plt.show()


        #Derive error bars
        rp_mcmc, t0_mcmc=map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
      
    #Lsq minimizer to fit the model - example 


    #find minimum in red noise array based on bred*sdnr
    #save gw, nbr, lc, time, orbparams etc.to file, 1 per aperture


    print(red_all)
    best=np.argmin(red_all,axis=0)
    best=best[2]
    print(best)

    return None


################################################################################
##########                  END MAIN                        ####################
################################################################################ 


def initialize_model(freeparams, orbparams, t, prisec, limbd):
    params = batman.TransitParams()
    params.t0 = freeparams[0]            
    params.per = orbparams[5]        
    params.rp = freeparams[1]      
    params.a = orbparams[0]              
    params.inc = orbparams[1]            
    params.ecc = 0.        
    params.w = 90.              
    params.u = ld.find_coeffs(orbparams[10], orbparams[8], orbparams[9], 2, limbd) 
    if "nonlinear" in limbd :                   
        params.limb_dark = "nonlinear"
    else:
        params.limb_dark = limbd        

    print (params.u)
    print (params.limb_dark)
    model = batman.TransitModel(params, t)  

    return params, model          #return parameters and model objects 

# def make_model(freeparams, orbparams, ti, prisec, limbd):

#     params = batman.TransitParams()                #object to store transit parameters

#     params.per = orbparams[5]                      #orbital period
#     params.rp = freeparams[1]                      #planet radius (in units of stellar radii)
#     params.a = orbparams[0]                        #semi-major axis (in units of stellar radii)
#     params.inc = orbparams[1]                      #orbital inclination (in degrees)
#     params.ecc = 0.                                #eccentricity
#     params.w = 90.                                 #longitude of periastron (in degrees)

#     if prisec=='secondary':
#         params.fp = freeparams[1]
#         params.t_secondary = freeparams[0]   
#         params.limb_dark = limbd  
#         params.u=[]  
#         md = batman.TransitModel(params, ti, transittype="secondary")    #insert limb darkening coeffs here, please! :) 
#         flux = md.light_curve(params)    
#     else: 
#         params.limb_dark = limbd                  #limb darkening model
#         params.u = ld.find_coeffs(orbparams[10], orbparams[8], orbparams[9], 2, 'q')                             #limb darkening coefficients    
#         params.t0 = freeparams[0]                      #time of inferior conjunction                                                                  
#         md = batman.TransitModel(params, ti)            #initializes model
#         flux = md.light_curve(params)
#     ret_flux = flux
#     #model = (model -  np.amax(model))/-(1.0 - np.amax(model))


#     return params, md, ret_flux

    
def filter_data(lc, cp, time, xpos, ypos, npix, dep_ind, err): 
    #print('In Filter')
    lch=lc
    th=time
    plt.close('all')
    for c in range(0,3):
        sigma=np.std(lc)
        bad_data=[i for i,v in enumerate(lc) if np.absolute(v-np.median(lc)) > 5.0*(sigma)] 
        lc=np.delete(lc, bad_data)
        time=np.delete(time, bad_data)
        cp=np.delete(cp, bad_data, axis=1)
        xpos=np.delete(xpos, bad_data)
        ypos=np.delete(ypos, bad_data)
        npix=np.delete(npix, bad_data)
        err=np.delete(err, bad_data)

    # for p in range(0,dep_ind):
    #    ts=cp[p,:]
    #    sigma=np.std(ts)
    #    bad_data=[i for i,v in enumerate(ts) if np.absolute(v-np.median(ts)) > 5.0*(sigma)] 
    #    lc=np.delete(lc, bad_data)
    #    t=np.delete(t, bad_data)
    #    cp=np.delete(cp, bad_data, axis=1)
    

    percent_trimmed=float((lch.size-lc.size))/float(lch.size)*100.0
    print('Percent trimmed:  ', percent_trimmed)


    return lc, cp, time, xpos, ypos, npix, err
def bin_anything(series, binsize):
    series=np.squeeze(series)
    numbins=int(np.floor(series.size/binsize))
    bin_series=np.zeros(numbins)
    for b in range(0,numbins):
        bin_series[b]=np.mean(series[b*binsize:(b+1)*binsize-1])      

    return bin_series
        
def bin_data(time, lc, cp, binsize, dep_ind):
    time=np.squeeze(time)
    numbins=int(np.floor(lc.size/binsize))
 
    blc=np.zeros(numbins)
    bt=np.zeros(numbins)
    bcp=np.zeros(shape=(dep_ind,numbins))
  
    for b in range(0,numbins):
        blc[b]=np.mean(lc[b*binsize:(b+1)*binsize-1])      
        bt[b]=np.mean(time[b*binsize:(b+1)*binsize-1])
    # for p in range(0,dep_ind):
    #     for bb in range(0,numbins):
    #         bcp[p, bb]=np.nanmean(cp[p, bb*binsize:(bb+1)*binsize-1])

    return bt, blc, bcp

def est_rednoise(resids, exptime):
    binsize=np.ones(5)+1
    bpow=np.arange(5)
    binsize=np.power(binsize, bpow)
    
    b2=np.arange(16, 2048, 16)
    binsize=np.append(binsize, b2)

    if exptime <1.: binsize*=4
    bin_std=np.zeros(binsize.size)
    bin_time=np.zeros(binsize.size)
    for b in range(0,binsize.size):
        if b ==0 : 
            bin_std[b]=np.std(resids)
            bin_time[b]=exptime
        else:
            bb=binsize[b]
            numbins=int(np.floor(resids.size/bb)) 
            bres=np.zeros(numbins)    
            for v in range(0,numbins):
                bres[v]=np.mean(resids[int(v*bb):int((v+1)*bb-1)]) 
            bin_std[b]=np.std(bres)
            bin_time[b]=exptime*bb
    y=np.log10(bin_std[0])-0.5*(np.log10(bin_time)-np.log10(bin_time[0]))
    broadb=np.sum(np.abs(np.log10(bin_std)-y)**0.5)
    # plt.figure()
    # plt.scatter(np.log10(bin_time), np.log10(bin_std), color='red')
    # plt.plot(np.log10(bin_time),y ,linewidth=3.3)
    # #plt.text(2,-2.1, 'Binned at:  '+str(np.log10(exptime*fitbin)))
    # plt.draw()
    # plt.pause(1.0)

    #plt.pause(0.2)
    
    return broadb, bin_std[0]#*bin_std[0]

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def get_pred_time(orbparams, t, prisec): 
    num_orbs=np.floor((t[0]+2400000.5-orbparams[6])/orbparams[5])

    if prisec == 'ecl':
        pred_ecl_time=(orbparams[6]-2400000.5)+(num_orbs+0.5)*orbparams[5]
    else:
        pred_ecl_time=(orbparams[6]-2400000.5)+num_orbs*orbparams[5]
            
    if (pred_ecl_time < np.amin(t) or pred_ecl_time > np.amax(t)): 
            raise ValueError('no eclipse in thie time period')
        
    return  pred_ecl_time

    
#prior
def lnprior(theta):
    return 0.                           #assumes all priors have uniform probability    

#likelihood function 
def lnlike(theta, params, model, t, flux, err):
    params.t0, params.rp  = theta[0], theta[1]          #update parameters
    lc = model.light_curve(params)
    residuals = flux - lc
    ln_likelihood = -0.5*(np.sum((residuals/err)**2 + np.log(2.0*np.pi*(err)**2)))

    return ln_likelihood

#posterior probability
def lnprob(theta,time, lc, err, gw, nbr, orbparams, params, m):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike_nnbr(theta, time, lc, err, gw, nbr, orbparams, params, m)

def lnlike_nnbr(freeparams,time, lc, err, gw, nbr, orbparams, params, m):
 
    res=nnbr_res(freeparams, time, lc, err, gw, nbr, orbparams, params, m)
    ln_likelihood=-0.5*(np.sum((res /err)**2+np.log(2.0*np.pi*(err)**2)))

    return ln_likelihood

def nnbr_res(coeffs, time, lc, err, gw, nbr, orbparams, params, m):
    params.t0 = coeffs[0]
    params.rp = coeffs[1]

    flx=m.light_curve(params)#update the 2 free params, generate a new lightcurve
    print (flx)
    lc2=lc/flx

    w1=lc2[nbr]

    w2=np.multiply(w1,gw)
    w3=np.sum(w2, 1)

    w4=np.divide(lc2,w3)
    resids=(w4-1.0)/err

    return resids

main()
