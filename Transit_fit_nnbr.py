import numpy as np
import sys
import pickle
import batman
import matplotlib.pyplot as plt
import os
import glob
import emcee
from scipy.optimize import leastsq
import corner
from find_nbr import *
import pandas as pd
from limb_darkening import find_coeffs


def main():
################################################################################
##########              READ IN THE RAW PHOTOMETRY          ####################
#################################################################################
    numecl=0
    plnm='WASP_79'
    verbose='false'
    fpath='/Users/rahuljayaraman/Documents/Miscellany/Research (Tucker Group)/Python (Transits)/' + plnm
    aorlist=os.listdir(fpath)
    
    #aorlist= [item for item in aorlist if not item.startswith('.')]
    #aorlist=aor_from_list(plnm, 1)
    #aorlist=[50494976]
    aorlist = ['62173184', '62173696']
    #aorlist=np.delete(aorlist, [0,1, len(aorlist)-1])
    for aor in aorlist:
        print (aor) 
        aor=str(aor)
        prisec='primary'
        ramp_style='exp'
        fpathout=fpath+aor+'/apr_fits/'+ramp_style +'/'
        directory = os.path.dirname(fpathout)
        if not os.path.exists(directory):
            os.makedirs(directory)

        #dd=np.load('/Users/Brian/Desktop/Tucker_Group/t_1/outputs/'+plnm+'/'+aor)
        dd=np.load(fpath+'/' + aor + 'extraction.npz') 
        t=dd['time']
        all_lc=dd['lc']
        #hp=dd['hp']
        cp=dd['cp']
        exptime=dd['exptime']
        framtime=0.1
        orbparams=dd['op']
        holdpos=dd['hold_pos']
        npix=dd['beta_np']
        chnum=dd['ch']
        red_all=[]


    ################################################################################
        pred_ecl_time=get_pred_time(orbparams, t, prisec)
        print (pred_ecl_time - t[0])
        
        freeparams=[pred_ecl_time-t[0], orbparams[2]]
    
        if prisec == 'secondary': 
            freeparams[1]=0.0011
            ldc=[]
        else: ldc=find_coeffs(orbparams[10], orbparams[9], orbparams[8], 2, 'quadratic')#(temp, log_g, metallicity, channel, type_limb)

        for apr in range(0,all_lc.shape[1]):
            
            directory = os.path.dirname(fpathout)
            if not os.path.exists(directory):
                os.makedirs(directory)
            lc=np.squeeze(all_lc[:,apr]*2.35481)
            time=(t-t[0])
            time=np.squeeze(time)
            norm=np.nanmedian(lc)
            #print('Photon Noise limit is: ',(np.sqrt(norm*1.002)/(norm*1.002)))
         
            err=1.1*lc**0.5
            lc=lc/norm
            err=err/norm
            err=np.ones(len(lc))*0.0045
       
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
            #fpathout='/Users/Brian/Desktop/Tucker_Group/Spitzer/mapping_files/outputs/'+plnm+'/'+aor+'/apr_fits/'
            filt_file = fpathout+'post_filter_'+str(apr)+'.npz'

            #print(filt_file)
            if os.path.isfile(filt_file):
                if verbose=='true': print('Found Filter File')
                ff=np.load(filt_file)
                lc=ff['lc']
                #cp3=ff['cp3']
                time=ff['time']
                xpos=ff['xpos']
                ypos=ff['ypos']
                npix=ff['npix']
                err=ff['err']
                found='true'

            else: 
                found='false'
                if verbose=='true': print('In Filter')  
                lc, cp3, time, xpos, ypos, npix, err=filter_data(lc, cp3, time, xpos, ypos, npix, dep_ind, err)
                if verbose=='true': print('Out of Filter')             
            
            plt.figure()
            plt.title(plnm+' Ch: '+str(chnum)+'\n'+str(aor)+'_'+str(apr))
            plt.axvline(x=pred_ecl_time-t[0])
            plt.axvline(x=pred_ecl_time-orbparams[4]*0.5-t[0],color = 'r', linestyle = 'dashed')
            plt.axvline(x=pred_ecl_time+orbparams[4]*0.5-t[0],color = 'r', linestyle = 'dashed')
            plt.scatter(time, lc, s=1)
            if prisec=='secondary': plt.ylim(0.95,1.05)
            else: plt.ylim(0.95, 1.03)

            #plt.xlim(time[0], np.amax(time))
            plt.savefig(fpathout+'raw_lc_plot_'+str(apr))
            if verbose=='true':
                plt.draw()
                plt.pause(1200)
            plt.close('all')

            # time2=np.multiply(time, time)
            # time=time[np.newaxis]
            # time2=time2[np.newaxis]
            # t2hours=time2*24.0**2.0
            # thours=time*24.0
            
    ################################################################################
    ##########                  TRIM THE DATA                 ####################
    ################################################################################     
            trim_time=0.  #in minutes
            if trim_time != 0.:
                trim_time=trim_time/(60.*24.0)  #convert to days
                start_index=int(trim_time/(exptime/86400.0))
                end_ind=np.squeeze(lc)
                end_ind=end_ind.size
            
                print (exptime)

                lc=lc[start_index:end_ind]
                time=np.squeeze(time[start_index:end_ind])
                xpos=xpos[start_index:end_ind]
                ypos=ypos[start_index:end_ind]
                npix=npix[start_index:end_ind]
                err=err[start_index:end_ind]
                plt.figure()
                plt.scatter(time, lc, s=1)
                plt.draw()
################################################################################
##########             FIND NEIGHBORS                ####################
################################################################################
            
            if found=='true':
                gw=ff['gw']
                nbr=ff['nbr']
            else:
                if verbose=='true':  print('In Find NBR')
                gw, nbr=find_nbr_qhull(xpos, ypos, npix, sm_num = 50, a = 1.0, b = 1.7777, c = 1.0, print_space = 10000.)
                if verbose=='true':  print('Out of Find NBR')
            np.savez(fpathout+'post_filter_'+str(apr), lc=lc, cp3=cp3, time=time, xpos=xpos, ypos=ypos, npix=npix, err=err, gw=gw, nbr=nbr, orbparams=orbparams, pred_ecl_time=pred_ecl_time)
################################################################################
##########                  FIT THE DATA                 ####################
################################################################################     

            if prisec=='secondary': freeparams=[pred_ecl_time-t[0], orbparams[2], 0.005, 0.05]  #the last 2 free params are ramp terms
            else:  
                if ramp_style=='linear': freeparams=[pred_ecl_time-t[0], orbparams[2], 0.00001, 1.000001] 
                if ramp_style=='exp': freeparams=[pred_ecl_time-t[0], orbparams[2], 0.005, 0.05] 
                if ramp_style=='none': freeparams=[pred_ecl_time-t[0], orbparams[2],1.0,1.0]
            params,m=initialize_model(np.squeeze(time), freeparams, orbparams, prisec, ldc)
            fluxcurve = m.light_curve(params)
            fit_params, pcov, infodict,flag, sucess=leastsq(nnbr_res, freeparams, args=(time, lc, err, gw, nbr, params, m, prisec, ramp_style), full_output=1)
            print('apr# '+str(apr), fit_params)
            file_name=fpathout+'apr_fit_'+str(apr)
            fileObject = open(file_name,'wb') 
            pickle.dump([lc, time, err, gw, nbr, fit_params],fileObject)
            fileObject.close()

################################################################################
##########                  PLOT THE FIT                ####################
################################################################################   
            if prisec=='secondary': 
                params.t_secondary=fit_params[0]
                params.fp=fit_params[1]
            else:
                params.t0=fit_params[0]
                params.rp=fit_params[1]
            eclipse_model=m.light_curve(params)
            ramp=ramp_model([fit_params[2], fit_params[3]], time, ramp_style)
            lc2=np.squeeze(lc/eclipse_model/ramp)

            w1=lc2[nbr]
            w2=np.multiply(w1,gw)
            w3=np.sum(w2, 1)
            w4=np.divide(lc2,w3)
            w5=w4*eclipse_model
            resids=(w4-1.)#/err
            res2=(lc/eclipse_model-1.0)/err

            pltbins=64

            blc=bin_anything(w5, pltbins)
            btime=bin_anything(time, pltbins)

            if prisec=='secondary': phase=0.5+(time+t[0]-pred_ecl_time)/orbparams[5]
            if prisec=='primary': phase=0.0+(time+t[0]-pred_ecl_time)/orbparams[5] 
            bphase=bin_anything(phase, pltbins)


            plt.figure()
            plt.title(plnm+' Ch: '+str(chnum)+'\n'+str(aor)+'_'+str(apr))
            plt.scatter(bphase, blc,s=10)
            #plt.scatter(time, lc, alpha=0.1, color='b', s=1)
            plt.plot(np.squeeze(phase),eclipse_model, color='r')
            if prisec=='secondary': 
                plt.ylim(0.9975, 1.0035)
                plt.text(0.47, 1.003, 'T_center O-C (s): '+str(round((fit_params[0]+t[0]-pred_ecl_time)*86400.,1))+'                   Depth: '+str(round(fit_params[1]*1.0e6, 0))+' ppm')
                plt.text(0.49, 1.0025, 'SDNR:  '+str(round(np.std(resids), 6)))
            else:
                plt.ylim(0.983, 1.005)
                plt.text(0.43, 0.9925, 'T_center O-C (s): '+str(round((fit_params[0]+t[0]-pred_ecl_time)*86400.,1)))
                plt.text(0.43, 0.990,'Transit Depth: '+str(round(fit_params[1]**2.*100, 4))+' %')
                plt.text(0.43, 0.9875, 'SDNR:  '+str(round(np.std(resids), 6)))
            plt.xlabel('Phase Units')
            plt.ylabel('Relative Flux')
            
            plt.savefig(fpathout+'apr_fit_plot_'+str(apr))
            if verbose=='true':
                plt.draw()
                plt.pause(1.2)
           
          
            
################################################################################
##########                 Get Red Noise                    ####################
################################################################################   
            sdnr, beta_red=est_rednoise(resids, framtime, fpathout, aor, apr, plnm, chnum, prisec)
            if red_all ==[]: red_all=np.ones(shape=(all_lc.shape[1], 5))*1000.
            red_all[apr,:]=[sdnr, beta_red*sdnr, beta_red, round(fit_params[1]*1.e6, 1), fit_params[0]]
          
             



        best=np.nanargmin(red_all,axis=0)
        best=best[1]


        np.save(fpathout+aor+'_summary', red_all)
        np.savetxt(fpathout+aor+'_summary', red_all)
        if verbose=='true': print(best)




################################################################################
##########                 Load the best apr results        ####################
################################################################################  

        filename=fpathout+'apr_fit_'+str(best)
        fileObject = open(filename,'rb') 
        lc, time, err, gw, nbr, fit_params=pickle.load(fileObject)
        err=err*red_all[best, 2]

        print('Best Beta_red', red_all[best,2])
        params,m=initialize_model(np.squeeze(time), freeparams, orbparams, prisec, ldc)

################################################################################
##########                        run_mcmc                 ####################
################################################################################  
        theta=fit_params
        ndim, nwalkers = len(theta), 20
        sampler=emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(time, lc, err, gw, nbr, params, m, prisec, ramp_style))
        pos=[theta+1.e-4*np.random.randn(ndim) for i in range(nwalkers)]
        sampler.run_mcmc(pos, 1500);

        samples=sampler.chain[:,50:,:].reshape((-1,ndim))
        np.save(fpathout+aor+'_samples', samples)
        if prisec == 'primary':
            fig=corner.corner(samples,labels=[ "t0", "rp", "a1", "a2"])#, "A/R", "inc"])
        else: fig=corner.corner(samples,labels=[ "t0", "Fp", "a1", "a2"])#, "A/R", "inc"])
        fig.savefig(fpathout+aor+'_corner_'+str(best)+'.png')
        #plt.show(block=False)
        #plt.pause(0.5)

 #Derive error bars
        t0_mcmc, rp_mcmc, a1_mcmc, a2_mcmc=map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
        print(rp_mcmc, t0_mcmc)
        np.savez(fpathout+aor+'_mcmc_results', rp_mcmc=rp_mcmc, t0_mcmc=t0_mcmc, a1_mcmc=a1_mcmc, a2_mcmc=a2_mcmc, best=best)
        phase=0.0+(time+t[0]-pred_ecl_time)/orbparams[5] 
        bphase=bin_anything(phase, pltbins)

        plt.figure()      
        for t0, rp, a1, a2 in samples[np.random.randint(len(samples), size=100)]:
            params.rp = rp
            params.t0 = t0
            ecl_mod = m.light_curve(params)

            plt.plot(phase, ecl_mod, color = 'k', alpha = 0.05)

            ramp=ramp_model([a1, a2], time, ramp_style)
            lc2=np.squeeze(lc/ecl_mod/ramp)
            w1=lc2[nbr]
            w2=np.multiply(w1,gw)
            w3=np.sum(w2, 1)
            w4=np.divide(lc2,w3)
            w5=w4*ecl_mod
            resids=(w4-1.)#/err
            res2=(lc/ecl_mod-1.0)/err

            blc=bin_anything(w5, pltbins)
            btime=bin_anything(time, pltbins)
        plt.scatter(bphase, blc, s=8, alpha=0.5)
        plt.xlabel("Phase Units")
        plt.ylabel("Relative Flux")
        plt.title(plnm+' Ch: '+ str(chnum))
        plt.show()
        #plt.savefig('/Users/Brian/Desktop/W79_summary/'+str(chnum)+'_mcmc_fit')

    return None


################################################################################
##########                  END MAIN                        ####################
################################################################################ 
def nnbr_res(coeffs, time, lc, err, gw, nbr, params, m, prisec, ramp_style):
    
    time=np.squeeze(time)
    if ramp_style == 'none': 
        ramp_coeffs=[-1.,-1.]
    else:
        ramp_coeffs=[coeffs[2], coeffs[3]]

    if prisec=='secondary':
        params.t_secondary=coeffs[0]
        params.fp=coeffs[1]

    else:  
        params.t0=coeffs[0]
        params.rp=coeffs[1]

    ecl_model=m.light_curve(params)
    
    ramp = ramp_model(ramp_coeffs, time, ramp_style)
    lc2=np.squeeze(lc/(ecl_model*ramp))
    w1=lc2[nbr]

    w2=np.multiply(w1,gw)
    w3=np.sum(w2, 1)

    w4=np.divide(lc2,w3)
    resids=(w4-1.0)/err
    
    return resids

def ramp_model(ramp_coeffs, time, ramp_style):

    if ramp_style != 'none':
        a1,a2=ramp_coeffs[0], ramp_coeffs[1]
  
    if ramp_style=='exp':  
        ramp=1.-a1*np.exp(-1.0*time/a2)
    elif ramp_style=='linear': ramp= a1*time+a2
    elif ramp_style=='none':  ramp=np.ones(len(time))

    return ramp

def initialize_model(t, freeparams, orbparams, prisec, ldc):
    params = batman.TransitParams()
    params.t0 = freeparams[0]              
    params.per = orbparams[5]        
    params.rp = orbparams[2]      
    params.a = orbparams[0]                
    params.inc = orbparams[1]            
    params.ecc = 0.        
    params.w = 90.  
    if prisec=='primary':
        params.limb_dark = "quadratic"
        params.u=ldc 
        model = batman.TransitModel(params, t) 
    else:            
        params.u = []                       
        params.limb_dark = 'uniform'  
        params.fp = 0.001
        params.t_secondary = freeparams[0]       
        model = batman.TransitModel(params, t, transittype="secondary") 

    return params, model            #return parameters and model objects 


    
def filter_data(lc, cp, t, xpos, ypos, npix, dep_ind, err): 
    #print('In Filter')
    lch=lc
    th=t
    plt.close('all')
    for c in range(0,3):
        sigma=np.std(lc)
        bad_data=[i for i,v in enumerate(lc) if np.absolute(v-np.median(lc)) > 3.0*(sigma)] 
        lc=np.delete(lc, bad_data)
        t=np.delete(t, bad_data)
        cp=np.delete(cp, bad_data, axis=1)
        xpos=np.delete(xpos, bad_data)
        ypos=np.delete(ypos, bad_data)
        npix=np.delete(npix, bad_data)
        err=np.delete(err, bad_data)
    for c in range(0,1):
        sigma=np.std(xpos)
        bad_data=[i for i,v in enumerate(xpos) if np.absolute(v-np.median(xpos)) > 5.0*(sigma)] 
        lc=np.delete(lc, bad_data)
        t=np.delete(t, bad_data)
        cp=np.delete(cp, bad_data, axis=1)
        xpos=np.delete(xpos, bad_data)
        ypos=np.delete(ypos, bad_data)
        npix=np.delete(npix, bad_data)
        err=np.delete(err, bad_data)
    for c in range(0,1):
        sigma=np.std(ypos)
        bad_data=[i for i,v in enumerate(ypos) if np.absolute(v-np.median(ypos)) > 5.0*(sigma)] 
        lc=np.delete(lc, bad_data)
        t=np.delete(t, bad_data)
        cp=np.delete(cp, bad_data, axis=1)
        xpos=np.delete(xpos, bad_data)
        ypos=np.delete(ypos, bad_data)
        npix=np.delete(npix, bad_data)
        err=np.delete(err, bad_data)

    percent_trimmed=float((lch.size-lc.size))/float(lch.size)*100.0
    #print('Percent trimmed:  ', percent_trimmed)
    #print('Out of Filter')
    return lc, cp,t, xpos, ypos, npix, err


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

def est_rednoise(resids, exptime, fpathout, aor, apr, plnm, chnum, prisec):
    binsize=np.ones(5)+1
    bpow=np.arange(5)
    binsize=np.power(binsize, bpow)
    
    b2=np.arange(16, 2048, 16)
    binsize=np.append(binsize, b2)

    if exptime <1.: binsize*=10.
    bin_std=np.zeros(binsize.size)
    bin_time=np.zeros(binsize.size)
    beta_red=np.zeros(binsize.size)
    numbins=binsize.size
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
        beta_red[b]=bin_std[b]/bin_std[0]*np.sqrt(binsize[b]*(numbins-1.)/numbins)
    y=np.log10(bin_std[0])-0.5*(np.log10(bin_time)-np.log10(bin_time[0]))
    broadb=np.sum(np.abs(np.log10(bin_std)-y)**0.5)

    tscales=np.where((np.log10(bin_time)>2.5))# & (np.log10(bin_time)<3.0))
    if prisec=='primary': tscales=np.where((np.log10(bin_time)>2.0) & (np.log10(bin_time)<2.5))#
    bred=np.nanmedian(beta_red[tscales])

    
    plt.figure()
    plt.title(plnm+' Ch: '+str(chnum)+'\n'+str(aor)+'_'+str(apr))
    plt.scatter(np.log10(bin_time), np.log10(bin_std), color='red')
    plt.plot(np.log10(bin_time),y ,linewidth=3.3)
    plt.text(1.5,-2.5, 'Beta_red:  '+str(round(np.amax(beta_red), 2)))
    plt.xlabel('Bin Size (log (s))')
    plt.ylabel('Log Std Dev of Residuals')
    #plt.savefig(fpathout+aor+'_red_noise_'+str(apr))

    return  bin_std[0], bred

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def get_pred_time(orbparams, t, prisec):
        
        num_orbs=np.floor((t[0]+2400000.5-orbparams[6])/orbparams[5])
        if prisec=='secondary':  pred_ecl_time=(orbparams[6]-2400000.5)+(num_orbs+0.5)*orbparams[5]
        else: pred_ecl_time=(orbparams[6]-2400000.5)+(num_orbs+1.)*orbparams[5]
        
        if pred_ecl_time < np.amin(t): 
            pred_ecl_time=pred_ecl_time+orbparams[5]
            if pred_ecl_time > np.amax(t): 
                raise ValueError('no event in this time period')
        
        return  pred_ecl_time  


    
def lnprior(theta, t):
    if theta[3] < 0.: return -np.inf

    return 0.



    
def lnlike_nnbr(freeparams,time, lc, err, gw, nbr,  params, m, prisec, ramp_style):
 
    res=nnbr_res(freeparams, time, lc, err, gw, nbr, params, m, prisec, ramp_style)
   
    ln_likelihood=-0.5*(np.sum((res)**2+np.log(2.0*np.pi*(err)**2)))
    
    #ln_likelihood=-0.5*(np.sum((res)**2+np.log(2.0*np.pi*(err)**2)))  #Here I removed err from denominator since nnbr res returns weighted res

    return ln_likelihood
    
#posterior probability
def lnprob(theta,time, lc, err, gw, nbr, params, m, prisec, ramp_style):
    lp=lnprior(theta, time)
    
    if not np.isfinite(lp):
        return -np.inf
    return lp+lnlike_nnbr(theta, time, lc, err, gw, nbr, params, m, prisec, ramp_style)

def aor_from_list(planet, ch):
    aor_arr=pd.read_csv('/Users/Brian/Desktop/Tucker_Group/Spitzer/mapping_files/aor_list.csv')
    chlist=aor_arr.CH
    aorlist=aor_arr.AOR
    namelist=aor_arr.Target
    ch_aors=np.where((namelist==planet) & (chlist==ch))
    aor_selected=np.array(aorlist[np.squeeze(ch_aors)])


    return aor_selected


main()
