import batman
import numpy as np
import matplotlib.pyplot as plt
from exoplanets_dot_org import *
import limb_darkening as ld
from scipy.optimize import leastsq
from find_nbr import *

def trial():
	#Loads the data#
	data = np.load('WASP_79/62173696.npz')
	lc = data['lc']
	time = data['time']
	holdpos=data['hold_pos']
	xpos=holdpos[:,0]
	ypos=holdpos[:,1]    
	npix=data['beta_np']

	n=250
	lcp = [sum(lc[i:i+n, 0])/n for i in range(0,len(lc[:,0]),n)]
	err=lc**0.5
	orbparams=get_orb_pars('WASP-79 b')
	pred_ecl_time=get_pred_time(orbparams, time, 'transit')
	print (pred_ecl_time)
	freeparams=[pred_ecl_time, orbparams[2]]

	pars, mdl = initialize_model(freeparams, orbparams, np.squeeze(time), 'transit', 'quadratic')
	np.set_printoptions(threshold=np.nan)
	#print (mod)
	mod = mdl.light_curve(pars)
	plt.figure()
	axes = plt.gca()
	axes.set_ylim([0.95,1.05])
	plt.scatter(time[:], lc[:,0]/np.mean(lc[:,0]), s=1)
	plt.scatter(time[:], mod, color='r', s=1)
	plt.show()

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


def get_pred_time(orbparams, t, prisec): 
        num_orbs=np.floor((t[0]+2400000.5-orbparams[6])/orbparams[5])

        if prisec == 'ecl':
        	pred_ecl_time=(orbparams[6]-2400000.5)+(num_orbs+0.5)*orbparams[5]
        else:
        	pred_ecl_time=(orbparams[6]-2400000.5)+num_orbs*orbparams[5]
        	
        if (pred_ecl_time < np.amin(t) or pred_ecl_time > np.amax(t)): 
                raise ValueError('no eclipse in thie time period')
        
        return  pred_ecl_time 

def main():
	trial()

main()