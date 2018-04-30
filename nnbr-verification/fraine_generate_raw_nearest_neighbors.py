# coding: utf-8
import pyhull

from functools       import partial
from multiprocessing import Pool, cpu_count
from pylab           import *;ion()
from time            import time
from tqdm            import tqdm

data_filename = 'xpos_ypos_npix_example.dat'
try:
    xpos, ypos, npix, flux = np.loadtxt(data_filename).T
    nPts = xpos.size
    
    print('loaded {} from system'.format(data_filename))
except:
    nPts = int(1e5)
    
    xpos = 0.35*sin(np.arange(0,nPts) / 1500 + 0.5) + 15 + np.random.normal(0,0.2,nPts)
    ypos = 0.35*sin(np.arange(0,nPts) / 2000 + 0.7) + 15 + np.random.normal(0,0.2,nPts)
    npix = 0.25*sin(np.arange(0,nPts) / 2500 + 0.4) + 15 + np.random.normal(0,0.2,nPts)
    flux = 1+0.01*(xpos - xpos.mean()) + 0.01*(ypos - ypos.mean()) + 0.01*(npix - npix.mean())
    
    np.savetxt(data_filename, np.transpose([xpos, ypos, npix, flux]))

plot_raw = False
if plot_raw:
    fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1)

    ax1.plot(xpos, '.', label='X-Position')
    ax2.plot(ypos, '.', label='Y-Position')
    ax3.plot(npix, '.', label='Noise Pixels')
    ax4.plot(flux, '.', label='Flux')

    ax1.set_ylabel('X-Positions')
    ax2.set_ylabel('Y-Position')
    ax3.set_ylabel('Noise Pixels')
    ax4.set_ylabel('Flux')

    fig.savefig('simulted_xpos_ypos_npix_flux.png')

def nnbr_raw_one_data(yp0, xp0, np0, ypos, xpos, npix, n_nnbr=50):
    yp_diff = (yp0 - ypos)**2
    xp_diff = (xp0 - xpos)**2
    np_diff = (np0 - npix)**2
    
    rad_diff= np.sqrt(yp_diff + xp_diff + np_diff)
    
    return rad_diff.argsort()[1:n_nnbr+1] # skip the original point, which is always zero distant

n_nnbr  = 50
nCores  = cpu_count()-2
pool    = Pool(nCores)

partial_nnbr    = partial(nnbr_raw_one_data, ypos=ypos, xpos=xpos, npix=npix, n_nnbr=n_nnbr)

start = time()
raw_sorted_NNBR = pool.starmap(partial_nnbr, zip(ypos, xpos, npix))
print('{} points took {} seconds'.format(nPts, time()-start))

pool.close()
pool.join()

np.savetxt('raw_sorted_NNBR_from_xpos_ypos_npix_example.dat', np.array(raw_sorted_NNBR))
