import matplotlib.pyplot as plt 
import numpy as np 
import os

def graph(plnm, aor, ramp):
	curr_dir = os.getcwd()
	k = np.load(os.getcwd() + '/' + plnm + '/' + aor + '/apr_fits/' + ramp + '/' + aor + '_summary.npy')
	plt.figure()
	plt.scatter(np.arange(20)/10.+1.4,k[:,3]/1e6)
	plt.title('Transit Depth vs Aperture Number for AOR ' + aor + ' with Ramp ' + ramp)
	plt.xlabel('Aperture #')
	plt.ylabel('Transit Depth (ppm)')
	plt.savefig('HAT_P_41/apertures/' + aor + ramp + '_.png')


graph('HAT_P_41', '62152704', 'exp')
graph('HAT_P_41', '62152704', 'none')
graph('HAT_P_41', '62152704', 'slant')



graph('HAT_P_41', '62153216', 'exp')
graph('HAT_P_41', '62153216', 'none')
graph('HAT_P_41', '62153216', 'slant')