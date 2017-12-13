import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def main(plnm, aor, ramp):
	data = np.load('../'+plnm + '/' + aor + '/apr_fits/' + ramp + '/' + aor + '_mcmc_results.npz')
	print (str(plnm) + ', ' + str(aor) + ', with ramp ' + str(ramp))
	print('\n')
	print('BIC: ' + str(data['bic']))
	print('\n')
	print('T_0: ' + str(data['t0_mcmc']))
	print('R_p: ' + str(data['rp_mcmc']))
	print('\n')
	if (ramp=='exp'):
		print ('a_1: ' + str(data['a1_mcmc']))
		print ('a_2: ' + str(data['a2_mcmc']))
	elif (ramp=='slant'):
		print ('a_1: ' + str(data['a1_mcmc']))

	print('\n\n\n\n')

main(sys.argv[1], sys.argv[2], sys.argv[3])
