import numpy as np
import matplotlib.pyplot as plt
import sys

def main(plnm, aor, ramp):
	data = np.load(plnm + '/' + aor + '/apr_fits/' + ramp + '/' + aor + '_mcmc_results.npz')
	print (str(plnm) + ', ' + str(aor) + ', with ramp ' + str(ramp))
	print('T0: ' + str(data['t0_mcmc']))
	print('Rp: ' + str(data['rp_mcmc']))
	print('BIC: ' + str(data['bic']))
	if (ramp=='exp'):
		print ('a1: ' + str(data['a1_mcmc']))
		print ('a2: ' + str(data['a2_mcmc']))
	elif (ramp=='slant'):
		print ('a1: ' + str(data['a1_mcmc']))

	print('\n\n\n\n')

main(sys.argv[1], sys.argv[2], sys.argv[3])
