import numpy as np
import os
import corner

def main():
	s = np.load(os.getcwd() + '/WASP_101/62159360/apr_fits/exp/62159360_samples.npy')
	print(len(s))
	fig=corner.corner(s,labels=[ "t0", 'Rp', "a1", "a2"])#, "A/R", "inc"])
	fig.savefig('/Users/rahuljayaraman/Desktop/62159360_corner.png')

main()