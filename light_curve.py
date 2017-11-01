import numpy as np
import matplotlib.pyplot as plt

def main():
	data = np.load('WASP_79/62173696.npz')
	lc = data['lc']
	time = data['time']
	print (time)
	n=250
	lcp = [sum(lc[i:i+n, 1])/n for i in range(0,len(lc[:,1]),n)]
	tp = [sum(time[i:i+n])/n for i in range(0,len(time),n)]
	plt.figure()
	axes = plt.gca()
	axes.set_ylim([0.95,1.05])
	r = len(lcp)
	plt.scatter(time[:], lc[:,1]/np.mean(lc[:,1]), s=1)
	plt.show()

main()