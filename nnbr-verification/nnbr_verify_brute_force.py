import numpy as np
import matplotlib.pyplot as plt
import time
import csv
import pandas as pd

def calc_distance(a, b, c, curr_pt, next_pt):
	x_dist = a*(curr_pt[0]-next_pt[0])**2
	y_dist = b*(curr_pt[1]-next_pt[1])**2
	np_dist = c*(curr_pt[2]-next_pt[2])**2

	return np.sqrt(x_dist + y_dist + np_dist)

def getKey(item):
	return item[1]

def calc_nnbrs(a,b,c,curr_pt,xpos,ypos,npix):
	dists = []
	for i in range(0,len(xpos)):
		dists.append((i,calc_distance(a,b,c,curr_pt,(xpos[i],ypos[i],npix[i]))))
	s = sorted(dists,key=getKey)[1:51]
	return [k[0] for k in s]


def brute_force(a, b, c):
	quant = True
	xpos = np.genfromtxt('test_csv/xpos.csv', delimiter=',')
	xpos = xpos[1:]
	print(len(xpos))

	ypos = np.genfromtxt('test_csv/ypos.csv', delimiter=',')
	ypos = ypos[1:]
	print(len(ypos))

	npix = np.genfromtxt('test_csv/npix.csv', delimiter=',')
	npix = npix[1:]
	print(len(npix))

	print('finished reading')
	if (quant):
		phots = 1 + 0.01 * (xpos - np.mean(xpos)) + 0.01 * (ypos - np.mean(ypos)) + 0.01 * (npix - np.mean(npix))

		gw = np.genfromtxt('test_csv/gw-idl.csv',delimiter=',')
		gw_idl = np.transpose(gw)
		print(gw_idl.shape)
		nbr = np.genfromtxt('test_csv/nbr-idl.csv',delimiter=',', dtype='int')
		nbr_idl = np.transpose(nbr)
		print(nbr_idl.shape)

		gw = np.genfromtxt('test_csv/gw-python.csv',delimiter=',')
		gw_pyth = np.transpose(gw)
		print(gw_pyth.shape)
		nbr = np.genfromtxt('test_csv/nbr-python.csv',delimiter=',', dtype='int')
		nbr_pyth = np.transpose(nbr)
		print(nbr_pyth.shape)

		tru_nbrs = np.genfromtxt('test_csv/true_nbrs.csv',delimiter=',',dtype='int')
		print(tru_nbrs.shape)

		idl_arr= list(range(len(phots)))
		for i in range(len(phots)):
			curr_nbrs = list(phots[nbr_idl[i]])
			curr_sum = curr_nbrs*gw_idl[i]
			idl_arr[i] = np.sum(curr_sum)

		pyth_arr = list(range(len(phots)))
		for i in range(len(phots)):
			curr_nbrs = list(phots[nbr_pyth[i]])
			curr_sum = curr_nbrs*gw_pyth[i]
			pyth_arr[i] = np.sum(curr_sum)

		tru_arr = list(range(len(phots)))
		for i in range(len(phots)):
			curr_nbrs = list(phots[tru_nbrs[i]])
			curr_sum = curr_nbrs*gw_idl[i]
			tru_arr[i] = np.sum(curr_sum)

		diff_pyth = [x1 - x2 for (x1, x2) in zip(pyth_arr, tru_arr)]
		diff_idl = [x1 - x2 for (x1, x2) in zip(idl_arr, tru_arr)]

		print(np.linalg.norm(diff_pyth)) #0.08547720592495461
		print(np.linalg.norm(diff_idl)) #0.0

	else:
		all_nbrs_list=[]
		csvfile = open('test_csv/true_nbrs_upd.csv', 'w')
		writer = csv.writer(csvfile,delimiter=',')

		for i in range(0,len(xpos)):
			curr_pt = (xpos[i],ypos[i],npix[i])
			curr_nbr = calc_nnbrs(a,b,c,curr_pt,xpos,ypos,npix)
			writer.writerow(curr_nbr)
			if (i%500 == 0): print(str(i) + 'th Iteration\n')

	if quant:
		print("Done quantifying error")
	else:
		print("Done generating true neighbors")

def verify():
	idl = False 
	true_csv = pd.read_csv('test_csv/true_nbrs.csv',header=None)
	if (idl):
		idl = pd.read_csv('test_csv/nbr-idl.csv',header=None)
		for i in range(idl.shape[1]):
			curr = idl[i]
			curr = np.array(np.transpose(curr))
			tru_row = true_csv.iloc[i]
			true = np.array(tru_row)
			if (i%500 == 0): print('yo ' + str(i))
			if (set(true) != set(curr)):
				print(str(i) + 'was the disagreeing index. Take a look')
				return False
	else:
		pyth = pd.read_csv('test_csv/nbr-python.csv',header=None)
		count=0
		for i in range(pyth.shape[1]):
			curr = pyth[i]
			curr = np.array(np.transpose(curr))
			tru_row = true_csv.iloc[i]
			true = np.array(tru_row)
			if (i%500 == 0): print('At index ' + str(i))
			inters = set(true).intersection(set(curr))
			print(50 - len(inters))
		print(count)

brute_force(1.,1.,1.)
#verify()