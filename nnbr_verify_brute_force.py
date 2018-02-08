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

	all_nbrs_list=[]
	csvfile = open('test_csv/true_nbrs.csv', 'w')
	writer = csv.writer(csvfile,delimiter=',')
	
	for i in range(0,len(xpos)):
		curr_pt = (xpos[i],ypos[i],npix[i])
		curr_nbr = calc_nnbrs(a,b,c,curr_pt,xpos,ypos,npix)
		writer.writerow(curr_nbr)
		if (i%500 == 0): print(str(i) + 'th Iteration\n')

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

#brute_force(1.,1.,1.)
verify()