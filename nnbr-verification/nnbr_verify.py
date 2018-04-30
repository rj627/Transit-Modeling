import numpy as np
from numpy import genfromtxt
import numpy.random as npr
import matplotlib.pyplot as plt
from find_nbr import *
import csv 
import os

def verify():
	arr = list(range(100000))
	read = 't'
	mode = 'i'

	if (read == 'f'):
		xpos = list(range(100000))
		ypos = list(range(100000))
		npix = list(range(100000))
		for i in arr:
			xpos[i] = 0.35*np.sin(i/1500 + 0.5) + 15
			ypos[i] = 0.35*np.sin(i/2000 + 0.7) + 15
			npix[i] = 0.25*np.sin(i/2500 + 0.4) + 15
		
		xpos = [x*(1+npr.randn()*0.01) for x in xpos]
		ypos = [y*(1+npr.randn()*0.01) for y in ypos]
		npix = [n*(1+npr.randn()*0.01) for n in npix]

		csvfile = open('xpos.csv', 'w')
		fieldnames = ['xpos']
		writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

		writer.writeheader()
		for i in xpos:
			writer.writerow({'xpos': str(i)})


		csvfile2 = open('ypos.csv', 'w')
		fieldnames = ['ypos']
		writer = csv.DictWriter(csvfile2, fieldnames=fieldnames)

		writer.writeheader()
		for i in ypos:
			writer.writerow({'ypos': str(i)})


		csvfile3 = open('npix.csv', 'w')
		fieldnames = ['npix']
		writer = csv.DictWriter(csvfile3, fieldnames=fieldnames)

		writer.writeheader()
		for i in npix:
			writer.writerow({'npix': str(i)})

	else:
		xpos = genfromtxt('test_csv/xpos.csv', delimiter=',')
		xpos = xpos[1:]
		print(len(xpos))

		ypos = genfromtxt('test_csv/ypos.csv', delimiter=',')
		ypos = ypos[1:]
		print(len(ypos))

		npix = genfromtxt('test_csv/npix.csv', delimiter=',')
		npix = npix[1:]
		print(len(npix))
		print('finished reading')

	phots = 1 + 0.01 * (xpos - np.mean(xpos)) + 0.01 * (ypos - np.mean(ypos)) + 0.01 * (npix - np.mean(npix))

	#plot(arr,xpos, ypos, npix,phots)
	#exit(0)

	############################################

	if (mode == 'p'):
		gw, nbr = find_nbr_qhull(xpos, ypos, npix, sm_num=50)
		print(gw.shape)
		print(nbr.shape)
		gw_csv = open('test_csv/gw-python.csv', 'x')
		writer = csv.writer(gw_csv)
		writer.writerows(np.transpose(gw))

		nbr_csv = open('test_csv/nbr-python.csv', 'x')
		writer = csv.writer(nbr_csv)
		writer.writerows(np.transpose(nbr))

	else:
		gw = np.genfromtxt('test_csv/gw-idl.csv',delimiter=',')
		gw = np.transpose(gw)
		print(gw.shape)
		nbr = np.genfromtxt('test_csv/nbr-idl.csv',delimiter=',', dtype='int')
		nbr = np.transpose(nbr)
		print(nbr.shape)
	
	arr2 = list(range(len(phots)))
	for i in range(len(phots)):
		curr_nbrs = list(phots[nbr[i]])
		curr_sum = curr_nbrs*gw[i]
		arr2[i] = np.sum(curr_sum)

	csvfile2 = open('old_arr2-i.csv', 'w')
	fieldnames = ['old_pos']
	writer = csv.DictWriter(csvfile2, fieldnames=fieldnames)

	writer.writeheader()
	for i in ypos:
		writer.writerow({'old_pos': str(i)})

	plt.scatter(bin_anything(arr, 100), bin_anything(phots, 100), s=1)
	plt.scatter(bin_anything(arr, 100), bin_anything(arr2, 100), s=1)
	plt.savefig(os.getcwd()+'/nnbr-verification/model-'+mode+'.png')
	plt.close()

	res = phots-arr2
	plt.figure()
	plt.scatter(arr, res, s=1)
	plt.savefig(os.getcwd()+'/nnbr-verification/resids-'+mode+'.png')

	est_rednoise(res, 3.5, mode)
	print('done with red noise')

def plot(t,x,y,n,p):
	plt.figure()
	plt.subplot(411)
	plt.title('X pos, Y pos, Noise Pixel, Light Curve - Model Data')
	plt.scatter(t, x, s=1)
	plt.ylim(14, 16)
	plt.ylabel('X position')
	plt.xticks([])
	plt.subplot(412)
	plt.scatter(t,y, s=1)
	plt.ylim(14, 16)
	plt.ylabel('Y position')
	plt.xticks([])
	plt.subplot(413)
	plt.scatter(t, n, s=1)
	plt.ylim(14, 16)
	plt.ylabel('Sqrt Noise Pixel')
	plt.xticks([])
	plt.subplot(414)
	plt.scatter(t, p, s=1)
	plt.ylim(0.97, 1.03)
	plt.ylabel('Light Curve (Flux)')
	plt.xlabel('Time')
	plt.savefig(os.getcwd()+'/nnbr-verification/gaussian-xyb_plot')
	print("saved plot")

def est_rednoise(resids, exptime, mode):
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

    tscales=np.where((np.log10(bin_time)>2.5) & (np.log10(bin_time)<3.0))
    
    bred=np.nanmedian(beta_red[tscales])
    if bred < 1.: bred=np.amax(beta_red[tscales])
    if bred < 1.: bred=1.
    print(bred)
    print(bin_std[0])
    
    plt.figure()
    if mode == 'p':
    	plt.title('Red Noise - Model Data (Python)')
    else:
    	plt.title('Red Noise - Model Data (IDL)')
    plt.scatter(np.log10(bin_time), np.log10(bin_std), color='red')
    plt.plot(np.log10(bin_time),y ,linewidth=3.3)
    plt.text(1.5,-2.5, 'Beta_red:  '+str(round(bred, 2)))
    plt.xlabel('Bin Size (log (s))')
    plt.ylabel('Log Std Dev of Residuals')
    plt.savefig(os.getcwd()+'/nnbr-verification/red-noise-' + mode + '.png')
    plt.draw()
    plt.pause(1.1)

    return  bin_std[0], bred

def bin_anything(series, binsize):
	series=np.squeeze(series)
	numbins=int(np.floor(series.size/binsize))
	bin_series=np.zeros(numbins)
	for b in range(0,numbins):
		bin_series[b]=np.mean(series[b*binsize:(b+1)*binsize-1])

	return bin_series


verify()