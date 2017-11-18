import scipy.interpolate
import csv
import pandas
import numpy as np
from scipy.interpolate import griddata

def find_coeffs(temp, log_g, metallicity, channel, type_limb):
	if channel==1: 
		filename='limb-darkening_ch1.csv'
		sep='\s*,\s*'
		enc='ascii'
	elif channel==2: 
		filename='limb-darkening_ch2.csv'
		sep = '\t'
		enc='utf-8-sig'
	else: 
		print('invalid channel')
		return
	csv = pandas.read_csv(filename, sep, header=0, encoding=enc, engine='python')
	temp_arr = csv['Temperature']
	log_g_arr = csv['Log G']
	met = csv['Metallicity']
	xi=np.array([temp, log_g, metallicity])
	points=np.transpose(np.array([temp_arr, log_g_arr, met]))
	if type_limb=='linear':  
		values=np.transpose(np.array([csv['Linear']]))
		values=np.squeeze(values)
		coeffs = griddata(points, values, xi, method='linear')
	if type_limb=='quadratic':
		coeffs=np.zeros(2)
		values=np.transpose(np.array([csv['Quadratic A']]))
		values=np.squeeze(values)
		coeffs[0] = griddata(points, values, xi, method='linear')
		values=np.transpose(np.array([csv['Quadratic B']]))
		values=np.squeeze(values)
		coeffs[1] = griddata(points, values, xi, method='linear')
	if type_limb=='nonlinear 3':
		coeffs=np.zeros(3)
		values=np.transpose(np.array([csv['NL3 c2']]))
		values=np.squeeze(values)
		coeffs[0] = griddata(points, values, xi, method='linear')
		values=np.transpose(np.array([csv['NL3 c3']]))
		values=np.squeeze(values)
		coeffs[1] = griddata(points, values, xi, method='linear')
		values=np.transpose(np.array([csv['NL3 c4']]))
		values=np.squeeze(values)
		coeffs[2] = griddata(points, values, xi, method='linear')
	if type_limb=='nonlinear 4':
		coeffs=np.zeros(4)
		values=np.transpose(np.array([csv['NL4 c1']]))
		values=np.squeeze(values)
		coeffs[0] = griddata(points, values, xi, method='linear')
		values=np.transpose(np.array([csv['NL4 c2']]))
		values=np.squeeze(values)
		coeffs[1] = griddata(points, values, xi, method='linear')
		values=np.transpose(np.array([csv['NL4 c3']]))
		values=np.squeeze(values)
		coeffs[2] = griddata(points, values, xi, method='linear')
		values=np.transpose(np.array([csv['NL4 c4']]))
		values=np.squeeze(values)
		coeffs[3] = griddata(points, values, xi, method='linear')

	return coeffs

def main():
	print(find_coeffs(5234., 4.460, 0.359, 2, 'NL3'))

#main()
