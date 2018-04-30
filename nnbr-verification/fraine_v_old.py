import numpy as np

fr_gkr= np.genfromtxt('fraine-gkr.csv',delimiter=',')
fr_gkr = np.transpose(fr_gkr)

print(fr_gkr.shape)

old_vals = np.genfromtxt('old_arr2-i.csv',delimiter=',')
old_vals = np.transpose(old_vals)
print(old_vals.shape)

ct = 0
for i in range(len(fr_gkr)):
	if (old_vals[i] != fr_gkr [i]):
		print(old_vals[i], fr_gkr[i])

print(ct)
