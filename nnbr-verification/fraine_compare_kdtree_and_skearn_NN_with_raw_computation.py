# coding: utf-8
import pyhull
from scipy.spatial     import cKDTree
from scipy.spatial     import qhull
from scipy.spatial     import qhull, Delaunay
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import KNeighborsClassifier

raw_sorted_fname= 'raw_sorted_NNBR_from_xpos_ypos_npix_example.dat'
raw_sorted_NNBR = np.loadtxt(raw_sorted_fname)

data_filename   = 'xpos_ypos_npix_example.dat'
xpos, ypos, npix, flux = np.loadtxt(data_filename).T

n_nbr   = 50
points  = np.array([xpos,ypos,npix])
neigh   = NearestNeighbors(n_neighbors=n_nbr)

neigh.fit(points.T)

kdtree  = cKDTree(points.T)

ind_kdtree  = kdtree.query(kdtree.data, n_nbr+1)[1][1:]
ind_sklearn = neigh.kneighbors(return_distance=False)

print('Do all KDTREE  indicies match RAW nn-indicies: {}'.format(not any(ind_kdtree[:,1:] - raw_sorted_NNBR[1:])))
print('Do all SKLEARN indicies match RAW nn-indicies: {}'.format(not any(ind_sklearn-raw_sorted_NNBR)))


def find_nbr_qhull(xpos, ypos, npix, sm_num = 100, a = 1.0, b = 1.0, c = 1.0, print_space = 10000.):
    '''
        Python Implimentation of N. Lewis method, described in Lewis etal 2012, Knutson etal 2012

        Taken from N. Lewis IDL code:

            Construct a 3D surface (or 2D if only using x and y) from the data
            using the qhull.pro routine.  Save the connectivity information for
            each data point that was used to construct the Delaunay triangles (DT)
            that form the grid.  The connectivity information allows us to only
            deal with a sub set of data points in determining nearest neighbors
            that are either directly connected to the point of interest or
            connected through a neighboring point

        Python Version:

        J. Fraine    first edition, direct translation from IDL 12.05.12
    '''
    neigh = NearestNeighbors(n_neighbors=sm_num) # need 1 more to exclude the test point itself
    #The surface fitting performs better if the data is scattered about zero
    
    npix    = np.sqrt(npix)
    
    x0  = (xpos - np.median(xpos))/a
    y0  = (ypos - np.median(ypos))/b
    
    if np.sum(npix) != 0.0 and c != 0:
        np0 = (npix - np.median(npix))/c
    else:
        if np.sum(npix) == 0.0:
            print('SKIPPING Noise Pixel Sections of Gaussian Kernel because Noise Pixels are Zero')
        if c == 0:
            print('SKIPPING Noise Pixel Sections of Gaussian Kernel because c == 0')
    
    k            = sm_num                           # This is the number of nearest neighbors you want
    n            = x0.size                          # This is the number of data points you have
    nearest      = np.zeros((k,n),dtype=np.int64)   # This stores the nearest neighbors for each data point
    
    #Multiplying by 1000.0 avoids precision problems
    if npix.sum() != 0.0 and c != 0:
        # kdtree  = cKDTree(np.transpose((y0*1000., x0*1000., np0*1000.)))
        neigh.fit(np.transpose((y0*1000., x0*1000., np0*1000.)))
    else:
        # kdtree  = cKDTree(np.transpose((y0*1000., x0*1000.)))
        neigh.fit(np.transpose((y0*1000., x0*1000.)))
    
    inds= neigh.kneighbors(return_distance=False)
    
    gw  = np.zeros((k,n),dtype=np.float64) # This is the gaussian weight for each data point determined from the nearest neighbors
    
    start   = time.time()
    for point in tqdm(range(n), total=n):
        # if np.round(point/print_space) == point/print_space: print_timer(point, n, start)
        
        # ind         = kdtree.query(kdtree.data[point],sm_num+1)[1][1:]
        ind         = inds[point]#neigh.kneighbors(points[point:point+1])[1].flatten()[1:]
        
        dx          = x0[ind] - x0[point]
        dy          = y0[ind] - y0[point]
        
        if npix.sum() != 0.0 and c != 0:
            dnp         = np0[ind] - np0[point]
        
        sigx        = np.std(dx )
        sigy        = np.std(dy )
        if npix.sum() != 0.0 and c != 0:
            signp       = np.std(dnp)
        if npix.sum() != 0.0 and c != 0:
            gw_temp     = np.exp(-dx**2./(2.0*sigx**2.)) * \
                          np.exp(-dy**2./(2.*sigy**2.))  * \
                          np.exp(-dnp**2./(2.*signp**2.))
        else:
            gw_temp     = np.exp(-dx**2./(2.0*sigx**2.)) * \
                          np.exp(-dy**2./(2.*sigy**2.))
        
        gw_sum      = gw_temp.sum()
        gw[:,point] = gw_temp/gw_sum
        
        if (gw_sum == 0.0) or ~np.isfinite(gw_sum):
            raise Exception('(gw_sum == 0.0) or ~isfinite(gw_temp))')
        
        nearest[:,point]  = ind
    
    return gw.transpose(), nearest.transpose() # nearest  == nbr_ind.transpose()

nPts = int(1e5)

xpos_nf = 0.35*sin(np.arange(0,nPts) / 1500 + 0.5) + 15
ypos_nf = 0.35*sin(np.arange(0,nPts) / 2000 + 0.7) + 15
npix_nf = 0.25*sin(np.arange(0,nPts) / 2500 + 0.4) + 15

flux_nf = 1+0.01*(xpos - xpos.mean()) + 0.01*(ypos - ypos.mean()) + 0.01*(npix - npix.mean())

gw   , nbr    = find_nbr_qhull(xpos   , ypos   , npix   ) 
gw_nf, nbr_nf = find_nbr_qhull(xpos_nf, ypos_nf, npix_nf)

gkr    = np.sum(flux[nbr]*gw,axis=1)
gkr_nf = np.sum(flux_nf[nbr_nf]*gw_nf,axis=1)

figure()
plot(flux_nf, '.', ms=1)
plot(gkr_nf, '.', ms=1)

figure()
plot(flux_nf - gkr_nf, '.', ms=1)

figure()
plot(flux - gkr, '.', ms=1)

figure()
plot(flux, '.', ms=1)
plot(gkr, '.', ms=1)