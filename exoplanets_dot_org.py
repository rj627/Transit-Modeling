import pandas as pd
import numpy as np
import sys
import os

def exoplanet_import_web():
  #Function that downloads and imports the latest CSV file from exoplanets.org:
  
  exoarr=pd.read_csv('http://exoplanets.org/csv-files/exoplanets.csv')
  
  return exoarr
  
def exoplanet_import_file(basepath):
  #Function that imports CSV file from HD:
  
  exoarr=pd.read_csv('/Users/rahuljayaraman/Documents/Miscellany/Research (Tucker Group)/exoplanet.csv')
  #exoarr=pd.read_csv(basepath+'exoplanets.csv')
  
  return exoarr
  

def get_orb_pars(plnm, basepath):
  #Now define the main, which will call all the above functions and determine 
  #exoplanets fit the criteria.

  
  #Import the latest exoplanets.org dataset and dictionary of units, definitions.

  exoarr=exoplanet_import_web()
  #exoarr=exoplanet_import_file(basepath)


  
  #Call needed columns from exoplanets.org dataset
  decs=exoarr.DEC #Declination of planets (Deg)
  ra=exoarr.RA #Right ascention of planets (Hours, decimal)
  tt=exoarr.TT #Time of transit center (JD)
  per=exoarr.PER #Period of transit (Days)
  trans=exoarr.TRANSIT #1: transit observed; 0: no transit observed
  dur=exoarr.T14 #Duration of transit, days
  name=exoarr.NAME #Exoplanet name
  depth=exoarr.DEPTH #Depth of transit, in (Rp/Rs)^2
  dist=exoarr.DIST #Distance to star, parsecs
  teff=exoarr.TEFF #Effective temp. of star, K
  rstar=exoarr.RSTAR #Radius of star, solar radii
  starmag=exoarr.V #Visual magnitude of star
  starhipp=exoarr.HIPP #hipparcos catalogue # of stars
  sep=exoarr.SEP #separation from planet to star (AU)
  ar=exoarr.AR #separation from planet to star (a/R*)
  inc=exoarr.I
  ecc=exoarr.ECC
  met=exoarr.FE
  logg=exoarr.LOGG
  index=np.where(name == plnm)
  index=index[0]

  orbpars=np.zeros(11)
  orbpars[5]=per[index]
  orbpars[6]=tt[index]
  orbpars[2]=depth[index]**0.5
  orbpars[3]=depth[index]
  orbpars[0]=ar[index]
  orbpars[1]=inc[index]
  orbpars[7]=ecc[index]
  orbpars[4]=dur[index]
  orbpars[8]=met[index]
  orbpars[9]=logg[index]
  orbpars[10]=teff[index]




  return orbpars

