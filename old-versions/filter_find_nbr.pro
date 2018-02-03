pro filter_find_nbr, apr, trim_time, noise_pix, ext, plnm, exptime, gain, fluxconv, photnoise


sfile='/Users/Brian/IDL/filter_nbr_out_'+ext+'.sav'


apr_sz=[2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.5,4.0,4.5,5.0]
apr_scale=[0.6,0.7,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.15,1.2]
apr_shift=[-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4]

  ext1=plnm + '_' 
  if noise_pix eq 1 then restore, '/Users/Brian/IDL/' +ext1 + 'var'
  if noise_pix eq 0 then restore, '/Users/Brian/IDL/' +ext1 + 'fixed'


;calculate photon noise limit
;channel 1
;CREATOR = 'S18.18.0'           / SW version used to create this FITS file
;FLUXCONV=               0.1253 ;/ Flux Conv. factor (MJy/sr per DN/sec)          
;GAIN    =                  3.7 ;/ e/DN conversion     
;EXPTIME =                 0.08 ;/ [sec] Effective integration time per pixel
;channel 2
;CREATOR = 'S18.18.0'           / SW version used to create this FITS file 
;FLUXCONV=               0.1469 ;/ Flux Conv. factor (MJy/sr per DN/sec)          
;GAIN    =                 3.71 ;/ e/DN conversion               
;EXPTIME =                 0.36 ;/ [sec] Effective integration time per pixel

;RAMP Correction to line


ftot=reform(flux(*,apr))
xpos=x_pos
ypos=y_pos
bkg_tot=bkd

bkg_tot=bkg_tot*!dpi*mean(apr_all(apr,*))^2
temp=ftot*exptime*gain/fluxconv
bkg_tot=bkg_tot*exptime*gain/fluxconv

photon_noise=sqrt(temp+bkg_tot)/temp
err_phot=median(photon_noise)
print, 'Photon noise limit', err_phot
photnoise=err_phot
nseg=1

;sort by bjd
bjd=bjd-2.455d6
ind=sort(bjd)
x_pos=x_pos(ind)
y_pos=y_pos(ind)
npix=npix(ind)
bjd=bjd(ind)
flux=flux(ind, apr)
;flux=flux(ind)

;;This section is clipping out an eclipse of T1b for testing 
ind=indgen(-13592+17342, start=13592)
x_pos=x_pos(ind)
y_pos=y_pos(ind)
npix=npix(ind)
bjd=bjd(ind)
flux=flux(ind, apr)


;clean any clear outliers
ind=where(abs(flux-median(flux)) lt median(flux))
x_pos=x_pos(ind)
y_pos=y_pos(ind)
npix=npix(ind)
bjd=bjd(ind)
flux=flux(ind)


;remove 'feature' in transit
;ind=where(bjd lt 575.94 or bjd gt 576.07, count)
;x_pos=x_pos(ind)
;y_pos=y_pos(ind)
;npix=npix(ind)
;bjd=bjd(ind)
;flux=flux(ind)

dt=bjd-min(bjd)
npoints=n_elements(bjd)

;find # of data segments
ftime=bjd-[bjd(npoints-1.0),bjd(0:npoints-2.0)]
ind_seg=where(abs(ftime) gt 10.d0/3600.d0,nseg)
ind_seg=[ind_seg, npoints]
print, 'apr #, # of segments', apr, nseg

for i=0,nseg-1 do begin
  nbeg=ind_seg(i)
  nend=ind_seg(i+1)-1.0
  fsub=flux(nbeg:nend)
  xsub=x_pos(nbeg:nend)
  ysub=y_pos(nbeg:nend)
  npsub=npix(nbeg:nend)
  bjdsub=bjd(nbeg:nend)
  dtsub=bjdsub-min(bjdsub)
;trim clear outliers
  ind=where(ABS(fsub-median(fsub)) le 6.d0*stddev(fsub,/NAN) and $
            ABS(xsub-median(xsub)) le 6.d0*stddev(xsub,/NAN) and $ 
            ABS(ysub-median(ysub)) le 6.d0*stddev(ysub,/NAN))
  print, i, stddev(fsub,/NAN), stddev(xsub,/NAN), stddev(ysub,/NAN)
  fsub=fsub(ind)
  xsub=xsub(ind)
  ysub=ysub(ind)
  npsub=npsub(ind)
  bjdsub=bjdsub(ind)
  dtsub=dtsub(ind)
;trim hot pixels
  nn=n_elements(fsub)
  mflux=median(fsub,50)
  mflux(0:24)=median(fsub(0:49))
  mflux(nn-25:nn-1)=median(fsub(nn-50:nn-1))
  ;fsdev=stddev(fsub(50:60))
  fsdev=100.d0
  mx=median(xsub,50)
  mx(0:24)=median(xsub(0:49))
  mx(nn-25:nn-1)=median(xsub(nn-50:nn-1))
  ;xsdev=stddev(xsub(50:60))
  xsdev=0.01d0
  my=median(ysub,50)
  my(0:24)=median(ysub(0:49))
  my(nn-25:nn-1)=median(ysub(nn-50:nn-1))
  ;ysdev=stddev(ysub(50:60))
  ysdev=0.01d0
  ind=where(ABS(fsub-mflux) le 3.d0*fsdev and $
            ABS(xsub-mx) le 3.d0*xsdev and $ 
            ABS(ysub-my) le 3.d0*ysdev and $ 
            dtsub gt trim_time)
          
  if (i eq 0) then begin
     xtot=xsub(ind)
     ytot=ysub(ind)
     nptot=npsub(ind)
     ftot=fsub(ind)
     bjd_tot=bjdsub(ind)
  endif else begin
;exclude transit segment
;     if (i ne 2) then begin
       xtot=[xtot,xsub(ind)]
       ytot=[ytot,ysub(ind)]
       nptot=[nptot,npsub(ind)]
       ftot=[ftot,fsub(ind)]
       bjd_tot=[bjd_tot,bjdsub(ind)]
;     endif
  endelse
endfor

err_tot=ftot*0.d0+stddev(ftot(50:100)/median(ftot))
time_tot=bjd_tot-min(bjd_tot)
flux_tot=ftot



;window,0
;plot, time_tot*24.0, sqrt(nptot), psym=3
;window,1
;plot, time_tot*24.0, xtot, psym=3, yrange=[14,16]
;window,2
;plot, time_tot*24.0, ytot, psym=3, yrange=[14,16]


print, '% trimmed', double(n_elements(flux)-n_elements(ftot))/double(n_elements(flux))*100.0

find_nbr_qhull, xtot, ytot, sqrt(nptot), 50.0, 1.0, 1.0, 1.0, nbr_ind, gw

save, xtot, ytot, nptot, flux_tot, bjd_tot, time_tot, nbr_ind, gw, err_tot, filename=sfile


end
