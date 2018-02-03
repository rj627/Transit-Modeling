pro find_nbr_qhull, xpos, ypos, np, sm_num, a, b, c, nbr_ind, gw

;construct a 3D surface (or 2D if only using x and y) from the data
;using the qhull.pro routine.  Save the connectivity information for
;each data point that was used to construct the Delaunay triangles (DT)
;that form the grid.  The connectivity information allows us to only
;deal with a sub set of data points in determining nearest neighbors
;that are either directly connected to the point of interest or
;connected through a neighboring point

;a=1.0d0
;b=0.75d0
;c=1.0d0

;The surface fitting performs better if the data is scattered about zero
x0=a*(xpos-median(xpos))
y0=b*(ypos-median(ypos))
np0=c*(np-median(np))

;Multiplying by 1000.0 avoids precision problems
qhull, x0*1000.d0, y0*1000.d0, np0*1000.d0, Tr, connectivity=c, /DELAUNAY  
;qhull, x0*1000.d0, y0*1000.d0, Tr, connectivity=c, /DELAUNAY

k=sm_num ;This is the number of nearest neighbors you want
n=n_elements(x0) ;This is the number of data points you have
nearest=lonarr(k,n,/NOZERO) ;This stores the nearest neighbors for each data point
nearest_d=fltarr(k,n,/NOZERO) ;This stores the distance to each of the nearest neighbors
gw=fltarr(n, sm_num) ;This is the gaussian weight for each data point determined from the nearest neighbors
for point=0L,n-1 do begin
     ;if ROUND(point/10000.) eq point/10000. then print,point
     p=c[c[point]:c[point+1]-1] ;start with this point's DT neighbors
     d=(x0[p]-x0[point])^2+(y0[p]-y0[point])^2+(np0[p]-np0[point])^2 ;calculate the distance to each of the DT neighbors
     for i=1,k do begin 
        s=sort(d)               ; sort the DT neighbors by distance
        p=p[s] & d=d[s]
        nearest[i-1,point]=p[0] ;for i=1 the nearest point is simply the nearest DT neighbor 
        nearest_d[i-1,point]=d[0]
        if i eq k then continue
        
        ;; Add all its neighbors not yet seen
        new=c[c[p[0]]:c[p[0]+1]-1]  ;use the ith nearest neighbor to search for other neighbors
        nnew=n_elements(new) 
        ;The already_got vector includes the nearest DT neighbors (p), the nearest i
        ;neighbors that have already been determined, and the point itself 
        already_got=[p,nearest[0:i-1,point],point] 
        ngot=n_elements(already_got) 
        ;here we compare the points that we already have as a possible nearest neighbors (already_got) with the 
        ;neighbors of the ith nearest neighbor (new).  This could also be done with a loop where each 
        ;element of new is compared with already_got, but this method is faster. 
        wh=where(long(total(rebin(already_got,ngot,nnew,/SAMPLE) ne $
                            rebin(transpose(new),ngot,nnew,/SAMPLE),1)) $
                 eq ngot, cnt)
        if cnt gt 0 then begin ;if any of the ith points nearest neighbors has not already been included as a possible nearest neighbor
           new=new[wh] 
           p=[p[1:*],new] ;form a new p vector and exclude the ith nearest neighbor p[0] and add in new possible nearest neighbors
           d=[d[1:*],(x0[new]-x0[point])^2+(y0[new]-y0[point])^2+(np0[new]-np0[point])^2]
        endif else begin 
           p=p[1:*] ;if no new possible nearest neighbors are found, simple exclude the ith nearest neighbor p[0] and continue
           d=d[1:*]
        endelse 
     endfor
   ind=nearest(*,point)
   dx=xpos(ind)-xpos(point)
   dy=ypos(ind)-ypos(point)
   dnp=np(ind)-np(point)
   sigx=stddev(dx,/DOUBLE,/NAN)
   sigy=stddev(dy,/DOUBLE,/NAN)
   signp=stddev(dnp,/DOUBLE,/NAN)
   gw_temp=exp(-dx^2./(2.d0*sigx^2))*exp(-dnp^2./(2.d0*signp^2))*exp(-dy^2./(2.d0*sigy^2)) 
   gw(point,*)=gw_temp/total(gw_temp)
   if (total(gw_temp) eq 0.0 or finite(total(gw_temp),/NAN)) then stop
  endfor 
  nearest_d=sqrt(nearest_d)
  nbr_ind=transpose(nearest)


end
