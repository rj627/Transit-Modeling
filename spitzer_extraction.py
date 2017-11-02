import numpy as np
from photutils import aperture_photometry
from photutils import CircularAperture
from astropy.io import fits
import sys
import glob
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy.ma as ma
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
import os
from photutils import centroid_com, centroid_1dg, centroid_2dg
from exoplanets_dot_org import *
import pandas as pd


def main(plnm, channel, aornum):

    
    if plnm == 'HD209': plnm2='HD 209458 b'
    if plnm == 'HD189': plnm2='HD 189733 b'
    if plnm== 'W16': 
        plnm2='WASP-16 b'
        basepath='/Users/Brian/Desktop/Tucker_Group/Spitzer/mapping_files/'
    
    #basepath='/Users/Brian/Desktop/Tucker_Group/t_1/'
    #basepath='/Users/Brian/Desktop/Tucker_Group/Spitzer/'+plnm+'/'
    #basepath='/home/bkilpatr/mapping_files/'

    if plnm== 'WASP_79': 
        plnm2='WASP-79 b'
        basepath='/Users/rahuljayaraman/Documents/Miscellany/Research (Tucker Group)/Python (Transits)/' + plnm + '/'

    chnum=channel#1
  
    extract='on'
    gaussiancent='on' 
    
    
    radii=np.arange(10)/5.+2.00
    var_rad=np.arange(10)/10.+1.0
    var_add=np.arange(10)/5-0.4
    if aornum==[]:aornum=aor_from_list(plnm, chnum, basepath)

    print(aornum)
    print (plnm)
   
    if extract=='on':  
        extraction(aornum, chnum, plnm, plnm2, radii, var_rad, var_add, gaussiancent, basepath)


def extraction(aornum, chnum, plnm, plnm2, radii, var_rad, var_add, gaussiancent, basepath):

    for aor in aornum:
        aor=str(aor)
        paths = glob.glob(basepath+'r'+aor+ '/ch*')
        fpathout='/Users/rahuljayaraman/Documents/Miscellany/Research (Tucker Group)/Python (Transits)/'+plnm+'/'+aor[-8:]
        directory = os.path.dirname(fpathout)
        if not os.path.exists(directory):
            os.makedirs(directory)

        fpath=basepath+'r' + aor + '/ch' + str(chnum) + '/bcd/'

        filenames=glob.glob(fpath+ '*bcd.fits')
        filenames.sort()
        for i in filenames:
            print (i)

        nframes=len(filenames)
        print (fpath)
     
        hold_pos=np.zeros(shape=(nframes*64, 2))
        central_pix=np.zeros(shape=(5,5, nframes*64))
        #all_rad=np.concatenate((radii, var_rad))
        lightcurve=np.zeros(shape=(nframes*64, len(radii)+len(var_rad)+len(var_add)))
        time=np.zeros(nframes*64)
        beta_np=np.zeros(nframes*64)
        
        
        for i in range(0,nframes):
            if i % 10==0: 
                os.system('clear')
                print(aor, i,' of ',str(nframes))
                # if i !=0:
                #     plt.close('all')
                #     plt.figure()
                #     plt.scatter((time[0:(i-1)*64]-time[0])*24.0*60., hold_pos[0:(i-1)*64,1], s=1, color='r')
                #     plt.scatter((time[0:(i-1)*64]-time[0])*24.0*60., hold_pos[0:(i-1)*64,0]-0.5, s=1)
                #     plt.ylim(14.2, 15.5)
                #     plt.draw()
                #     plt.pause(0.1)

            hdulist = fits.open(filenames[i])
            channel=str(hdulist[0].header['CHNLNUM'])
            
            cube=hdulist[0].data
            exptime=hdulist[0].header['EXPTIME']
            framtime=hdulist[0].header['FRAMTIME']
            mmbjd=hdulist[0].header['BMJD_OBS']
            for j in range(0,64):    
                scidata=cube[j,:,:]
                bkgd= backgr(scidata)
                data=scidata-bkgd  
                data=ma.masked_invalid(data)
                
                bnp1=np.sum(data)**2
                bnp2=np.sum(np.multiply(data,data))
                bnp=bnp1/bnp2
        
                xc,yc=centroid(data, gaussiancent)
                position=[xc,yc]
                beta_np[64*i+j]=bnp
                hold_pos[64*i+j, :]=position  
                vrad1=var_rad*np.sqrt(bnp)
                vrad2=var_add+np.sqrt(bnp)
                vrad=np.concatenate((vrad1, vrad2))
                all_rad=np.concatenate((radii, vrad))
                apertures = [CircularAperture(position, r=r) for r in all_rad]

                phot_table = aperture_photometry(scidata, apertures)

                for q in range (0,len(all_rad)):
                    if q ==0: phot=np.zeros(len(all_rad))
                    phot[q]=phot_table.columns[q+3]  
                time[64*i+j]=mmbjd+framtime*j/86400.0
                lightcurve[64*i+j, :]=phot
                central_pix[:,:,64*i+j]=data[13:18, 13:18]
            hdulist.close()
        orbparams=get_orb_pars(plnm2, basepath)

        #pmap=pmap_corr(hold_pos, channel)

     
        np.savez(fpathout+'extraction',  ch=channel, time=time, lc=lightcurve, cp=central_pix,  op=orbparams, exptime=exptime, framtime=framtime, beta_np=beta_np, hold_pos=hold_pos, all_rad=all_rad)
        t=time
        npix=beta_np

        plt.figure()
        plt.subplot(311)
        plt.title(plnm+' Ch: '+str(chnum)+'\n'+aor)
        plt.scatter(t, hold_pos[:,0], s=1)
        plt.ylim(14.5, 15.5)
        plt.ylabel('X position')
        plt.xticks([])
        plt.subplot(312)
        plt.scatter(t, hold_pos[:,1], s=1)
        plt.ylim(14.5, 15.5)
        plt.ylabel('Y position')
        plt.xticks([])
        plt.subplot(313)
        plt.scatter(t, np.sqrt(npix), s=1)
        plt.ylim(2, 3)
        plt.ylabel('Sqrt Noise Pixel')
        plt.xlabel('Time')
        plt.savefig(fpathout+'xyb_plot')


        #send_mail('xyb_plot.png', fpathout, aor)

    return None


def backgr(a):
    
    backmask=np.zeros(shape=(32,32))
    backmask[10:20, 10:20]=1
    mean, median, std = sigma_clipped_stats(a, sigma=5.0, mask=backmask)

    return median

def centroid(a, cent):
    
    top=17
    bot=14
    
    a=sigma_clip(a, sigma=7, iters=1)    
    a=a[bot:top+1, bot:top+1]
    a2=np.multiply(a,a)
    beta_np=np.sum(a)**2/np.sum(a2)
    if cent=='on':
        xc, yc = centroid_2dg(a)+bot
    else:
        xc, yc = centroid_com(a)+bot
    return (xc,yc)    

# def get_aor(plnm, plnm2, t1, basepath):
#     #a=find_all_dir('/Users/Brian/Desktop/Tucker_Group/Spitzer/'+plnm2)
#     if t1 != 'true':
#         #a=os.listdir('/Users/Brian/Desktop/Tucker_Group/Spitzer/'+plnm)
#         basepath='/Users/Brian/Desktop/Tucker_Group/Spitzer/'+plnm+'/'
#         a=glob.glob(basepath+ 'r*')
        

      
#     else: 
#         #basepath='/Users/Brian/Desktop/Tucker_Group/t_1/'
#         a=glob.glob(basepath+ 'r*')
#         nn=[ name for name in os.listdir(basepath) if os.path.isdir(os.path.join(basepath, name)) ]
#         del nn[0]    

        
#     print(a)
#     sys.exit()
#     return a
def pmap_corr(position,  channel):


    from scipy.interpolate import griddata

    fpath='/Users/Brian/Documents/Python_Scripts/Spitzer/pmap_fits/'
    if channel == '1':
        hdulist=fits.open(fpath+'xgrid_ch1_500x500_0043_120828.fits')
        xgrid=hdulist[0].data
        hdulist.close()
        hdulist=fits.open(fpath+'ygrid_ch1_500x500_0043_120828.fits')
        ygrid=hdulist[0].data
        hdulist.close()
        hdulist=fits.open(fpath+'pmap_ch1_500x500_0043_120828.fits')
        pmap=hdulist[0].data
        hdulist.close()
    
    if channel=='2':
        hdulist=fits.open(fpath+'xgrid_ch2_0p1s_x4_500x500_0043_120124.fits')
        xgrid=hdulist[0].data
        hdulist.close()
        hdulist=fits.open(fpath+'ygrid_ch2_0p1s_x4_500x500_0043_120124.fits')
        ygrid=hdulist[0].data
        hdulist.close()
        hdulist=fits.open(fpath+'pmap_ch2_0p1s_x4_500x500_0043_120124.fits')
        pmap=hdulist[0].data
        hdulist.close()


    x=np.ravel(xgrid)
    y=np.ravel(ygrid)
    values=np.ravel(pmap)
    points=np.transpose(np.array([x,y]))
    if position.shape[0] != 2: position=np.transpose(position)
    corr=np.zeros(position.shape[1])
  
    corr=griddata(points, values, np.transpose(position), method='linear')
    
    return corr

def aor_from_list(planet, ch, basepath):
    filename=basepath+'aor_list.csv'
    #filename='/home/bkilpatr/Spitzer_Routines/aor_list.csv'
    aor_arr=pd.read_csv(filename)
    chlist=aor_arr.CH
    aorlist=aor_arr.AOR
    namelist=aor_arr.Target
    ch_aors=np.where((namelist==planet) & (chlist==ch))
    aor_selected=np.array(aorlist[np.squeeze(ch_aors)])


    return aor_selected

def send_mail(filename, filepath, aor, servername, email, pwd):
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart
    from email.mime.base import MIMEBase
    from email import encoders



    fromx = 'brian.m.kilpatrick@gmail.com'
    to  = 'brian_kilpatrick@brown.edu'
    msg = MIMEMultipart()
    #msg = MIMEText('Python test')
    msg['Subject'] = 'Extraction Complete: '+aor
    msg['From'] = fromx
    msg['To'] = to


    body='The extraction for AOR '+aor+' is now complete.  Plots are attached.\n'

    msg.attach(MIMEText(body, 'plain'))
    #filename="t1_sys.png"
    #filepath="/home/bkilpatr/mapping_files/"
    attachment= open(filepath+filename, 'rb')
    p = MIMEBase('application', 'octet-stream')
     
    # To change the payload into encoded form
    p.set_payload((attachment).read())
     
    # encode into base64
    encoders.encode_base64(p)
      
    p.add_header('Content-Disposition', "attachment; filename= %s" % filename)
     
    # attach the instance 'p' to instance 'msg'
    msg.attach(p)

    server = smtplib.SMTP(servername)
    server.starttls()
    server.ehlo()
    server.login(email, pwd)
    server.sendmail(fromx, to, msg.as_string())
    server.quit()

    return None



main(sys.argv[1], int(sys.argv[2]), [sys.argv[3]])
