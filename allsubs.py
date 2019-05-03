#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 16:46:05 2018

@author: zhitu
"""

import numpy as np
from obspy.geodetics import locations2degrees
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq
from scipy.fftpack import rfft,rfftfreq,irfft

def getsacname(network,station,starttime,endtime,channel):
    return network+'.'+station+'.'+starttime+'.'+endtime+'.'+channel+'.SAC'


def Getwins(inventory,xlat1,xlon1,t1,U0=4.0):
    
    """ calculate the windows that does not contain an earthquake in the catalog
    Input
      an obspy inventory
      xlat1,xlon1 for the station to be examined
      UTCtime for the start time of the series t1"""
    
    win_len=2000
    wins=np.arange(0,86400-2000,win_len)
    nwin=len(wins)
    idx=np.ones(nwin,dtype=bool)
    eq_labels=[]
    for ev in inventory:
        t0=ev.origins[0].time
        xlat0=ev.origins[0].latitude
        xlon0=ev.origins[0].longitude
        dist=locations2degrees(xlat0,xlon0,xlat1,xlon1)
        time=dist*111.1949/U0
        if (t0+time < t1):
            continue
        if (t0+time >= t1+86400):
            break
        ijk=np.floor((t0+time-t1)/win_len).astype(int)
        if (ijk <= nwin-1):
            idx[ijk]=False
        idx[np.maximum(1,ijk-1)]=False
        idx[np.minimum(nwin-1,ijk+1)]=False
        eq_labels.append(t0+time-t1)
        print(ev.origins[0].time,dist)
    
    return wins[idx],eq_labels


def Caltransfer(y1,y2,wins,nlen=2000,iopt=1):
    """ calculate the transfer function from y1 to y2
    return the coherence, admittance, phase and their corresponding error 
    if iopt==0, then only coherence is returned """
    
    coh_debug=[]
    win_debug=[]
    
    for ijk,win in enumerate(wins):
        y1tmp=y1[win:win+nlen]
        y2tmp=y2[win:win+nlen]
        hann=np.hanning(nlen)
        y1_fft=np.split(fft(hann*y1tmp),2)[0]
        y2_fft=np.split(fft(hann*y2tmp),2)[0]
        if (ijk == 0):
            Gxy=np.conj(y1_fft)*y2_fft
            Gxx=np.conj(y1_fft)*y1_fft
            Gyy=np.conj(y2_fft)*y2_fft
        else:
            Gxy=Gxy+np.conj(y1_fft)*y2_fft
            Gxx=Gxx+np.conj(y1_fft)*y1_fft
            Gyy=Gyy+np.conj(y2_fft)*y2_fft
    
        ff=np.split(fftfreq(nlen,1.0),2)[0]   
        idx=(ff>0.005) & (ff<0.010)
        cohtmp=np.abs(Gxy)**2/np.real(Gxx)/np.real(Gyy)
        cohtmp=np.sqrt(cohtmp)
        coh_debug.append(np.mean(cohtmp[idx]))
        win_debug.append(win)

            
    coh=np.abs(Gxy)**2/np.real(Gxx)/np.real(Gyy)
    coh=np.sqrt(coh)
    if (iopt == 0):
        adm=0.
        phs=0.
        adm_err=0.
        phs_err=0.
    else:
        adm=np.abs(Gxy)/np.real(Gxx)
        phs=np.angle(Gxy)
        nd=len(wins)
        adm_err=np.sqrt(1.-coh**2)/coh/np.sqrt(2*nd)
        adm_err=adm*adm_err
        phs_err=adm_err
    ff=np.split(fftfreq(nlen,1.0),2)[0]
    
    plt.plot(win_debug,coh_debug,'o')
    
    
    return ff,coh,adm,phs,adm_err,phs_err


def Remove(tr1,tr2,adm,adm_err,phs,phs_err,f1,f2,ff,iplot=0):
    """ calculate a quadratic fit to adm and phs
    use this information to predict from tr1, then remove this from tr2
    returning two trace (obspy class), one is the prediction
    one is this prediction removed from tr2 """
    idx=(ff>f1) & (ff<f2)
    ff_select=ff[idx]
    adm_select=adm[idx]
    adm_err_select=adm_err[idx]
    w=1./adm_err_select
    apol=np.polyfit(ff_select,adm_select,2,w=w)
    phs_select=phs[idx]
    phs_err_select=phs_err[idx]
    w=1./phs_err_select
    ppol=np.polyfit(ff_select,phs_select,2,w=w)
    
    if (iplot==1):
        plt.subplot(1,2,1)
        adm_fit=apol[0]*ff_select**2+apol[1]*ff_select+apol[2]
        plt.plot(ff_select,adm_select)
        plt.plot(ff_select,adm_fit)
        plt.subplot(1,2,2)
        phs_fit=ppol[0]*ff_select**2+ppol[1]*ff_select+ppol[2]
        plt.plot(ff_select,phs_select)
        plt.plot(ff_select,phs_fit)
        plt.show()
        plt.close()
    
    ffr=rfftfreq(len(tr1.data),1.0)
    tr_pred=tr1.copy()
    tr_left=tr1.copy()
    Htmp_spec=rfft(tr1.data)
    Htmp_spec[0]=0
    Htmp_spec[-1]=0
    for i in np.arange(1,len(ffr)-1,2):
        rp=Htmp_spec[i]
        ip=Htmp_spec[i+1]
        if(ffr[i]>f2 or ffr[i]<f1):
            Htmp_spec[i]=0.
            Htmp_spec[i+1]=0.
            continue
        amp=apol[0]*ffr[i]**2+apol[1]*ffr[i]+apol[2]
        phs=ppol[0]*ffr[i]**2+ppol[1]*ffr[i]+ppol[2]
        c=amp*np.cos(phs)
        d=amp*np.sin(phs)
        Htmp_spec[i]=rp*c-ip*d
        Htmp_spec[i+1]=ip*c+rp*d
    Htmp=irfft(Htmp_spec)
    tr_pred.data=Htmp
    tr_left.data=tr2.data-Htmp
    return tr_pred,tr_left

def Plot_Trace(tr_list,labels=[],eq_labels=[],title=[],outfile='test.ps'):
    plt.figure(figsize=(7,9))
    ntr=len(tr_list)
    fac=[1e+6,1e+3,1e+3,1e+3,1e+6,1e+6,1,1e+6,1e+6]
    for itr,tr in enumerate(tr_list,1):
        tt=tr.times()
        tc=tr.copy()
        tc.filter('bandpass',freqmin=0.01,freqmax=0.05)
        ax=plt.subplot(ntr,1,itr)
        plt.plot(tt,tc.data*fac[itr-1]);
        ax.ticklabel_format(style='plain')
        if itr < len(tr_list):
            ax.tick_params(labelbottom=False)
#        if (itr in [1,5,6,8,9]):
#            plt.ylim((-0.00003,0.00003))
        if (len(labels)>0):
            plt.ylabel(labels[itr-1])
        if (title and itr == 1):
            plt.title(title)
        if (itr in [1,6,9]):
            ymax=np.max(tc.data)*fac[itr-1]    
            for x in eq_labels:
                plt.plot(x,ymax,'rv')
            
    plt.savefig(outfile,orientation='landscape')
    # plt.show()
    # plt.close()
    
    # tr=tr_list[-1]
    # plt.plot(tr.times(),tr.data);plt.ylim((-0.0001,0.0001))
    # plt.show()

