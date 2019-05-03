#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 11:55:30 2018

@author: zhitu
"""

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.core.event import read_events
from remove_noise_oneday import DoDay
import numpy as np

def Doit():
    from sys import argv
    
    nstas=[argv[1]]
#    nstas=['FS05B']
    print("doing nstas",nstas)    

    
    client=Client("IRIS")
    t1=UTCDateTime(year=2012,julday=1)
    t2=UTCDateTime(year=2012,julday=365)
    
    inventory = client.get_stations(network='7D',starttime=t1,endtime=t2)
    
    allstas={}
    for sta in inventory[0].stations:
        nsta=sta.code
        xlat=sta.latitude
        xlon=sta.longitude
        xdep=sta.elevation
        allstas[nsta]=(xlat,xlon,np.abs(xdep))


    cat=read_events("cmt2012.ndk")
    cat=cat.filter("magnitude >= 5.0")
    cat_large=cat.filter("magnitude >= 5.5")
    t1_all=[]
    for item in cat_large:
       t1_all.append(item.origins[0].time)
    print("number of evenst:",len(t1_all))
    # t1_all.append(cat_large[100].origins[0].time)
    # print(t1_all) 
    
    for nsta in nstas:
        outfile='./'+nsta+'_result_new.txt'
#        for t1 in t1_all[100:101]:
        for t1 in t1_all:
            fobj=open(outfile,'a')
            t2=t1+86400
            print("doing",nsta,t1)
            ierr,angle,coh,adm,phs,coh2,adm2,phs2=DoDay(nsta,t1,t2,cat,allstas,iopt=2)
            fobj.write("%3d %3d %10.3f %10.3f %12.5g %12.5g %10.3f %12.5g %12.5g\n" \
                       % (ierr,t1.julday,np.rad2deg(angle),coh,adm,phs,coh2,adm2,phs2))
            fobj.close()


if __name__ == "__main__":
    Doit()
