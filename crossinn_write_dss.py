#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 10:08:33 2018
Write daily scan schedule (dss) for CROSSINN scan scenario and create involved
scan pattern

Here, the module write_scan_file.py is required. 
It can be found in the GITHub repositories:
    marenha/doppler_wind_lidar_toolbox/write_scan_file
@author: maren
"""
import sys
import os
# import numpy as np
import pandas as pd
import datetime as dt

# Import own modules
path_parent = os.path.abspath('..')
path_scan_files = os.path.join(path_parent,'doppler_wind_lidar_toolbox','write_scan_file')
sys.path.append(path_scan_files)
import write_scan_file as wsf

class scan_pattern():
    def __init__(self,sd,rays_avg,name):
        self.sd=sd #scan duration in minutes
        self.name=name #name of scan file
        self.rays_avg=rays_avg  #refers to average configuarion in SL software
        self.start_delay=0 #in case delay is necessary to syncronize single scans

#%% name of the used lidar system and of the scan scenario
# SL_88, SLXR_142
lidar_id='SLXR_142'
scan_shedule_name='scenario2'


#%% StreaLine Lidars parameters 
focus=7     #sl xr does not have a focus any longer
wait=0 
if lidar_id=='SL_88':
    pulse_frequency=15000
    bearing=1.7
elif lidar_id=='SLXR_142': 
    pulse_frequency=10000
    bearing=18.4  #bearing during the test campaing in autumn 2018 SLXR_142:194.2

#%% define scan pattern which should be performed in the dss
stare=scan_pattern(28,1,'stare')
ppi70=scan_pattern(2,1,'ppi70')
#rhi_ppi=scan_pattern(28,1,'rhi_ppi')
ppi_el=scan_pattern(28,1,'rhi_ppi')

#%% Create scan pattern and return scan file name
ppi70.scan_file=wsf.write_ppi(lidar_id,0,360,el=70,s=3,c=1,bearing=bearing)
#rhi_ppi.scan_file=wsf.write_ppi_rhi(lidar_id,4,65,d=28,s=3,bearing=bearing)
ppi_el.scan_file=wsf.write_ppi_el(lidar_id,0,360,el=[4,7],d=28,s=3,bearing=bearing)
#ppi_el_change.scan_file=wsf.write_ppi_el(lidar_id,0,360,el=[4,7],n=6,s=3,bearing=bearing)
stare.scan_file='stare'

#%% here, scan starts can be tuned; e.g., in case of coordinated scans
ppi70.start_delay=0
#rhi_ppi.start_delay=0
stare.start_delay=0
ppi_el.start_delay=0


#%% order in the list corresponds order in scan schedule
#scan_shedule=[stare,ppi70,rhi_ppi,ppi70]
scan_shedule=[stare,ppi70,ppi_el,ppi70]

#%% write scan scedule
time_delta=[scan.sd for scan in scan_shedule]
date_array=pd.date_range('20170101 00:00','20170102 00:00',freq=('%imin' % sum(time_delta)))
#date_array=pd.date_range('00:00','24:00',freq=('%imin' % sum(time_delta)))
scans_n=len(scan_shedule)

text_file=open(os.path.join(lidar_id,'%s.dss' % (scan_shedule_name)),'w')

for date in date_array:
    date_temp=date
    for si in range(scans_n):
        scan=scan_shedule[si]
        scan_name=scan.scan_file
        if date_temp>=date_array[-1]:
            break
        if (scan_name!='stare'):
            date_write=date_temp+dt.timedelta(0,scan.start_delay)
            if date_write<date_array[0]:
                rays_av_temp=(scan.rays_avg*pulse_frequency)/1000
                str_end='%s\t%s\t%s\t%s\t%s\r\n' % (date_write.strftime('%H%M%S'),scan_name,int(rays_av_temp),'C',focus)
                continue                
            if scan_name[0:3]=='csm':
                rays_av_temp=(scan.rays_avg*pulse_frequency)/1000
                text_file.write('%s\t%s\t%s\t%s\t%s\r\n' % (date_write.strftime('%H%M%S'),scan_name,int(rays_av_temp),'C',focus))
            elif scan_name[0:2]=='ss':
                text_file.write('%s\t%s\t%s\t%s\t%s\r\n' % (date_write.strftime('%H%M%S'),scan_name,scan.rays_av,'S',focus))
        date_temp=date_temp+dt.timedelta(minutes=scan.sd)
#text_file.write(str_end)
text_file.close()
sys.exit()




