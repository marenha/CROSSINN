#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 09:33:34 2019
Create level1 (l1) netCDF files for data collected during the CROSSINN campaign 
(https://www.imk-tro.kit.edu/english/844_8306.php) out of level0 (l0) netCDF files
    - l1 data is corrected for not perfect alignment of the lidars to north (bearing)
    - l1 data is corrected for wrong assignment of the range gate centre (gc_offset). 
    - l1 data give infromation about lidar location (latitude, longitude, altitude)
    - filesname of l1 files matches l0/hpl filenames
    - l0 data was converted from .hpl files into netcdf using the module 
        hpl2NetCDF.py (GitHub: marenha/doppler_wind_lidar_toolbox/2NetCDF)

SL88:
    - bearing correction (-1.7 deg)
    - gate_centers= gate_centers - 0 m (no correction needed)
SLXR142:
    - bearing correction (-18.4 deg)
    - gate_centers= gate_centers - 18 m
    

Bering angles were determined by performing hart target scans of free standing 
church towers and power poles; The SLXR142 lidar was installed in a trailer, 
therefore no better alignment to north could be achieved in the field;

Used directories:
    path_l0     - input directory of *_l0.nc files; structure: path_l0_in/[lidar_str]/netCDF/yyyy/yyyymm/yyyymmdd/*_l0.nc
    path_l1     - output directory of *_l1.nc files; structure: path_l1/[lidar_str]_data/data_corrected/yyyy/yyyymm/yyyymmdd/*_l1.nc

@author: maren
"""
import numpy as np
import sys,os
import xarray as xr
import matplotlib.dates as mdates

#%% paths
path_l0=os.path.join('/mnt','lidar_new','data_raw')
path_l1=os.path.join('/mnt','crossinn','Lidar','acinn_data_publish')

#%% loop through lidars and campaign period
# start_str = '20190614'
# end_str = '20191010'
# lidars_str = ['SLXR_142','SL_88']

start_str = '20190710'
end_str = '20191011'
lidars_str = ['SL_88']

date_num_start,date_num_end = mdates.datestr2num(start_str),mdates.datestr2num(end_str)
dates_num=np.arange(date_num_start,date_num_end+1,1)

for lidar_str in lidars_str: #loop through lidars
    for date_num in dates_num: #loop through days
        date_str=mdates.num2datestr(date_num,'%Y%m%d')
        
        #directory structure matches structure created by the StreamLine software
        path_l0_lidar_date=os.path.join(path_l0,lidar_str,'netCDF',date_str[0:4],date_str[0:6],date_str)
        path_l1_lidar_date=os.path.join(path_l1,'%s_data' %lidar_str,'data_corrected',date_str[0:4],date_str[0:6],date_str)
 
        if not os.path.exists(path_l1_lidar_date): os.makedirs(path_l1_lidar_date)

        if not os.path.exists(path_l0_lidar_date): continue
    
        # list of data files
        files_all=sorted(os.listdir(path_l0_lidar_date))

        for file in files_all:
            ds_temp=xr.open_dataset(os.path.join(path_l0_lidar_date,file))
        
            if lidar_str=='SL_88':
                ds_temp.azimuth.values = ds_temp.azimuth.values-1.7
                ds_temp.azimuth.attrs['comment'] = 'corrected azimuth angle (az_corrected=az_measured-1.7 deg) --> az_corrected = 0 deg points to geographical North'
                ds_temp.attrs['lat'] = 47.305205
                ds_temp.attrs['lon'] = 11.622219
                ds_temp.attrs['alt'] = 546+.8 #height MSL of ground level + height of scanner AGL
                ds_temp.attrs['description'] = 'corrected data of Halo Photonics Streamline, corrected variables: azimuth'
            if lidar_str=='SLXR_142':
                ds_temp.azimuth.values = ds_temp.azimuth.values-18.4
                ds_temp.gate_centers.values = ds_temp.gate_centers.values-18
                ds_temp.azimuth.attrs['comment'] = 'corrected azimuth angle (az_corrected=az_measured-18.4 deg) --> az_corrected = 0 deg points to geographical North'
                ds_temp.gate_centers.attrs['comment'] = 'corrected gate center (gate_centers_corrected = gate_centers_measured - 18 m --> gate_centers_corrected represent distance to lidar'
                ds_temp.attrs['lat'] = 47.305177
                ds_temp.attrs['lon'] = 11.622245
                ds_temp.attrs['alt'] = 546+2.5 #height MSL of ground level + height of scanner AGL
                ds_temp.attrs['description'] = 'corrected data of Halo Photonics Streamline, corrected variables: azimuth, gate_centers'
        
            '''
            Sometimes timesteps of the next day are assigned to the current day
            this means that the decimal_time starts with zero again
            this has to be considered when calculating dn_time_temp
            '''
            #use ffile name to estimate datenum
            dn_file = mdates.datestr2num('%s %s' %(file.split('_')[2],file.split('_')[3]))
            
            dec_time_temp=ds_temp.decimal_time.values
            dn_time_temp=date_num+dec_time_temp/24
            
            # here 12 hours are used to find affected files
            if any(np.abs((dn_time_temp-dn_file)*24) > 12):
                dn_abs=np.abs((dn_time_temp-dn_file)*24)
                #assign data to next day
                dn_time_temp[dn_abs>12] = dn_time_temp[dn_abs>12]+1

            
            dn_time_temp_var=xr.Variable(['NUMBER_OF_RAYS'],dn_time_temp,attrs={'units': 'number of days since 0001-01-01','long_name':'start time number of each ray'})
            ds_temp=ds_temp.assign(datenum_time=dn_time_temp_var)
            
            file_out=file.replace('l0','l1')
            path_out_temp=os.path.join(path_l1_lidar_date,file_out)

            if os.path.isfile(path_out_temp): os.remove(path_out_temp)

            ds_temp.to_netcdf(path_out_temp)
            ds_temp.close()