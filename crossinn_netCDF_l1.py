#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 09:33:34 2019
Create level1 (l1) netCDF (.nc) files from level0 (l0) netCDF files collected 
during the CROSSINN campaign (https://www.imk-tro.kit.edu/english/844_8306.php) 
    - l1 data is corrected for a misalignment to north of the lidars scanner (bearing)
    - l1 data is corrected for wrong assignment of the range gate centres (gc_offset). 
    - l1 data give information about lidar location (latitude, longitude, altitude)
    - filesname of l1 files matches the l0 .nc filenames
    - l0 data was converted from .hpl files into netcdf using the module 
        hpl2NetCDF.py (GitHub: marenha/doppler_wind_lidar_toolbox/2NetCDF)

SL88:
    - bearing correction (-1.7 deg)
    - gate_centers = gate_centers - 0 m (no correction needed)
    
SLXR142:
    - bearing correction (-18.4 deg)
    - gate_centers = gate_centers - 18 m
    

Bering angles were determined by performing hart target scans of free standing 
church towers and power poles; The SLXR142 lidar was installed in a trailer, 
therefore no better alignment to north was achieved in the field.

Used directories:
    path_l0     - input directory of *_l0.nc files; structure: path_l0/[lidar_str]/netCDF/yyyy/yyyymm/yyyymmdd/*_l0.nc
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

#%% loop through lidars and campaign period (installation period)
start_str = '20190614'
end_str = '20191010'
lidars_str = ['SLXR_142','SL_88']

date_num_start,date_num_end = mdates.datestr2num(start_str),mdates.datestr2num(end_str)
dates_num=np.arange(date_num_start,date_num_end+1,1)

for lidar_str in lidars_str: #loop lidars

    # site-specific information about the Doppler wind lidars
    if lidar_str=='SL_88':
        bearing = -1.7
        gc_offset = 0
        #height MSL of ground level + height of scanner AGL
        lidar_lat, lidar_lon, lidar_alt = 47.305205, 11.622219,  546+.8
        
    elif lidar_str=='SLXR_142':
        bearing = -18.4
        gc_offset = -18
        #height MSL of ground level + height of scanner AGL
        lidar_lat, lidar_lon, lidar_alt = 47.305177, 11.622245,  546+2.5 
    else:
        ('%s was not part of CROSSINN' %lidar_str)
        continue

    for date_num in dates_num: #loop days
        date_str=mdates.num2datestr(date_num,'%Y%m%d')
        
        #directory structure corresponds to structure created by the StreamLine software
        path_l0_lidar_date=os.path.join(path_l0,lidar_str,'netCDF',date_str[0:4],date_str[0:6],date_str)
        path_l1_lidar_date=os.path.join(path_l1,'%s_data' %lidar_str,'data_corrected',date_str[0:4],date_str[0:6],date_str)

        # check if input data folder exist
        if not os.path.exists(path_l0_lidar_date): continue
    
        # create output folder 
        if not os.path.exists(path_l1_lidar_date): os.makedirs(path_l1_lidar_date)

        # list of all data files in input data folder
        files_all=sorted(os.listdir(path_l0_lidar_date))

        for file in files_all: # loop data files
            #%% Create .nc file
            ds_temp=xr.open_dataset(os.path.join(path_l0_lidar_date,file))
            
            ds_temp.azimuth.values = ds_temp.azimuth.values + bearing
            ds_temp.gate_centers.values = ds_temp.gate_centers.values + gc_offset
            ds_temp.azimuth.attrs['comment'] = 'corrected azimuth angle (az_corrected=az_measured + (%.2f)) --> az_corrected = 0 deg points to geographical North' %bearing
            ds_temp.gate_centers.attrs['comment'] = 'corrected gate center (gate_centers_corrected = gate_centers_measured + (%i) m --> gate_centers_corrected represent distance to lidar' %gc_offset
            ds_temp.attrs['lat'] = lidar_lat
            ds_temp.attrs['lon'] = lidar_lon
            ds_temp.attrs['alt'] = lidar_alt
            ds_temp.attrs['description'] = 'corrected data of Halo Photonics Streamline, corrected variables: azimuth, gate_centers'
   
            
            '''
            Sometimes timesteps of the next day are stored in a data file which 
            means that that the decimal_time starts with zero again.
            This has to be considered when calculating dn_time_temp
            '''
            #use filename to estimate datenum
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