#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 11:40:08 2019
Estimate vertical profiles of horizontal wind from corrected data (*_l1.nc) 
measured during the CROSSINN measurement campaign. During the campaign regular  
conical PPI scans an elevation angle of 70 deg (PPI70) were performed for the SL_88 
and SLXR_142 lidars. For each PPI70 scan, one vertical profile is estimated. 
These wind profiles are collected for one day and stored into daily .nc files.
Additionally, quicklooks of the retrieved data are created. 
Used directories:
    path_lidar_in           - input directory of *_l1.nc files; structure: path_lidar_in/yyyy/yyyymm/yyyymmdd/*_l1.nc
    path_lidar_out          - output directory of retrieved wind profiles: path_lidar_out/yyyy/yyyymmdd/*_vad.nc
    path_lidar_plot_out     - output directory of figures depicting the produced data: path_lidar_plot_out/*_vad.png

Corrected .nc data files are created with crossinn_netCDF_l1.py

To run this code, the modules calc_vad.py and plot_vad.py are required. 
They can be found in the GITHub repositories:
    marenha/VAD_retrieval
    marenha/doppler_wind_lidar_toolbox/quicklooks

@author: maren
"""
import os,sys
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from datetime import date
import matplotlib.dates as mdates

# Import own modules
path_parent = os.path.abspath('..')
# path_2NetCDF = os.path.join(path_parent,'doppler_wind_lidar_toolbox','2NetCDF')
path_VAD = os.path.join(path_parent,'doppler_wind_lidar_toolbox','VAD_retrieval')
path_quicklooks = os.path.join(path_parent,'doppler_wind_lidar_toolbox','quicklooks')
# sys.path.append(path_2NetCDF)
sys.path.append(path_VAD)
sys.path.append(path_quicklooks)
import calc_vad as calc_vad
import plot_vad as plot_vad


#%% path parent 
path_data = os.path.join('/mnt','crossinn','Lidar','acinn_data_publish')

#%% lidar information and define time period
lidars_str = ['SL_88','SLXR_142']
start_str = '20190614'
end_str = '20191010'

dtn = 24*60*60 #second per day
for lidar_str in lidars_str: #loop lidars
    
    # specific paramters for the lidars can be defined here
    if lidar_str == 'SLXR_142':
        lidar_id = 'slxr142'
        snr_threshold = -24
    elif lidar_str == 'SL_88':
        lidar_id = 'sl88'
        snr_threshold = -24
    
    # paths
    path_lidar_in = os.path.join(path_data,'%s_data' %lidar_str,'data_corrected')
    path_lidar_out = os.path.join(path_data,'%s_data' %lidar_str,'data_products','VAD')
    path_lidar_plot_out = os.path.join(path_data,'%s_quicklooks' %lidar_str,'VAD')
    
    # time array containing each day
    start_num,end_num = mdates.datestr2num(start_str),mdates.datestr2num(end_str)
    time_num_array = np.arange(start_num,end_num+1,1)
    for date_num in time_num_array: #loop days
        date_str=mdates.num2datestr(date_num,'%Y%m%d')
        
        path_date = os.path.join(path_lidar_in,date_str[0:4],date_str[0:6],date_str)
        
        # list of all files in directory
        files_all=os.listdir(path_date)  
        # Only "User" files contain PPI scans
        files_user=[file for file in files_all if 'User' in file]
        
        # if no User file is available, no PPI scan was performed for this day
        if len(files_user) == 0: continue
        
        #%% Find files containing VAD PPI scans
        '''
        Find VAD PPI scans based on file size. This method only works for the 
        CROSSINN dataset and has to be adapted for other scan scenarios
        '''
        files_ppi_vad=[]
        for file_temp in files_user:
            statinfo = os.stat(os.path.join(path_date,file_temp))
            size_temp=statinfo.st_size
            if size_temp<1e6:
                files_ppi_vad.append(file_temp)
        
        # sort VAD PPI scans by time
        files_ppi_vad_sorted=sorted(files_ppi_vad)
        
        #%% collect all retrieved data in lists 
        '''
        Three types of vertical profiles are estimated :
            - Estimation of (u,v) assuming vertical veclocity w = 0 for filtered data: v_list,u_list,wd_list,ws_list
            - Estimation of (u,v) assuming vertical veclocity w = 0 for unfiltered data: v_nf_list,u_nf_list,wd_nf_list,ws_nf_list
            - Estimation of (u,v,w) by solving overdetermined system of linear equations: v_lin_list,u_lin_list,w_lin_list,wd_lin_list,ws_lin_list ,rv_fluc_list
        '''
        dn_list=[] # dn - mdates.datenum
        v_list,u_list,wd_list,ws_list = [],[],[],[] # (u,v) - 2D wind vector, wd - wind direction, ws - wind speed
        v_nf_list,u_nf_list,wd_nf_list,ws_nf_list = [],[],[],[] # nf - no filter
        v_lin_list,u_lin_list,w_lin_list,wd_lin_list,ws_lin_list = [],[],[],[],[]
        rv_fluc_list = [] # deviation of radial velocity from homogenous wind (u,v,w)
        snr_list = [] # averaged signal-to-noise ratio for each range gate
        an_list = [] # number of involved azimuth angles for each range gate
        for file_temp in files_ppi_vad_sorted: # loop PPI70 files
    
            #%% Import corrected data and extract necessary variables
            with xr.open_dataset(os.path.join(path_date,file_temp),decode_times=False) as data_temp:
            
                dn_temp=data_temp.datenum_time.values
                rv_temp=data_temp.radial_velocity.values
                
                az_temp=data_temp.azimuth.values
                el_temp=data_temp.elevation.values
                
                el_deg = np.unique(np.round(el_temp))
                
                # code online works for PPI scans at a constant elevation angle
                if len(el_deg) > 1 : continue
                
                
                gz_temp=data_temp.gate_centers.values*np.sin(np.deg2rad(el_deg))
                snr_temp=data_temp.intensity.values-1
                snr_db_temp=10*np.log10(snr_temp)
                
                lat=data_temp.lat
                lon=data_temp.lon
                alt=data_temp.alt
                
                range_gate_length=data_temp.range_gate_length
                
            # TODO: User5_142_20190801_113000_l1.nc had only 0 as azimuth angle
            if np.unique(az_temp).size < 5:
                continue
              
            #%% Apply filter to measured data 
            '''
            The data is filtered using a SNR threshold; the lowest range gates 
            contian wrong data due to laser reflection at the lidar lense
            '''
            snr_mean_temp=10*np.log10(np.mean(snr_temp,axis=1))
            rv_temp[:3,:]=np.nan
            rv_temp_nf=np.copy(rv_temp)
            rv_temp[np.isnan(snr_db_temp)]=np.nan
            rv_temp[snr_db_temp<snr_threshold]=np.nan
            rn_temp=el_temp.size
            rv_nan_sum=np.sum(np.isnan(rv_temp),axis=1)/rn_temp
    
            # at least 80 % of the data has to be available
            nan_th=.2
            rv_temp[rv_nan_sum>nan_th]=np.nan
            
            #%% loop through range gate heights --> apply VAD algorithms for each height seperately
            # create empty profiles for each variable
            gn=rv_temp.shape[0] #number of range gates
            u_lin_temp,v_lin_temp,w_lin_temp = np.full(gn,np.nan),np.full(gn,np.nan),np.full(gn,np.nan)
            ws_lin_temp,wd_lin_temp = np.full(gn,np.nan),np.full(gn,np.nan)
            rv_fluc_temp = np.full(gn,np.nan)
            ws_temp,wd_temp,u_temp,v_temp = np.full(gn,np.nan),np.full(gn,np.nan),np.full(gn,np.nan),np.full(gn,np.nan)
            ws_nf_temp,wd_nf_temp,u_nf_temp,v_nf_temp = np.full(gn,np.nan),np.full(gn,np.nan),np.full(gn,np.nan),np.full(gn,np.nan)
            
            
            '''
            since for csm scans the start azimuth angle is stored, the azimuth 
            angle has to be corrected to the centre of the angle sector
            '''
            az_delta = np.median(np.diff(az_temp))
            az_temp = az_temp + az_delta/2
            az_rad_temp,el_rad_temp=np.deg2rad(az_temp),np.deg2rad(el_temp)
            
            for gi in range(gn): # loop range gates
                rv_tempp = rv_temp[gi,:]
                ind_nan = np.where(~np.isnan(rv_tempp))[0]
    
                rn=len(ind_nan)
                if rn/rv_tempp.size<.8: continue
            
                #estimate three versions of horizontal wind 
                u_lin_temp[gi],v_lin_temp[gi],w_lin_temp[gi],ws_lin_temp[gi],wd_lin_temp[gi],rv_fluc_temp[gi] = calc_vad.calc_vad_3d(rv_tempp[ind_nan],el_rad_temp[ind_nan],az_rad_temp[ind_nan])
                ws_temp[gi],wd_temp[gi],u_temp[gi],v_temp[gi] = calc_vad.calc_vad_2d(rv_tempp[ind_nan],az_temp[ind_nan],el_temp[ind_nan])
                ws_nf_temp[gi],wd_nf_temp[gi],u_nf_temp[gi],v_nf_temp[gi] = calc_vad.calc_vad_2d(rv_temp_nf[gi,:],az_temp,el_temp)
            
            #%% collect data in list 
            an_list.append(el_temp.size)
            snr_list.append(snr_mean_temp)

            dn_list.append(dn_temp[0])
            v_list.append(v_temp)
            u_list.append(u_temp)
            wd_list.append(wd_temp)
            ws_list.append(ws_temp)
            
            v_lin_list.append(v_lin_temp)
            u_lin_list.append(u_lin_temp)
            w_lin_list.append(w_lin_temp)
            wd_lin_list.append(wd_lin_temp)
            ws_lin_list.append(ws_lin_temp)
            
            rv_fluc_list.append(rv_fluc_temp)
            
            v_nf_list.append(v_nf_temp)
            u_nf_list.append(u_nf_temp)
            wd_nf_list.append(wd_nf_temp)
            ws_nf_list.append(ws_nf_temp)
        
            az_temp[az_temp>360] = az_temp[az_temp>360]-360
        
        #%% create matrices out of the lists
        dn = np.array(dn_list)
        u = np.vstack(u_list).T
        v = np.vstack(v_list).T
        ws = np.vstack(ws_list).T
        wd = np.vstack(wd_list).T
        
        u_nf = np.vstack(u_nf_list).T
        v_nf = np.vstack(v_nf_list).T
        
        u_lin = np.vstack(u_lin_list).T
        v_lin = np.vstack(v_lin_list).T
        w_lin = np.vstack(w_lin_list).T
        ws_lin = np.vstack(ws_lin_list).T
        wd_lin = np.vstack(wd_lin_list).T
        rv_fluc = np.vstack(rv_fluc_list).T
        
        snr = np.vstack(snr_list).T
    
        
        an = np.median(an_list)
        gn = gz_temp.size
        
        datenum_array = dn
        tn = datenum_array.shape[0]
            
        #%% Create netCDF file containing horizontal wind 
        # TODO: add this part as module to marenha/doppler_wind_lidar_toolbox/2NetCDF
        
        path_out = os.path.join(path_lidar_out,date_str[0:4],date_str[0:6])
        # create output directory if non-existant
        if not os.path.exists(path_out): os.makedirs(path_out)  
            
        file_name='%s_%s_vad.nc' %(lidar_id,mdates.num2datestr(date_num,'%Y%m%d'))
        file_path=os.path.join(path_out,file_name)
    
        # overrite already existing file
        if os.path.isfile(file_path): os.remove(file_path)
        
        dataset_temp=Dataset(file_path,'w',format ='NETCDF4')
        # define dimensions
        dataset_temp.createDimension('NUMBER_OF_GATES',gn)
        dataset_temp.createDimension('NUMBER_OF_SCANS',tn)
        dataset_temp.createDimension('STATION_KEY',1)
        #    dataset_temp.createDimension('TEXT',1)
        
        # Metadata
        dataset_temp.description = "Profiles of horizontal wind speed and direction"
        dataset_temp.institution = "Department of Atmospheric and Cryospheric sciences (ACINN), University of Innsbruck, AUSTRIA",
        dataset_temp.contact = "Alexander Gohm (alexander.gohm@uibk.ac.at),Maren Haid, Lukas Lehner"
        dataset_temp.range_gate_length = range_gate_length
        dataset_temp.system_id = lidar_id
        dataset_temp.history = 'File created %s by M. Haid' %date.today().strftime('%d %b %Y')
        dataset_temp.lat = lat
        dataset_temp.lon = lon
        dataset_temp.alt = alt
        dataset_temp.snr_threshold = '%.2f dB' %snr_threshold
        dataset_temp.location_information = 'SLX142 located in trailer next to i-Box station in Kolsass during CROSSINN field campaign. Location coordinates are taken from tiris map'
        
        elevation = dataset_temp.createVariable('elevation',np.int64, ('STATION_KEY'))
        elevation.units = 'degrees'
        elevation.long_name = 'elevation'
        elevation.description = 'elevation of rays (number of rays specified in variable: rays) for VAD scan'
        elevation[:] = el_deg
        
        rays = dataset_temp.createVariable('rays',np.int64, ('STATION_KEY'))
        rays.units = 'unitless'
        rays.long_name = 'number of rays'
        rays.description = 'number of rays used to calculate mean wind profile (VAD algorithm) within the interval'
        rays[:] = an
        
        datenum = dataset_temp.createVariable('datenum',np.float64, ('NUMBER_OF_SCANS'))
        datenum.units = 'Number of days from January 0, 0000'
        datenum.long_name = 'start time of each conical scan'
        datenum.description = 'datenum (matlab) timestamp'
        datenum[:] = datenum_array
        
        unixtime = (datenum_array-mdates.datestr2num('19700101'))*(24*60*60)
        
        time = dataset_temp.createVariable('time',np.int64, ('NUMBER_OF_SCANS'))
        time.units = 'Seconds since 01-01-1970 00:00:00'
        time.long_name = 'start time of each conical scan'
        time.description = 'UNIX timestamp'
        time[:] = unixtime
        
        ff = dataset_temp.createVariable('ff', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        ff.units  ='m s-1'
        ff.long_name = 'mean horizontal wind speed'
        ff.description = 'wind speed filtered with snr threshold'
        ff[:,:] = ws_lin
        
        dd = dataset_temp.createVariable('dd', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        dd.units = 'degrees'
        dd.long_name = 'mean horizontal wind direction'
        dd.description = 'wind direction filtered with snr threshold'
        dd[:,:] = wd_lin
        
        ucomp = dataset_temp.createVariable('ucomp', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        ucomp.units = 'm s-1'
        ucomp.long_name = 'u component of horizontal wind vector'
        ucomp.description = 'u component filtered with snr threshold'
        ucomp[:,:] = u_lin
        
        vcomp = dataset_temp.createVariable('vcomp', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        vcomp.units = 'm s-1'
        vcomp.long_name = 'v component of horizontal wind vector'
        vcomp.description = 'v component filtered with snr threshold of'
        vcomp[:,:] = v_lin
        
        wcomp = dataset_temp.createVariable('wcomp', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        wcomp.units = 'm s-1'
        wcomp.long_name = 'vertical velocity component'
        wcomp.description = 'v component filtered with snr threshold'
        wcomp[:,:] = w_lin
        
        vr_fluc_var = dataset_temp.createVariable('vr_fluc_var', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        vr_fluc_var.units = 'm2 s-2'
        vr_fluc_var.long_name = 'variance of radial velocity fluctuations'
        vr_fluc_var.description = 'variance of radial velocity variations around mean (u,v,w) retrieval'
        vr_fluc_var[:,:] = rv_fluc
        
        ucomp_unfiltered = dataset_temp.createVariable('ucomp_unfiltered', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        ucomp_unfiltered.units = 'm s-1'
        ucomp_unfiltered.long_name = 'u component of horizontal wind vector'
        ucomp_unfiltered.description = 'non filtered u component'
        ucomp_unfiltered[:,:] = u_nf
        
        vcomp_unfiltered = dataset_temp.createVariable('vcomp_unfiltered', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        vcomp_unfiltered.units = 'm s-1'
        vcomp_unfiltered.long_name = 'v component of horizontal wind vector'
        vcomp_unfiltered.description = 'non filtered v component'
        vcomp_unfiltered[:,:] = v_nf
        
        snr_var = dataset_temp.createVariable('snr', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'),)
        snr_var.units = 'unitless'
        snr_var.long_name = 'signal to noise ratio (SNR)'
        snr_var.description = 'averaged profiles of snr'
        snr_var[:,:] = snr
        
        height = dataset_temp.createVariable('height',np.float32,('NUMBER_OF_GATES'))
        height.units = 'm'
        height.long_name = 'heigth of range gate centers'
        height.description = 'height of range gate centers above ground: gate_centers = (range_gate + 0.5) * range_gate_length * sin(elevation)'
        height[:] = gz_temp
            
        dataset_temp.close()
        
        #%% quicklook for created vad netCDF file
        z_ref = 546 #height of ground level at i-Box Kolsass station
        location = 'Kolsass, i-Box station'
        plot_vad.plot_VAD_day(file_path,path_lidar_plot_out,lidar_str,date_str,z_ref,location)


