#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 11:40:08 2019
Estimate horizontal wind vectors for SL88 PPI70 scans from l0 data and combine
them into daily files: [lidar_str]_vad_yyyymmdd.nc

@author: maren
"""
import os,sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from datetime import date
import matplotlib.dates as mdates

# Import modules
path_parent = os.path.abspath('..')
path_2NetCDF = os.path.join(path_parent,'doppler_wind_lidar_toolbox','VAD_retrievals')
path_2NetCDF = os.path.join(path_parent,'doppler_wind_lidar_toolbox','quicklooks')
sys.path.append(path_2NetCDF)
import calc_vad as calc_vad
import plot_VAD as plot_VAD


#%% path of input level0 files and level1 output directory 
path_data = os.path.join('/mnt','crossinn','Lidar','acinn_data_publish')

#%% lidar information and define time period
lidars_str = ['SL_88','SLXR_142']
start_str='20190704'
end_str='20191010'


el_deg=70.

dtn=24*60*60
for lidar_str in lidars_str:
    
    # specific paramters for the lidars can be defined here
    if lidar_str == 'SLXR_142':
        lidar_id='sl88'
        snr_threshold=-24
    elif lidar_str == 'SL_88':
        lidar_id='slxr142'
        snr_threshold=-24
        
    path_lidar_in = os.path.join(path_data,'%s_data' %lidar_str,'data_corrected')
    path_lidar_out = os.path.join(path_data,'%s_data' %lidar_str,'data_products','VAD')
    path_lidar_plot_out = os.path.join(path_data,'%s_quicklooks' %lidar_str,'VAD')
    
    start_num,end_num=mdates.datestr2num(start_str),mdates.datestr2num(end_str)
    time_num_array=np.arange(start_num,end_num+1,1)
    for date_num in time_num_array:
        date_str=mdates.num2datestr(date_num,'%Y%m%d')
        
        '''
        coordiantes are taken from tiris map and are just a rough determination 
        same is valid for height of the instrument
        '''
        path_date = os.path.join(path_lidar_in,date_str[0:4],date_str[0:6],date_str)
        
        files_all=os.listdir(path_date)  
        files_user=[file for file in files_all if 'User' in file]
        
        if len(files_user)==0: continue
        
        
        files_ppi_vad=[]
        for file_temp in files_user:
            statinfo = os.stat(os.path.join(path_date,file_temp))
            size_temp=statinfo.st_size
            if size_temp<1e6:
                files_ppi_vad.append(file_temp)
                
        
        files_ppi_vad_dates=[mdates.datestr2num('%s %s' %(file.split('_')[2],file.split('_')[3].split('.')[0])) for file in files_ppi_vad]
        ind_sort=np.argsort(files_ppi_vad_dates)
        files_ppi_vad_dates_sorted=np.array(files_ppi_vad_dates)[ind_sort]
        files_ppi_vad_sorted=np.array(files_ppi_vad)[ind_sort]
        dn_list=[]
        v_list,u_list,wd_list,ws_list=[],[],[],[]
        v_nf_list,u_nf_list,wd_nf_list,ws_nf_list=[],[],[],[]
        v_lin_list,u_lin_list,w_lin_list=[],[],[]
        wd_lin_list,ws_lin_list=[],[]
        rv_fluc_list=[]
        snr_list=[]
        an_list=[]
        for file_temp in files_ppi_vad_sorted: # go through files for one day
    
            with xr.open_dataset(os.path.join(path_date,file_temp),decode_times=False) as data_temp:
            
                dn_temp=data_temp.datenum_time.values
                rv_temp=data_temp.radial_velocity.values
                
                az_temp=data_temp.azimuth.values
                el_temp=data_temp.elevation.values
                gz_temp=data_temp.gate_centers.values*np.sin(np.deg2rad(el_deg))
                snr_temp=data_temp.intensity.values-1
                snr_db_temp=10*np.log10(snr_temp)
                
                lat=data_temp.lat
                lon=data_temp.lon
                alt=data_temp.alt
                
                range_gate_length=data_temp.range_gate_length
              
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
            
            gn=rv_temp.shape[0] #number of range gates
            u_lin_temp,v_lin_temp,w_lin_temp=np.full(gn,np.nan),np.full(gn,np.nan),np.full(gn,np.nan)
            ws_lin_temp,wd_lin_temp=np.full(gn,np.nan),np.full(gn,np.nan)
            rv_fluc_temp=np.full(gn,np.nan)
            ws_temp,wd_temp,u_temp,v_temp=np.full(gn,np.nan),np.full(gn,np.nan),np.full(gn,np.nan),np.full(gn,np.nan)
            ws_nf_temp,wd_nf_temp,u_nf_temp,v_nf_temp=np.full(gn,np.nan),np.full(gn,np.nan),np.full(gn,np.nan),np.full(gn,np.nan)
            
            az_rad_temp,el_rad_temp=np.deg2rad(az_temp),np.deg2rad(el_temp)
            for gi in range(gn):
                rv_tempp=rv_temp[gi,:]
                ind_nan=np.where(~np.isnan(rv_tempp))[0]
    
                rn=len(ind_nan)
                if rn/rv_tempp.size<.8: continue
            
                #estimate three versions of horizontal wind 
                u_lin_temp[gi],v_lin_temp[gi],w_lin_temp[gi],ws_lin_temp[gi],wd_lin_temp[gi],rv_fluc_temp[gi] = calc_vad.calc_vad_3d(rv_tempp[ind_nan],el_rad_temp[ind_nan],az_rad_temp[ind_nan])
                ws_temp[gi],wd_temp[gi],u_temp[gi],v_temp[gi]=calc_vad.calc_vad_2d(rv_temp[ind_nan],az_temp[ind_nan],el_temp[ind_nan])
                ws_nf_temp[gi],wd_nf_temp[gi],u_nf_temp[gi],v_nf_temp[gi]=calc_vad.calc_vad(rv_temp_nf[gi,:],az_temp,el_temp)
            
            an_list.append(el_temp.size)
            snr_list.append(snr_mean_temp)
            
            # colltect each profile in list 
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
        
            az_temp[az_temp>360]=az_temp[az_temp>360]-360
        
        #create matriced our of lists
        dn=np.array(dn_list)
        u=np.vstack(u_list).T
        v=np.vstack(v_list).T
        ws=np.vstack(ws_list).T
        wd=np.vstack(wd_list).T
        
        u_nf=np.vstack(u_nf_list).T
        v_nf=np.vstack(v_nf_list).T
        
        u_lin=np.vstack(u_lin_list).T
        v_lin=np.vstack(v_lin_list).T
        w_lin=np.vstack(w_lin_list).T
        ws_lin=np.vstack(ws_lin_list).T
        wd_lin=np.vstack(wd_lin_list).T
        rv_fluc=np.vstack(rv_fluc_list).T
        
        snr=np.vstack(snr_list).T
    
        
        an=np.median(an_list)
        gn=gz_temp.size
        
        datenum_array=dn
        tn=datenum_array.shape[0]
            
        path_out=os.path.join(path_lidar_out,date_str[0:4],date_str[0:6])
        print(path_out)
        # create output directory if non-existant
        if not os.path.exists(path_out):
            os.makedirs(path_out)  
        file_name='%s_%s_vad.nc' %(lidar_id,mdates.num2datestr(date_num,'%Y%m%d'))
        file_path=os.path.join(path_out,file_name)
    
        if os.path.isfile(file_path):
            os.remove(file_path)
        
        #%% Create netCDF file containing hroizontal wind 
        
        dataset_temp=Dataset(file_path,'w',format ='NETCDF4')
        # define dimensions
        dataset_temp.createDimension('NUMBER_OF_GATES',gn)
        dataset_temp.createDimension('NUMBER_OF_SCANS',tn)
        dataset_temp.createDimension('STATION_KEY',1)
        #    dataset_temp.createDimension('TEXT',1)
        
        # Metadata
        dataset_temp.description="Profiles of horizontal wind speed and direction"
        dataset_temp.institution="Department of Atmospheric and Cryospheric sciences (ACINN), University of Innsbruck, AUSTRIA",
        dataset_temp.contact="Alexander Gohm (alexander.gohm@uibk.ac.at),Maren Haid, Lukas Lehner"
        dataset_temp.range_gate_length=range_gate_length
        dataset_temp.system_id=lidar_id
        dataset_temp.history='File created %s by M. Haid' %date.today().strftime('%d %b %Y')
        dataset_temp.lat=lat
        dataset_temp.lon=lon
        dataset_temp.alt=alt
        dataset_temp.snr_threshold='%.2f dB' %snr_threshold
        dataset_temp.location_information='SLX142 located in trailer next to i-Box station in Kolsass during CROSSINN field campaign. Location coordinates are taken from tiris map'
        
        elevation=dataset_temp.createVariable('elevation',np.int64, ('STATION_KEY'))
        elevation.units='degrees'
        elevation.long_name='elevation'
        elevation.description='elevation of rays (number of rays specified in variable: rays) for VAD scan'
        elevation[:]=el_deg
        
        rays=dataset_temp.createVariable('rays',np.int64, ('STATION_KEY'))
        rays.units='unitless'
        rays.long_name='number of rays'
        rays.description='number of rays used to calculate mean wind profile (VAD algorithm) within the interval'
        rays[:]=an
        
        datenum=dataset_temp.createVariable('datenum',np.float64, ('NUMBER_OF_SCANS'))
        datenum.units='Number of days from January 0, 0000'
        datenum.long_name='start time of each conical scan'
        datenum.description='datenum (matlab) timestamp'
        datenum[:]=datenum_array
        
        unixtime=(datenum_array-mdates.datestr2num('19700101'))*(24*60*60)
        
        time=dataset_temp.createVariable('time',np.int64, ('NUMBER_OF_SCANS'))
        time.units='Seconds since 01-01-1970 00:00:00'
        time.long_name='start time of each conical scan'
        time.description='UNIX timestamp'
        time[:]=unixtime
        
        ff=dataset_temp.createVariable('ff', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        ff.units='m s-1'
        ff.long_name='mean horizontal wind speed'
        ff.description='wind speed filtered with snr threshold'
        ff[:,:]=ws_lin
        
        dd=dataset_temp.createVariable('dd', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        dd.units='degrees'
        dd.long_name='mean horizontal wind direction'
        dd.description='wind direction filtered with snr threshold'
        dd[:,:]=wd_lin
        
        ucomp=dataset_temp.createVariable('ucomp', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        ucomp.units='m s-1'
        ucomp.long_name='u component of horizontal wind vector'
        ucomp.description='u component filtered with snr threshold'
        ucomp[:,:]=u_lin
        
        vcomp=dataset_temp.createVariable('vcomp', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        vcomp.units='m s-1'
        vcomp.long_name='v component of horizontal wind vector'
        vcomp.description='v component filtered with snr threshold of'
        vcomp[:,:]=v_lin
        
        wcomp=dataset_temp.createVariable('wcomp', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        wcomp.units='m s-1'
        wcomp.long_name='vertical velocity component'
        wcomp.description='v component filtered with snr threshold'
        wcomp[:,:]=w_lin
        
        vr_fluc_var=dataset_temp.createVariable('vr_fluc_var', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        vr_fluc_var.units='m2 s-2'
        vr_fluc_var.long_name='variance of radial velocity fluctuations'
        vr_fluc_var.description='variance of radial velocity variations around mean (u,v,w) retrieval'
        vr_fluc_var[:,:]=rv_fluc
        
        ucomp_unfiltered=dataset_temp.createVariable('ucomp_unfiltered', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        ucomp_unfiltered.units='m s-1'
        ucomp_unfiltered.long_name='u component of horizontal wind vector'
        ucomp_unfiltered.description='non filtered u component'
        ucomp_unfiltered[:,:]=u_nf
        
        vcomp_unfiltered=dataset_temp.createVariable('vcomp_unfiltered', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        vcomp_unfiltered.units='m s-1'
        vcomp_unfiltered.long_name='v component of horizontal wind vector'
        vcomp_unfiltered.description='non filtered v component'
        vcomp_unfiltered[:,:]=v_nf
        
        snr_var=dataset_temp.createVariable('snr', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'),)
        snr_var.units='unitless'
        snr_var.long_name='signal to noise ratio (SNR)'
        snr_var.description='averaged profiles of snr'
        snr_var[:,:]=snr
        
        height=dataset_temp.createVariable('height',np.float32,('NUMBER_OF_GATES'))
        height.units='m'
        height.long_name='heigth of range gate centers'
        height.description='height of range gate centers above ground: gate_centers = (range_gate + 0.5) * range_gate_length * sin(elevation)'
        height[:]=gz_temp
            
        dataset_temp.close()
        
        #%% quicklook for created vad netCDF file
        z_ref = 546 #height of ground level at i-Box Kolsass station
        location = 'Kolsass, i-Box station'
        plot_VAD.plot_VAD_day(file_path,path_lidar_plot_out,lidar_str,date_str,z_ref,location)