#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 6-June-2020

@author: Heather Guy, heather.guy@ncas.ac.uk
"""


# INPUTS --- EDIT THESE #

import datetime as dt
from netCDF4 import Dataset, MFDataset,num2date,date2num

infile = '/Users/heather/Desktop/Everest/'
outfile = '/Users/heather/Desktop/Everest/Everest_Storm_stats_%s.csv'%(str(dt.date.today()))

################################################################################################

# Import packages.
import sys
import numpy as np       
import pandas as pd
import warnings
import math
import glob
warnings.filterwarnings("ignore")

# Functions to process wind directions
def uv2deg(uwind,vwind):
    deg = 180 + ((np.arctan2(uwind/(np.sqrt(uwind**2 + vwind**2)), vwind/(np.sqrt(uwind**2 + vwind**2)))) * 180/np.pi) # degrees
    return deg

def deg2uv(deg):
    rad = math.radians(deg)
    uwind = -np.sin(rad)
    vwind = -np.cos(rad)
    return uwind,vwind

# Process synop codes
def synop(synop,lower_inclusive,upper_inclusive):
    tot_hours = float(len(synop))
    nhours = ((synop<=upper_inclusive) & (synop>=lower_inclusive)).sum()
    if nhours !=0: 
        percent  = nhours/tot_hours * 100
    elif nhours==0:
        percent = 0 
    if synop.isna().sum()>0:
        percent=np.nan
    return nhours,percent

# Retrieve and process infile
df = pd.read_csv(infile+'phortse.csv')#,parse_dates=[0],index_col=0,infer_datetime_format=True)

# sort dates
df.index = pd.to_datetime(df['TIMESTAMP'],dayfirst=True)
df = df.sort_index()
df = df[~df.index.duplicated()]

# Get data from other stations
basecamp = pd.read_csv(infile+'basecamp.csv')#,parse_dates=[0],index_col=0)
basecamp .index = pd.to_datetime(basecamp['TIMESTAMP'],dayfirst=True)
df = df.sort_index()
df = df[~df.index.duplicated()]

southcol = pd.read_csv(infile+'southcol.csv')#,parse_dates=[0],index_col=0)
southcol.index = pd.to_datetime(southcol['Timestamp'])
southcol = southcol.sort_index()
southcol = southcol[~southcol.index.duplicated()]

camp2 = pd.read_csv(infile+'camp2.csv')#,parse_dates=[0],index_col=0)
camp2.index = pd.to_datetime(camp2['TIMESTAMP'])
camp2 = camp2.sort_index()
camp2 = camp2[~camp2.index.duplicated()]

balcony = pd.read_csv(infile+'balcony.csv')#',parse_dates=[0],index_col=0)
balcony.index = pd.to_datetime(balcony['Timestamp'])
balcony = balcony.sort_index()
balcony = balcony[~balcony.index.duplicated()]

pyramid = pd.read_csv(infile+'Pyramid_2019.csv')#',parse_dates=[0],index_col=0)
pyramid['Timestamp'] = [pd.datetime(pyramid['YR'].iloc[i],pyramid['MN'].iloc[i],pyramid['DAY'].iloc[i],pyramid['HR'].iloc[i]) for i in range(0,len(pyramid))]
pyramid.index = pd.to_datetime(pyramid['Timestamp'])
pyramid = pyramid.sort_index()
pyramid = pyramid[~pyramid.index.duplicated()]

pheriche = pd.read_csv(infile+'Pheriche_2019.csv')#',parse_dates=[0],index_col=0)
pheriche['Timestamp'] = [pd.datetime(pheriche['YR'].iloc[i],pheriche['MN'].iloc[i],pheriche['DAY'].iloc[i],pheriche['HR'].iloc[i]) for i in range(0,len(pheriche))]
pheriche.index = pd.to_datetime(pheriche['Timestamp'])
pheriche = pheriche.sort_index()
pheriche = pheriche[~pheriche.index.duplicated()]






# Get era5 data

era_surface_f = '/Users/heather/jasmin_ncas_vol1/heather/ERA5/Phortse_surface.nc'
era_surface = Dataset(era_surface_f,'r')

era_level_fils = glob.glob('/Users/heather/jasmin_ncas_vol1/heather/ERA5/Phortse_20*.nc')
era_levels = MFDataset(era_level_fils,'r',aggdim='time')

times = [dt.datetime(1900,1,1,0) + dt.timedelta(hours=int(era_levels.variables['time'][:][i])) for i in range(0,len(era_levels.variables['time'][:]))]


# To get surface orography, see https://confluence.ecmwf.int/display/CKB/ERA5%3A+surface+elevation+and+orography
# Surface geopotential in m2 s-2
surface_geopotential = np.mean(era_surface.variables['z'][:,0,1,1])
# divide my g to get value in m
magl = surface_geopotential / 9.80665

#surface_times = [dt.datetime(1900,1,1,0) + dt.timedelta(hours=int(era_surface.variables['time'][:][i])) for i in range(0,len(era_surface.variables['time'][:]))]

era_cols = ['era5_500hPa_T', 'era5_400hPa_T', 'era5_300hPa_T', 'era5_500hPa_rh', 'era5_400hPa_rh', 'era5_300hPa_rh', 'era5_500hPa_sp', 'era5_400hPa_sp', 'era5_300hPa_sp', 'era5_500hPa_ws', 'era5_400hPa_ws', 'era5_300hPa_ws', 'era5_500hPa_wd', 'era5_400hPa_wd', 'era5_300hPa_wd', 'era5_500hPa_vv', 'era5_total_column_wv', 'era5_total_precipitation', 'era5_zero_degree_level']
era_pdf = pd.DataFrame(columns=era_cols,index=times)

# Level indices: 0: 300hpa, 1: 400 hPa, 2: 500 hPa
# Grid is averaged over the following box:
# longitude: 86.75
# Latitude: 27.75
# variable format: time, level, lat, lon

era_pdf['era5_500hPa_T'] = era_levels.variables['t'][:,2,1,1] - 273.15 
era_pdf['era5_400hPa_T'] = era_levels.variables['t'][:,1,1,1] - 273.15 
era_pdf['era5_300hPa_T'] = era_levels.variables['t'][:,0,1,1] - 273.15 
era_pdf['era5_500hPa_rh'] = era_levels.variables['r'][:,2,1,1]
era_pdf['era5_400hPa_rh'] = era_levels.variables['r'][:,1,1,1]
era_pdf['era5_300hPa_rh'] = era_levels.variables['r'][:,0,1,1]
era_pdf['era5_500hPa_sp'] = era_levels.variables['q'][:,2,1,1]
era_pdf['era5_400hPa_sp'] = era_levels.variables['q'][:,1,1,1]
era_pdf['era5_300hPa_sp'] = era_levels.variables['q'][:,0,1,1]
era_pdf['era5_500hPa_vv'] = era_levels.variables['w'][:,2,1,1]

# varible format: time, expver, latitude, longitude
# Has 9 additional hours compared to levels
# Expver is 1 or 5. These are different versions because accumulated variables are derived from forecast fields
# exper 5 is near real time gap filling preliminary data
# exper 1 is normal era5
# for this range should just be able to use expver=1 (index=0)

era_pdf['era5_total_column_wv'] = era_surface.variables['tcwv'][:-57,0,1,1]
era_pdf['era5_total_precipitation'] = era_surface.variables['tp'][:-57,0,1,1]

# notes for zero deg level
# The height above the Earth's surface where the temperature passes from positive to negative values, 
# corresponding to the top of a warm layer, at the specified time. This parameter can be used to help forecast snow.
# If more than one warm layer is encountered, then the zero degree level corresponds to the top of the second atmospheric layer.
# This parameter is set to zero when the temperature in the whole atmosphere is below 0℃.

era_pdf['era5_zero_degree_level'] = era_surface.variables['deg0l'][:-57,0,1,1]

# Sort era5 winds
 
era_pdf['era5_500hPa_ws'] = np.sqrt(era_levels.variables['u'][:,2,1,1]**2 + era_levels.variables['v'][:,0,1,1]**2)
era_pdf['era5_400hPa_ws'] = np.sqrt(era_levels.variables['u'][:,1,1,1]**2 + era_levels.variables['v'][:,0,1,1]**2)
era_pdf['era5_300hPa_ws'] = np.sqrt(era_levels.variables['u'][:,0,1,1]**2 + era_levels.variables['v'][:,0,1,1]**2)

era_pdf['era5_500hPa_wd'] = uv2deg(era_levels.variables['u'][:,2,1,1],era_levels.variables['v'][:,0,1,1])    
era_pdf['era5_400hPa_wd'] = uv2deg(era_levels.variables['u'][:,1,1,1],era_levels.variables['v'][:,0,1,1])    
era_pdf['era5_300hPa_wd'] = uv2deg(era_levels.variables['u'][:,0,1,1],era_levels.variables['v'][:,0,1,1])    

era_pdf['era5_500hPa_vv'] = era_levels.variables['w'][:,2,1,1]


# Group individual storms
df_precip = df[df.Precip_Tot != 0]                                              # Only times with precip
df_precip['temp1'] =df_precip.index                                             # Seperate times 
df_precip['Storm'] = df_precip['temp1'].diff() > pd.Timedelta(hours=6)          # True where difference greater than x hours.
df_precip['Storm'] = df_precip['Storm'].apply(lambda x: 1 if x else 0).cumsum() # Storm number = cumulative sum of trues. 
del df_precip['temp1']                                                          # delete time date
storms = df_precip.groupby('Storm')                                             # group by storms
nstorms = len(storms)                                                           # Total number of storms

# Open write file
f = open(outfile,'w')
columns = 'StormID,Start_year,Start_month,Start_day,Start_hour,Maturation_year,Maturation_month,Maturation_day,Maturation_hour,End_year,End_month,End_day,End_hour,Duration,MeanT,MeanRH,TotalPrecip,MeanWindS,MeanWindD,MaxWindS,SnowAccumulation,synop_nhours_60,synop_percent_60,synop_nhours_61_63,synop_percent_61_63,synop_nhours_67_68,synop_percent_67_68,synop_nhours_71,synop_percent_71,synop_nhours_72,synop_percent_72,synop_nhours_73,synop_percent_73,synop_nhours_88_89,synop_percent_88_89, basecamp_precip, basecamp_Tmean,camp2_Tmean,southcol_Tmean,balcony_Tmean, basecamp_rh,camp2_rh,southcol_rh,balcony_rh,camp2_ws,southcol_ws,balcony_ws, camp2_wd,southcol_wd,balcony_wd,pyramid_Tmean,pyramid_rh,pyramid_ws,pyramid_precip,pheriche_Tmean,pheriche_rh,pheriche_ws,pheriche_precip, era5_500hPa_T, era5_400hPa_T, era5_300hPa_T, era5_500hPa_rh, era5_400hPa_rh, era5_300hPa_rh, era5_500hPa_sp, era5_400hPa_sp, era5_300hPa_sp, era5_500hPa_ws, era5_400hPa_ws, era5_300hPa_ws, era5_500hPa_wd, era5_400hPa_wd, era5_300hPa_wd, era5_500hPa_vv, era5_total_column_wv, era5_total_precipitation, era5_zero_degree_level \n'
units = '#,,,,,,,,,,,,,,C,%,Pecip_Tot,WS,degrees,WS_Max,SnDep,hours,%,hours,%,hours,percent,hours,percent,hours,percent,hours,percent,hours,percent, PRECIP, C,C,C,C, %,%,%,%, ,,,degrees,degrees,degrees,C,%,,,C,%,,, C, C, C, %, %, %, kg/kg, kg/kg, kg/kg, m/s,  m/s,  m/s, degrees, degrees, degrees, Pa/s, kg/m3, m, m above surface\n'
f.write(columns) # Write headers. 
f.write(units) # Write headers. 

# Loop through storms to calculate stats
for i in range(0,nstorms):
    storm = storms.get_group(i)
    start = storm.index[0]
    end = storm.index[-1]
    duration = (end-start).days * 24*60*60 + (end - start).seconds
    if duration == 0: 
        duration = 3600
    duration = duration/3600 # hours
    

    meanT = storm['AirTC_Avg'].astype(float).mean()
    meanRH = storm['RH'].astype(float).mean()
    totP = storm['Precip_Tot'].astype(float).sum()
    meanWS = storm['WS'].astype(float).mean()
    maxWS = storm['WS_Max'].astype(float).max() 
    # Calculate mean wind direction. 
    us=[]
    vs=[]
    for j in range(0,len(storm)):
        u,v = deg2uv(storm['WD'].astype(int)[j])
        us.append(u)
        vs.append(v) 
    umean = np.nanmean(us)
    vmean = np.nanmean(vs)
    wdmean =uv2deg(umean,vmean)
       
    snd = storm['SnDep'].astype(float)[-1] - storm['SnDep'].astype(float)[0]
        
    if snd <0:   # Only record snow depth change if positive
        snd = 0

    mat = storm['Precip_Tot'].idxmax() - dt.timedelta(hours=1)  
    
    #Number of hours with synop code <60, 61-63, 67-68, 71, 72, 73, 88-89 in separate columns
    #Percent of hours with synop code <60, 61-63, 67-68, 71, 72, 73, 88-89 in separate columns
    synop_nhours_60,synop_percent_60 = synop(storm['synop'],0,59)  
    synop_nhours_61_63,synop_percent_61_63 = synop(storm['synop'],61,63) 
    synop_nhours_67_68,synop_percent_67_68 = synop(storm['synop'],67,68) 
    synop_nhours_71,synop_percent_71 = synop(storm['synop'],71,71) 
    synop_nhours_72,synop_percent_72 = synop(storm['synop'],72,72) 
    synop_nhours_73,synop_percent_73 = synop(storm['synop'],73,73) 
    synop_nhours_88_89,synop_percent_88_89 = synop(storm['synop'],88,89) 
    
    
    # Other station stats
    try:
        basecamp_precip = basecamp[start:end]['PRECIP'].sum()
        basecamp_Tmean = basecamp[start:end]['AirTC_Avg'].mean()
        basecamp_rh = basecamp[start:end]['RH'].mean()  
    except:
        basecamp_precip = np.nan
        basecamp_Tmean = np.nan
        basecamp_rh = np.nan

    try:
        camp2_Tmean = camp2[start:end]['AirTC'].mean()
        camp2_rh = camp2[start:end]['RH'].mean()  
        camp2_ws = camp2[start:end]['WS'].mean()
        # Calculate mean wind direction. 
        us=[]
        vs=[]
        for j in range(0,len(camp2[start:end])):
            u,v = deg2uv(camp2[start:end]['WD'].astype(int)[j])
            us.append(u)
            vs.append(v) 
        umean = np.nanmean(us)
        vmean = np.nanmean(vs)
        camp2_wd = uv2deg(umean,vmean)            
    except:
        camp2_Tmean = np.nan
        camp2_rh = np.nan    
        camp2_ws = np.nan
        camp2_wd = np.nan
        
    try:
        southcol_Tmean = southcol[start:end]['T_HMP'].mean()
        southcol_rh = southcol[start:end]['RH'].mean()  
        southcol_ws = southcol[start:end]['WS_AVG'].mean()
        # Calculate mean wind direction. 
        us=[]
        vs=[]
        for j in range(0,len(southcol[start:end])):
            u,v = deg2uv(southcol[start:end]['WDIR'].astype(int)[j])
            us.append(u)
            vs.append(v) 
        umean = np.nanmean(us)
        vmean = np.nanmean(vs)
        southcol_wd = uv2deg(umean,vmean)            
    except:
        southcol_Tmean = np.nan
        southcol_rh = np.nan    
        southcol_ws = np.nan
        southcol_wd = np.nan        
     
    try:
        balcony_Tmean = balcony[start:end]['T_HMP'].astype('float').mean()
        balcony_rh = balcony[start:end]['RH'].astype('float').mean()  
        balcony_ws = balcony[start:end]['WS_AVG_1'].astype('float').mean()
        # Calculate mean wind direction. 
        us=[]
        vs=[]
        for j in range(0,len(balcony[start:end])):
            u,v = deg2uv(balcony[start:end]['WDIR_1'].astype(int)[j])
            us.append(u)
            vs.append(v) 
        umean = np.nanmean(us)
        vmean = np.nanmean(vs)
        balcony_wd = uv2deg(umean,vmean)            
    except:
        balcony_Tmean = np.nan
        balcony_rh = np.nan    
        balcony_ws = np.nan
        balcony_wd = np.nan      
        
    try:
        pyramid_Tmean = pyramid[start:end]['AirTC'].astype('float').mean()
        pyramid_rh = pyramid[start:end]['RH'].astype('float').mean()  
        pyramid_ws = pyramid[start:end]['WS_u'].astype('float').mean()
        pyramid_precip = pyramid[start:end]['Precip'].astype('float').sum()
    except:
        pyramid_Tmean = np.nan
        pyramid_rh = np.nan    
        pyramid_ws = np.nan
        pyramid_precip = np.nan
        
    try:
        pheriche_Tmean = pheriche[start:end]['AirTC'].astype('float').mean()
        pheriche_rh = pheriche[start:end]['RH'].astype('float').mean()  
        pheriche_ws = pheriche[start:end]['WS_u'].astype('float').mean()
        pheriche_precip = pheriche[start:end]['Precip'].astype('float').sum()
    except:
        pheriche_Tmean = np.nan
        pheriche_rh = np.nan    
        pheriche_ws = np.nan
        pheriche_precip = np.nan
 
                
    # ERA stats
    
    era5_500hPa_T = era_pdf['era5_500hPa_T'][start:end].mean()
    era5_400hPa_T = era_pdf['era5_400hPa_T'][start:end].mean()
    era5_300hPa_T = era_pdf['era5_300hPa_T'][start:end].mean()
    era5_500hPa_rh = era_pdf['era5_500hPa_rh'][start:end].mean()
    era5_400hPa_rh = era_pdf['era5_400hPa_rh'][start:end].mean()
    era5_300hPa_rh = era_pdf['era5_300hPa_rh'][start:end].mean()
    era5_500hPa_sp = era_pdf['era5_500hPa_sp'][start:end].mean()
    era5_400hPa_sp = era_pdf['era5_400hPa_sp'][start:end].mean()
    era5_300hPa_sp = era_pdf['era5_300hPa_sp'][start:end].mean()
    era5_500hPa_ws = era_pdf['era5_500hPa_ws'][start:end].mean()
    era5_400hPa_ws = era_pdf['era5_400hPa_ws'][start:end].mean()
    era5_300hPa_ws = era_pdf['era5_300hPa_ws'][start:end].mean()
    era5_500hPa_wd = era_pdf['era5_500hPa_wd'][start:end].mean()
    era5_400hPa_wd = era_pdf['era5_400hPa_wd'][start:end].mean()
    era5_300hPa_wd = era_pdf['era5_300hPa_wd'][start:end].mean()
    era5_500hPa_vv = era_pdf['era5_500hPa_vv'][start:end].mean()
    
    era5_total_column_wv = era_pdf['era5_total_column_wv'][start:end].mean()
    era5_total_precipitation = era_pdf['era5_total_precipitation'][start:end].mean()
    era5_zero_degree_level = era_pdf['era5_zero_degree_level'][start:end].mean()
    
    if era5_total_column_wv <= 0 :
        era5_total_column_wv=np.nan
        era5_total_precipitation=np.nan
        era5_zero_degree_level=np.nan
    
    # Write stats to outfile
    if type(mat) != pd._libs.tslibs.timestamps.Timestamp:
        
        row = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s, %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'%(i,start.year,start.month,start.day,start.hour,mat,mat,mat,mat,end.year,end.month,end.day,end.hour, duration, meanT,meanRH,totP,meanWS,wdmean,maxWS,snd,synop_nhours_60,synop_percent_60,synop_nhours_61_63,synop_percent_61_63,synop_nhours_67_68,synop_percent_67_68,synop_nhours_71,synop_percent_71,synop_nhours_72,synop_percent_72,synop_nhours_73,synop_percent_73,synop_nhours_88_89,synop_percent_88_89,basecamp_precip, basecamp_Tmean,camp2_Tmean,southcol_Tmean,balcony_Tmean, basecamp_rh,camp2_rh,southcol_rh,balcony_rh, camp2_ws,southcol_ws,balcony_ws, camp2_wd,southcol_wd,balcony_wd,pyramid_Tmean,pyramid_rh,pyramid_ws,pyramid_precip,pheriche_Tmean,pheriche_rh,pheriche_ws,pheriche_precip, era5_500hPa_T, era5_400hPa_T, era5_300hPa_T, era5_500hPa_rh, era5_400hPa_rh, era5_300hPa_rh, era5_500hPa_sp, era5_400hPa_sp, era5_300hPa_sp, era5_500hPa_ws, era5_400hPa_ws, era5_300hPa_ws, era5_500hPa_wd, era5_400hPa_wd, era5_300hPa_wd, era5_500hPa_vv, era5_total_column_wv, era5_total_precipitation, era5_zero_degree_level)
    else:
        row = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s, %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'%(i,start.year,start.month,start.day,start.hour,mat.year,mat.month,mat.day,mat.hour,end.year,end.month,end.day,end.hour, duration, meanT,meanRH,totP,meanWS,wdmean,maxWS,snd,synop_nhours_60,synop_percent_60,synop_nhours_61_63,synop_percent_61_63,synop_nhours_67_68,synop_percent_67_68,synop_nhours_71,synop_percent_71,synop_nhours_72,synop_percent_72,synop_nhours_73,synop_percent_73,synop_nhours_88_89,synop_percent_88_89,basecamp_precip, basecamp_Tmean,camp2_Tmean,southcol_Tmean,balcony_Tmean, basecamp_rh,camp2_rh,southcol_rh,balcony_rh, camp2_ws,southcol_ws,balcony_ws, camp2_wd,southcol_wd,balcony_wd,pyramid_Tmean,pyramid_rh,pyramid_ws,pyramid_precip,pheriche_Tmean,pheriche_rh,pheriche_ws,pheriche_precip, era5_500hPa_T, era5_400hPa_T, era5_300hPa_T, era5_500hPa_rh, era5_400hPa_rh, era5_300hPa_rh, era5_500hPa_sp, era5_400hPa_sp, era5_300hPa_sp, era5_500hPa_ws, era5_400hPa_ws, era5_300hPa_ws, era5_500hPa_wd, era5_400hPa_wd, era5_300hPa_wd, era5_500hPa_vv, era5_total_column_wv, era5_total_precipitation, era5_zero_degree_level)
    
    f.write(row) # Write row.
f.close()
