#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 21:25:39 2020

@author: heather
"""

import pandas as pd
import numpy as np
import datetime as dt
import glob as glob
import os
import warnings
warnings.filterwarnings("ignore")

run_length = 72 # hours
sum_lat = 27.850
sum_lon = 86.750
col_headers = ['time','meanLon','meanLat','meanZ','meanTopo','meanPBL','meanTropo','meanPv','rmsHBefore','rmsHAfter','rmsVBefore','rmsVAfter','pblFract','pv2Fract','tropoFract']
nclusters=5
mean_widths = [5,8,9,9,8,8,8,8,8,8,8,8,8,6,6,6]
cluster_widths = [8,8,7,6,8]
widths = mean_widths + cluster_widths * nclusters




# Get event filesbase_dir='/Users/heather/Desktop/Everest/everest-flexpart/'
#out_dir='/Users/heather/Desktop/Everest/everest-flexpart/endpoints/'
base_dir='/nobackup/eehgu/everest-flexpart/events_masl/'
out_dir='/nobackup/eehgu/everest-flexpart/events_masl/endpoints/'

#release_times = list(set([x[-37:-23] for x in glob.glob(base_dir +'20*')]))
release_times = list(set([x[-38:-24] for x in glob.glob(base_dir +'20*')]))
#end_times = list(set([x[-22:-8] for x in glob.glob(base_dir +'20*')]))
end_times = [dt.datetime.strftime((dt.datetime.strptime(rt,'%Y%m%d%H%M%S') - dt.timedelta(hours=72)),'%Y%m%d%H%M%S') for rt in release_times]

pressure=3730

for i in range(0,len(release_times)):
    rt = release_times[i]
    rt_dt = dt.datetime.strptime(rt,'%Y%m%d%H%M%S')
    et = end_times[i]
    
    fil = base_dir + '%s-%s_REL_%s/'%(rt,et,pressure) + 'output/trajectories.txt'
    
    # Get trajectory
    
    all_dat = pd.read_fwf(fil,skiprows=6,header=None,widths=widths)
    del all_dat[0]

    
    # Extract retroplume centroid - mean trajectory data
    mean_df = all_dat[all_dat.columns[0:len(col_headers)]]
    mean_df.columns = col_headers
    
    f = open(out_dir + 'tdump%s'%(rt),'w')
    
    # number of meteorological grids used in calculation
    f.write('     1     1\n')    

    # Meteorological Model identification, Date file starting year, month, day , hour, forecast hour
    f.write('    ERA5%6s%6s%6s%6s     0\n'%(str(rt_dt.year)[-2:],str(rt_dt.month),str(rt_dt.day),str(rt_dt.hour) ))
            
    # Number of trajectories in file, direction, vertical motion calculation method
    f.write('     1 BACKWARD  OMEGA  \n')
    
    # starting year, month, day, hour, latitude, longitude, level above ground (m)
    f.write('%6s%6s%6s%6s%8s%8s%8s\n'%(str(rt_dt.year)[-2:],str(rt_dt.month),str(rt_dt.day),str(rt_dt.hour),sum_lat,sum_lon,'4000'))
            
    # number of output variables, identification of each variable
    f.write('     1 PRESSURE\n')#('     1 PRESSURE')
    # Traj number, met grid number, year, month, day, hour, minute, forecast hour, age of traj in hours, latitude, longitude, height in meters above ground, n diagnostic variables (1st is always pressure)
    f.write( '     1     1%6s%6s%6s%6s%6s     0%8s%9.3f%9.3f%9s        0\n'%(str(rt_dt.year)[-2:],str(rt_dt.month),str(rt_dt.day),str(rt_dt.hour),str(rt_dt.minute),'0.0',sum_lat,sum_lon,mean_df['meanZ'].iloc[0]))
    
    for j in range(0,len(mean_df)):
        traj_time = rt_dt + dt.timedelta(seconds=int(mean_df['time'].iloc[j]))
        f.write( '     1     1%6s%6s%6s%6s%6s     0%8s%9.3f%9.3f%9s        0\n'%(str(traj_time.year)[-2:],str(traj_time.month),str(traj_time.day),str(traj_time.hour),str(traj_time.minute),np.abs(mean_df['time'].iloc[j])/60/60,mean_df['meanLat'].iloc[j],mean_df['meanLon'].iloc[j],mean_df['meanZ'].iloc[j]))
    
    f.close()
     
