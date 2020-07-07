#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 21:25:39 2020

@author: heather
"""

import pandas as pd
import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as plt
import datetime as dt
import glob as glob
import os
import iris
import warnings
warnings.filterwarnings("ignore")


import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


# Plotting preferences:
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams.update({'font.size': 22}) 
rcParams['axes.titlepad'] = 22 
rcParams['xtick.major.pad']='8'
rcParams['ytick.major.pad']='8'


run_length = 72 # hours
#sum_lat = 27.85
#sum_lon = 86.75 
sum_lat = 27.9952
sum_lon = 86.8406


col_headers = ['time','meanLon','meanLat','meanZ','meanTopo','meanPBL','meanTropo','meanPv','rmsHBefore','rmsHAfter','rmsVBefore','rmsVAfter','pblFract','pv2Fract','tropoFract']
nclusters=5
mean_widths = [5,8,9,9,8,8,8,8,8,8,8,8,8,6,6,6]
cluster_widths = [8,8,7,6,8]
widths = mean_widths + cluster_widths * nclusters
levs = [0.1,1,2,5,9,14,20,100]
col_map = [plt.cm.jet((i / len(levs))) for i in range(len(levs))]


# Get event files
#base_dir='/Users/heather/Desktop/Everest/everest-flexpart/'
base_dir='/nobackup/eehgu/everest-flexpart/April17/'

release_times = list(set([x[-37:-23] for x in glob.glob(base_dir +'20*')]))
#end_times = list(set([x[-22:-8] for x in glob.glob(base_dir +'20*')]))
end_times = [dt.datetime.strftime((dt.datetime.strptime(rt,'%Y%m%d%H%M%S') - dt.timedelta(hours=72)),'%Y%m%d%H%M%S') for rt in release_times]

pressures=[500, 400, 300]

for i in range(0,len(release_times)):
    rt = release_times[i]
    et = end_times[i]

    for pressure in pressures:   
        fil = base_dir + '%s-%s_REL_%s/'%(rt,et,pressure) + 'output/trajectories.txt'
    
        # Get trajectory
        all_dat = pd.read_fwf(fil,skiprows=6,header=None,widths=widths)
        del all_dat[0]

        # Extract retroplume centroid - mean trajectory data
        mean_df = all_dat[all_dat.columns[0:len(col_headers)]]
        mean_df.columns = col_headers
        
        # Calculate heigh above topography (magl)
        meanZagl = mean_df['meanZ'] - mean_df['meanTopo']

        # Calculate tropopause height above topography (magl)
        meanTropoagl = mean_df['meanTropo'] - mean_df['meanTopo']
    
        # Get gridded
        gfil = base_dir + '%s-%s_REL_%s/'%(rt,et,pressure) + 'output/grid_time_%s.nc'%rt
        cube = iris.load_cube(gfil, 'spec001_mr')
        #except:
        #print('cant find output file, %s'%rt)
        #continue
    
        # get height values from data:
        height_coord = cube.coord('height')
        height_vals = height_coord.points
           
        # sum the values:
        cube_sum = cube.collapsed('grid_latitude', iris.analysis.SUM)
        cube_sum = cube_sum.collapsed('grid_longitude', iris.analysis.SUM)
        
        cube_data = cube_sum[0,0,:,:].data
        
        # convert to percentages before plotting
        cube_data = (cube_data / np.max(cube_data))*100

        # Get gridded time
        gridded_time = cube_sum.coord('time').points/60/60
    


        # Make figure

        fig = plt.figure(figsize=(14,5))
        ax1 = fig.add_subplot(111)
        ax1.grid(True)
        
        cp = ax1.contourf(-gridded_time,height_vals,np.transpose(cube_data),zorder=1,alpha=0.5,levels=levs,colors=col_map)#,vmin=0.01,vmax=100,cmap='jet')#,colors=col_map,extend='max')#norm=matplotlib.colors.LogNorm(),levels=log_levs,zorder=10)

        ax1.plot(-mean_df['time'].astype(int)/60/60,meanZagl,label='Mean altitude',lw=4,zorder=20,c='k')
        ax1.plot(-mean_df['time'].astype(int)/60/60,mean_df['meanPBL'],label='PBL',lw=1,ls='--',c='k',zorder=20)
        #ax1.plot(-mean_df['time'].astype(int)/60/60,meanTropoagl,label='Tropopause',lw=1,ls=':',c='k',zorder=20)

        ax1.set_ylim(0,10000)
        ax1.set_xlim(1,72)
        ax1.set_ylabel('m agl') 
        ax1.set_xlabel('Hours before release') 
        ax1.legend(loc='upper right',facecolor='white',fontsize=16)
        
        ax1.set_xticks([12,24,36,48,60,72])  # Set label locations.

        cbar = plt.colorbar(cp,ax=ax1,shrink=0.8)#,ticks=log_levs)
        cbar.set_label('%')
 
        fig.tight_layout()
        
        # Save figure

        save_loc = base_dir + 'out_figures/gridded_altitudes_%s_%s.png'%(rt,pressure) 
        print('Saving %s, %s'%(rt,pressure))
        fig.savefig(save_loc)

        plt.close(fig)
        fig.clf()



