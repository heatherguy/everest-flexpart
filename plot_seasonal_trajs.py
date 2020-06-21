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
sum_lat = 27.85
sum_lon = 86.75 
col_headers = ['time','meanLon','meanLat','meanZ','meanTopo','meanPBL','meanTropo','meanPv','rmsHBefore','rmsHAfter','rmsVBefore','rmsVAfter','pblFract','pv2Fract','tropoFract']
nclusters=5
mean_widths = [5,8,9,9,8,8,8,8,8,8,8,8,8,6,6,6]
cluster_widths = [8,8,7,6,8]
widths = mean_widths + cluster_widths * nclusters



land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])


# Get event files
#base_dir='/Users/heather/Desktop/Everest/everest-flexpart/'
base_dir='/nobackup/eehgu/everest-flexpart/events/'

release_times = list(set([x[-37:-23] for x in glob.glob(base_dir +'20*')]))
#end_times = list(set([x[-22:-8] for x in glob.glob(base_dir +'20*')]))
end_times = [dt.datetime.strftime((dt.datetime.strptime(rt,'%Y%m%d%H%M%S') - dt.timedelta(hours=72)),'%Y%m%d%H%M%S') for rt in release_times]

pressure=650
all_cubes=[]
JJAS_cubes=[]
SO_cubes=[]
DJF_cubes=[]
M_cubes=[]

for i in range(0,len(release_times)):
    rt = release_times[i]
    et = end_times[i]
    
    fil = base_dir + '%s-%s_REL_%s/'%(rt,et,pressure) + 'output/trajectories.txt'
    
    # Get trajectory
    
    all_dat = pd.read_fwf(fil,skiprows=6,header=None,widths=widths)
    del all_dat[0]

    
    # Extract retroplume centroid - mean trajectory data
    mean_df = all_dat[all_dat.columns[0:len(col_headers)]]
    mean_df.columns = col_headers
    
    all_cubes.append(mean_df)
    
    if rt[4:6] in ['06','07','08','09']:
        JJAS_cubes.append(mean_df)
    elif rt[4:6] in ['09','10']:
        SO_cubes.append(mean_df)
    elif rt[4:6] in ['12','01','02']:
        DJF_cubes.append(mean_df)
    elif rt[4:6] in ['03']:
        M_cubes.append(mean_df)




# Make figure - all plumes

fig = plt.figure(figsize=[10, 10])

ax1 = plt.subplot(1, 1, 1, projection=ccrs.Mercator(central_longitude=86.9))
ax1.set_extent([sum_lon-40, 100, sum_lat-20, sum_lat + 15], ccrs.PlateCarree())
ax1.add_feature(land_50m)
ax1.coastlines(resolution='50m',zorder=20)
ax1.add_feature(cfeature.BORDERS.with_scale('50m'))

gl1 = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl1.xlabels_top = False
gl1.ylabels_right = False
gl1.xlocator = mticker.FixedLocator([50,60,70,80,90,100,110])
gl1.ylocator = mticker.FixedLocator([0,10,20,30,40,50])
gl1.xformatter = LONGITUDE_FORMATTER
gl1.yformatter = LATITUDE_FORMATTER

#contour_plot1 = ax1.contourf(lon_vals, lat_vals,sumall,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
ax1.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)

# Plot retroplume centroid (mean trajectory)
for i in range(0,len(all_cubes)):
    mean_df=all_cubes[i]
    ax1.plot(mean_df['meanLon'], mean_df['meanLat'], transform=ccrs.PlateCarree(),zorder=20,lw=2,c='b')

#cb1 = plt.colorbar(contour_plot1, ax=ax1, shrink=0.8)
#cb1.ax.set_ylabel('Emission sensitivity (%)')
ax1.set_title('All events, n=%s'%len(all_cubes))

fig.tight_layout()


# Save plot
save_loc = base_dir + 'out_figures/AllEvents_traj.png' 
print('Saving all plumes')
fig.savefig(save_loc)

plt.close(fig)
fig.clf()



fig = plt.figure(figsize=[10, 10])

ax1 = plt.subplot(1, 1, 1, projection=ccrs.Mercator(central_longitude=86.9))
ax1.set_extent([sum_lon-40, 100, sum_lat-20, sum_lat + 15], ccrs.PlateCarree())
ax1.add_feature(land_50m)
ax1.coastlines(resolution='50m',zorder=20)
ax1.add_feature(cfeature.BORDERS.with_scale('50m'))

gl1 = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl1.xlabels_top = False
gl1.ylabels_right = False
gl1.xlocator = mticker.FixedLocator([50,60,70,80,90,100,110])
gl1.ylocator = mticker.FixedLocator([0,10,20,30,40,50])
gl1.xformatter = LONGITUDE_FORMATTER
gl1.yformatter = LATITUDE_FORMATTER

#contour_plot1 = ax1.contourf(lon_vals, lat_vals,sumall,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
ax1.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)

# Plot retroplume centroid (mean trajectory)
for i in range(0,len(JJAS_cubes)):
    mean_df=JJAS_cubes[i]
    ax1.plot(mean_df['meanLon'], mean_df['meanLat'], transform=ccrs.PlateCarree(),zorder=20,lw=2,c='b')

#cb1 = plt.colorbar(contour_plot1, ax=ax1, shrink=0.8)
#cb1.ax.set_ylabel('Emission sensitivity (%)')
ax1.set_title('JJAS, n=%s'%len(JJAS_cubes))

fig.tight_layout()


# Save plot
save_loc = base_dir + 'out_figures/JJAS_traj.png' 
print('Saving JJAS plumes')
fig.savefig(save_loc)

plt.close(fig)
fig.clf()




fig = plt.figure(figsize=[10, 10])

ax1 = plt.subplot(1, 1, 1, projection=ccrs.Mercator(central_longitude=86.9))
ax1.set_extent([sum_lon-40, 100, sum_lat-20, sum_lat + 15], ccrs.PlateCarree())
ax1.add_feature(land_50m)
ax1.coastlines(resolution='50m',zorder=20)
ax1.add_feature(cfeature.BORDERS.with_scale('50m'))

gl1 = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl1.xlabels_top = False
gl1.ylabels_right = False
gl1.xlocator = mticker.FixedLocator([50,60,70,80,90,100,110])
gl1.ylocator = mticker.FixedLocator([0,10,20,30,40,50])
gl1.xformatter = LONGITUDE_FORMATTER
gl1.yformatter = LATITUDE_FORMATTER

#contour_plot1 = ax1.contourf(lon_vals, lat_vals,sumall,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
ax1.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)

# Plot retroplume centroid (mean trajectory)
for i in range(0,len(DJF_cubes)):
    mean_df=DJF_cubes[i]
    ax1.plot(mean_df['meanLon'], mean_df['meanLat'], transform=ccrs.PlateCarree(),zorder=20,lw=2,c='b')#alpha=0.8)

#cb1 = plt.colorbar(contour_plot1, ax=ax1, shrink=0.8)
#cb1.ax.set_ylabel('Emission sensitivity (%)')
ax1.set_title('DJF, n=%s'%len(DJF_cubes))

fig.tight_layout()


# Save plot
save_loc = base_dir + 'out_figures/DJF_traj.png' 
print('Saving DJF plumes')
fig.savefig(save_loc)

plt.close(fig)
fig.clf()




fig = plt.figure(figsize=[10, 10])

ax1 = plt.subplot(1, 1, 1, projection=ccrs.Mercator(central_longitude=86.9))
ax1.set_extent([sum_lon-40, 100, sum_lat-20, sum_lat + 15], ccrs.PlateCarree())
ax1.add_feature(land_50m)
ax1.coastlines(resolution='50m',zorder=20)
ax1.add_feature(cfeature.BORDERS.with_scale('50m'))

gl1 = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl1.xlabels_top = False
gl1.ylabels_right = False
gl1.xlocator = mticker.FixedLocator([50,60,70,80,90,100,110])
gl1.ylocator = mticker.FixedLocator([0,10,20,30,40,50])
gl1.xformatter = LONGITUDE_FORMATTER
gl1.yformatter = LATITUDE_FORMATTER

#contour_plot1 = ax1.contourf(lon_vals, lat_vals,sumall,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
ax1.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)

# Plot retroplume centroid (mean trajectory)
for i in range(0,len(SO_cubes)):
    mean_df=SO_cubes[i]
    ax1.plot(mean_df['meanLon'], mean_df['meanLat'], transform=ccrs.PlateCarree(),zorder=20,lw=2,c='b')

#cb1 = plt.colorbar(contour_plot1, ax=ax1, shrink=0.8)
#cb1.ax.set_ylabel('Emission sensitivity (%)')
ax1.set_title('SO, n=%s'%len(SO_cubes))

fig.tight_layout()


# Save plot
save_loc = base_dir + 'out_figures/SO_traj.png' 
print('Saving SO plumes')
fig.savefig(save_loc)

plt.close(fig)
fig.clf()




fig = plt.figure(figsize=[10, 10])

ax1 = plt.subplot(1, 1, 1, projection=ccrs.Mercator(central_longitude=86.9))
ax1.set_extent([sum_lon-40, 100, sum_lat-20, sum_lat + 15], ccrs.PlateCarree())
ax1.add_feature(land_50m)
ax1.coastlines(resolution='50m',zorder=20)
ax1.add_feature(cfeature.BORDERS.with_scale('50m'))

gl1 = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl1.xlabels_top = False
gl1.ylabels_right = False
gl1.xlocator = mticker.FixedLocator([50,60,70,80,90,100,110])
gl1.ylocator = mticker.FixedLocator([0,10,20,30,40,50])
gl1.xformatter = LONGITUDE_FORMATTER
gl1.yformatter = LATITUDE_FORMATTER

#contour_plot1 = ax1.contourf(lon_vals, lat_vals,sumall,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
ax1.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)

# Plot retroplume centroid (mean trajectory)
for i in range(0,len(M_cubes)):
    mean_df=M_cubes[i]
    ax1.plot(mean_df['meanLon'], mean_df['meanLat'], transform=ccrs.PlateCarree(),zorder=20,lw=2,c='b')

#cb1 = plt.colorbar(contour_plot1, ax=ax1, shrink=0.8)
#cb1.ax.set_ylabel('Emission sensitivity (%)')
ax1.set_title('M, n=%s'%len(M_cubes))

fig.tight_layout()


# Save plot
save_loc = base_dir + 'out_figures/M_traj.png' 
print('Saving M plumes')
fig.savefig(save_loc)

plt.close(fig)
fig.clf()