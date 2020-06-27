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

#levels = [0.1,0.3,0.8,2,6,20,50,150,400]
levels = [0.001,0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1, 5, 10, 25, 50, 100]
# create a color map for these values, based on the number of values:
col_map = [plt.cm.jet((i / len(levels))) for i in range(len(levels))]
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])


# Get event files
#base_dir='/Users/heather/Desktop/Everest/everest-flexpart/'
base_dir='/nobackup/eehgu/everest-flexpart/events/'

release_times = list(set([x[-37:-23] for x in glob.glob(base_dir +'20*')]))
#end_times = list(set([x[-22:-8] for x in glob.glob(base_dir +'20*')]))
end_times = [dt.datetime.strftime((dt.datetime.strptime(rt,'%Y%m%d%H%M%S') - dt.timedelta(hours=72)),'%Y%m%d%H%M%S') for rt in release_times]

pressure=525
all_cubes=[]
JJAS_cubes=[]
ON_cubes=[]
DJF_cubes=[]
M_cubes=[]

for i in range(0,len(release_times)):
    rt = release_times[i]
    et = end_times[i]
    
    fil = base_dir + '%s-%s_REL_%s/'%(rt,et,pressure) + 'output/grid_time_%s.nc'%rt
    
    # load netcdf data using iris:
    try:
        cube = iris.load_cube(fil, 'spec001_mr')
    except:
        print('cant find output file, %s'%rt)
        continue
    
    # get lat and lon coordinates:
    lat_coord = cube.coord('grid_latitude')
    lon_coord = cube.coord('grid_longitude')
    
    # get lat and lon values:
    lat_vals = lat_coord.points
    lon_vals = lon_coord.points

    # get height values from data:
    height_coord = cube.coord('height')
    height_vals = height_coord.points
    
    # sum the values:
    try:
        sumcube = cube[0, 0, :, :].collapsed('time', iris.analysis.SUM)
        sumcube = sumcube.collapsed('height', iris.analysis.SUM)

        # mask any 0 values:
        sumcube_mask = np.ones(sumcube.shape)
        sumcube_mask[sumcube.data > 0] = 0
        sumcube.data = np.ma.array(sumcube.data)
        sumcube.data.mask = sumcube_mask
    
        temp=sumcube.data
        all_cubes.append(temp)
    
        if rt[4:6] in ['06','07','08','09']:
            JJAS_cubes.append(temp)
        elif rt[4:6] in ['10','11']:
            ON_cubes.append(temp)
        elif rt[4:6] in ['12','01','02']:
            DJF_cubes.append(temp)
        elif rt[4:6] in ['03']:
            M_cubes.append(temp)
    except:
        print('Cant parse %s, empty cube?'%rt)
        continue    



# Make figure - all plumes

# Sum all plumes
sumall = np.sum(all_cubes,axis=0)
#convert to percentage of max
sumall = (sumall / np.max(sumall)) * 100

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

contour_plot1 = ax1.contourf(lon_vals, lat_vals,sumall,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
ax1.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)
cb1 = plt.colorbar(contour_plot1, ax=ax1, shrink=0.8)
cb1.ax.set_ylabel('Emission sensitivity (%)')
ax1.set_title('All events, n=%s,end level= %s hPa'%(len(all_cubes),pressure))

fig.tight_layout()


# Save plot
save_loc = base_dir + 'out_figures/AllEvents_gridded_%s.png'%pressure 
print('Saving all plumes')
fig.savefig(save_loc)

plt.close(fig)
fig.clf()






# Make figure - JJAS

# Sum all plumes
sumjjas = np.sum(JJAS_cubes,axis=0)
#convert to percentage of max
sumjjas = (sumjjas / np.max(sumjjas)) * 100

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

contour_plot1 = ax1.contourf(lon_vals, lat_vals,sumjjas,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
ax1.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)
cb1 = plt.colorbar(contour_plot1, ax=ax1, shrink=0.8)
cb1.ax.set_ylabel('Emission sensitivity (%)')
ax1.set_title('JJAS events, n=%s, end level= %s hPa'%(len(JJAS_cubes),pressure))

fig.tight_layout()


# Save plot
save_loc = base_dir + 'out_figures/JJAS_gridded_%s.png'%pressure 
print('Saving JJAS')
fig.savefig(save_loc)

plt.close(fig)
fig.clf()





# Make figure - SO

# Sum all plumes
sumso = np.sum(ON_cubes,axis=0)
#convert to percentage of max
sumso = (sumso / np.max(sumso)) * 100

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

contour_plot1 = ax1.contourf(lon_vals, lat_vals,sumso,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
ax1.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)
cb1 = plt.colorbar(contour_plot1, ax=ax1, shrink=0.8)
cb1.ax.set_ylabel('Emission sensitivity (%)')
ax1.set_title('ON events, n=%s, end level= %s hPa'%(len(ON_cubes),pressure))

fig.tight_layout()


# Save plot
save_loc = base_dir + 'out_figures/ON_gridded_%s.png'%pressure  
print('Saving ON')
fig.savefig(save_loc)

plt.close(fig)
fig.clf()





# Make figure - DJF

# Sum all plumes
sumdjf = np.sum(DJF_cubes,axis=0)
#convert to percentage of max
sumdjf = (sumdjf / np.max(sumdjf)) * 100

fig = plt.figure(figsize=[15, 10])

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

contour_plot1 = ax1.contourf(lon_vals, lat_vals,sumdjf,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
ax1.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)
cb1 = plt.colorbar(contour_plot1, ax=ax1, shrink=0.8)
cb1.ax.set_ylabel('Emission sensitivity (%)')
ax1.set_title('DJF events, n=%s, end level= %s hPa'%(len(DJF_cubes),pressure))

fig.tight_layout()


# Save plot
save_loc = base_dir + 'out_figures/DJF_gridded_%s.png'%pressure 
print('Saving DJF')
fig.savefig(save_loc)

plt.close(fig)
fig.clf()


# Make figure - M

# Sum all plumes
summ = np.sum(M_cubes,axis=0)
#convert to percentage of max
summ = (summ / np.max(summ)) * 100

fig = plt.figure(figsize=[15, 10])

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

contour_plot1 = ax1.contourf(lon_vals, lat_vals,summ,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
ax1.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)
cb1 = plt.colorbar(contour_plot1, ax=ax1, shrink=0.8)
cb1.ax.set_ylabel('Emission sensitivity (%)')
ax1.set_title('M events, n=%s, end level= %s hPa'%(len(M_cubes),pressure))

fig.tight_layout()


# Save plot
save_loc = base_dir + 'out_figures/M_gridded_%s.png'%pressure 
print('Saving M')
fig.savefig(save_loc)

plt.close(fig)
fig.clf()
