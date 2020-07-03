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
rcParams.update({'font.size': 30}) 
rcParams['axes.titlepad'] = 42 
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

pressure=650
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
sumjjas = np.sum(JJAS_cubes,axis=0)
sumjjas = (sumjjas / np.max(sumjjas)) * 100

sumso = np.sum(ON_cubes,axis=0)
sumso = (sumso / np.max(sumso)) * 100

sumdjf = np.sum(DJF_cubes,axis=0)
sumdjf = (sumdjf / np.max(sumdjf)) * 100

summ = np.sum(M_cubes,axis=0)
summ = (summ / np.max(summ)) * 100





fig = plt.figure(figsize=[35, 30])

ax1 = plt.subplot(2, 2, 1, projection=ccrs.Mercator(central_longitude=86.9))
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
ax1.text(88,41,'JJAS, n=%s'%len(JJAS_cubes),transform=ccrs.PlateCarree(),bbox=dict(facecolor='white', alpha=0.9))



ax2 = plt.subplot(2, 2, 2, projection=ccrs.Mercator(central_longitude=86.9))
ax2.set_extent([sum_lon-40, 100, sum_lat-20, sum_lat + 15], ccrs.PlateCarree())
ax2.add_feature(land_50m)
ax2.coastlines(resolution='50m',zorder=20)
ax2.add_feature(cfeature.BORDERS.with_scale('50m'))

gl1 = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl1.xlabels_top = False
gl1.ylabels_right = False
gl1.xlocator = mticker.FixedLocator([50,60,70,80,90,100,110])
gl1.ylocator = mticker.FixedLocator([0,10,20,30,40,50])
gl1.xformatter = LONGITUDE_FORMATTER
gl1.yformatter = LATITUDE_FORMATTER

#contour_plot1 = ax2.contourf(lon_vals, lat_vals,sumon,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
ax2.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)
ax2.text(88,41,'ON, n=%s'%len(ON_cubes),transform=ccrs.PlateCarree(),bbox=dict(facecolor='white', alpha=0.9))


ax3 = plt.subplot(2, 2, 3, projection=ccrs.Mercator(central_longitude=86.9))
ax3.set_extent([sum_lon-40, 100, sum_lat-20, sum_lat + 15], ccrs.PlateCarree())
ax3.add_feature(land_50m)
ax3.coastlines(resolution='50m',zorder=20)
ax3.add_feature(cfeature.BORDERS.with_scale('50m'))

gl1 = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl1.xlabels_top = False
gl1.ylabels_right = False
gl1.xlocator = mticker.FixedLocator([50,60,70,80,90,100,110])
gl1.ylocator = mticker.FixedLocator([0,10,20,30,40,50])
gl1.xformatter = LONGITUDE_FORMATTER
gl1.yformatter = LATITUDE_FORMATTER

#contour_plot1 = ax3.contourf(lon_vals, lat_vals,sumdjf,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
ax3.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)
ax3.text(88,41,'DJF, n=%s'%len(DJF_cubes),transform=ccrs.PlateCarree(),bbox=dict(facecolor='white', alpha=0.9))


ax4 = plt.subplot(2, 2, 4, projection=ccrs.Mercator(central_longitude=86.9))
ax4.set_extent([sum_lon-40, 100, sum_lat-20, sum_lat + 15], ccrs.PlateCarree())
ax4.add_feature(land_50m)
ax4.coastlines(resolution='50m',zorder=20)
ax4.add_feature(cfeature.BORDERS.with_scale('50m'))

gl1 = ax4.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl1.xlabels_top = False
gl1.ylabels_right = False
gl1.xlocator = mticker.FixedLocator([50,60,70,80,90,100,110])
gl1.ylocator = mticker.FixedLocator([0,10,20,30,40,50])
gl1.xformatter = LONGITUDE_FORMATTER
gl1.yformatter = LATITUDE_FORMATTER

#contour_plot1 = ax4.contourf(lon_vals, lat_vals,summ,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
ax4.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)
ax4.text(88,41,'M, n=%s'%len(M_cubes),transform=ccrs.PlateCarree(),bbox=dict(facecolor='white', alpha=0.9))



fig.subplots_adjust(top=0.8)
cbar_ax = fig.add_axes([0.15, 0.05, 0.75, 0.02])
cb1=fig.colorbar(contour_plot1, cax=cbar_ax,orientation='horizontal')
cb1.ax.set_xlabel('Emission sensitivity (%)')

fig.tight_layout()


# Save plot
save_loc = base_dir + 'out_figures/seasonal_fourpanel.png' 
print('Saving Jfour panel')
fig.savefig(save_loc)

plt.close(fig)
fig.clf()




