#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 21:25:39 2020

@author: heather
"""

import pandas as pd
import numpy as np
from matplotlib import rcParams, rc
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
hfont = {'fontname':'Helvetica'}


# Plotting preferences:
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams.update({'font.size': 10}) 
rcParams['axes.titlepad'] = 42 
rcParams['xtick.major.pad']='8'
rcParams['ytick.major.pad']='8'
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


run_length = 72 # hours
#sum_lat = 27.85
#sum_lon = 86.75 

sum_lat = 27.9952
sum_lon = 86.8406


#levels = [0.1,0.3,0.8,2,6,20,50,150,400]
levels = [0.001,0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1, 5, 10, 25, 50, 100]


# create a color map for these values, based on the number of values:
col_map = [plt.cm.jet((i / len(levels))) for i in range(len(levels))]
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])


# Get event files
#base_dir='/Users/heather/Desktop/Everest/everest-flexpart/'
base_dir='/nobackup/eehgu/everest-flexpart/Fani/'

release_times = list(set([x[-37:-23] for x in glob.glob(base_dir +'20*')]))
#end_times = list(set([x[-22:-8] for x in glob.glob(base_dir +'20*')]))
end_times = [dt.datetime.strftime((dt.datetime.strptime(rt,'%Y%m%d%H%M%S') - dt.timedelta(hours=72)),'%Y%m%d%H%M%S') for rt in release_times]
#end_times = [dt.datetime.strftime((dt.datetime.strptime(rt,'%Y%m%d%H%M%S') - dt.timedelta(hours=48)),'%Y%m%d%H%M%S') for rt in release_times]


#pressures=[650,525,400]
pressures=[500,400,300]

for i in range(0,len(release_times)):
    rt = release_times[i]
    et = end_times[i]
    
    gf_400 = base_dir + '%s-%s_REL_%s/'%(rt,et,300) + 'output/grid_time_%s.nc'%rt
    gf_525 = base_dir + '%s-%s_REL_%s/'%(rt,et,400) + 'output/grid_time_%s.nc'%rt
    gf_650 = base_dir + '%s-%s_REL_%s/'%(rt,et,500) + 'output/grid_time_%s.nc'%rt
        
    # load netcdf data using iris:
    try:
        cube400 = iris.load_cube(gf_400, 'spec001_mr')
        cube525 = iris.load_cube(gf_525, 'spec001_mr')
        cube650 = iris.load_cube(gf_650, 'spec001_mr')
    except:
        print('cant find output file, %s'%rt)
        continue
    
    # get lat and lon coordinates:
    lat_coord = cube400.coord('grid_latitude')
    lon_coord = cube400.coord('grid_longitude')
    
    # get lat and lon values:
    lat_vals = lat_coord.points
    lon_vals = lon_coord.points

    # get height values from data:
    height_coord = cube400.coord('height')
    height_vals = height_coord.points
    
    # sum the values:
    sum_400 = cube400[0, 0, :, :].collapsed('time', iris.analysis.SUM)
    sum_400 = sum_400.collapsed('height', iris.analysis.SUM)
    
    sum_525 = cube525[0, 0, :, :].collapsed('time', iris.analysis.SUM)
    sum_525 = sum_525.collapsed('height', iris.analysis.SUM)

    sum_650 = cube650[0, 0, :, :].collapsed('time', iris.analysis.SUM)
    sum_650 = sum_650.collapsed('height', iris.analysis.SUM)

    # mask any 0 values:
    sum_400_mask = np.ones(sum_400.shape)
    sum_400_mask[sum_400.data > 0] = 0
    sum_400.data = np.ma.array(sum_400.data)
    sum_400.data.mask = sum_400_mask
    sum_400 = (sum_400.data / np.max(sum_400.data)) * 100
    
    sum_525_mask = np.ones(sum_525.shape)
    sum_525_mask[sum_525.data > 0] = 0
    sum_525.data = np.ma.array(sum_525.data)
    sum_525.data.mask = sum_525_mask
    sum_525 = (sum_525.data / np.max(sum_525.data)) * 100
    
    sum_650_mask = np.ones(sum_650.shape)
    sum_650_mask[sum_650.data > 0] = 0
    sum_650.data = np.ma.array(sum_650.data)
    sum_650.data.mask = sum_650_mask
    sum_650 = (sum_650.data / np.max(sum_650.data)) * 100

    # Make figure

    fig = plt.figure(figsize=[12, 5])

    ax1 = plt.subplot(1, 3, 1, projection=ccrs.Mercator(central_longitude=86.9))
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

    contour_plot1 = ax1.contourf(lon_vals, lat_vals,sum_650,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
    ax1.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)
    #cb1 = plt.colorbar(contour_plot1, ax=ax1, shrink=0.8)
    #cb1.ax.set_ylabel('Emission sensitivity (%)')
    ax1.set_title('%s: 500 hPa release'%rt)


    ax2 = plt.subplot(1, 3, 2, projection=ccrs.Mercator(central_longitude=86.9))
    ax2.set_extent([sum_lon-40, 100, sum_lat-20, sum_lat + 15], ccrs.PlateCarree())
    ax2.add_feature(land_50m)
    ax2.coastlines(resolution='50m',zorder=20)
    ax2.add_feature(cfeature.BORDERS.with_scale('50m'))

    gl2 = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl2.xlabels_top = False
    gl2.ylabels_right = False
    gl2.xlocator = mticker.FixedLocator([50,60,70,80,90,100,110])
    gl2.ylocator = mticker.FixedLocator([0,10,20,30,40,50])
    gl2.xformatter = LONGITUDE_FORMATTER
    gl2.yformatter = LATITUDE_FORMATTER

    contour_plot2 = ax2.contourf(lon_vals, lat_vals,sum_525,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
    ax2.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)
    #cb2 = plt.colorbar(contour_plot2, ax=ax2, shrink=0.8)
    #cb2.ax.set_ylabel('Emission sensitivity (%)')
    ax2.set_title('%s: 400 hPa release'%rt)


    ax3 = plt.subplot(1, 3, 3, projection=ccrs.Mercator(central_longitude=86.9))
    ax3.set_extent([sum_lon-40, 100, sum_lat-20, sum_lat + 15], ccrs.PlateCarree())
    ax3.add_feature(land_50m)
    ax3.coastlines(resolution='50m',zorder=20)
    ax3.add_feature(cfeature.BORDERS.with_scale('50m'))

    gl3 = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl3.xlabels_top = False
    gl3.ylabels_right = False
    gl3.xlocator = mticker.FixedLocator([50,60,70,80,90,100,110])
    gl3.ylocator = mticker.FixedLocator([0,10,20,30,40,50])
    gl3.xformatter = LONGITUDE_FORMATTER
    gl3.yformatter = LATITUDE_FORMATTER

    contour_plot3 = ax3.contourf(lon_vals, lat_vals,sum_400,transform=ccrs.PlateCarree(),zorder=10,alpha=0.5,levels=levels,extend='max',colors=col_map)#,vmin=0, vmax=100, zorder=10, alpha=0.9,extend='min')
    ax3.plot(sum_lon, sum_lat, 'kx',markersize=10, transform=ccrs.PlateCarree(),zorder=30)
    #cb3 = plt.colorbar(contour_plot3, ax=ax3, shrink=0.8)
    #cb3.ax.set_ylabel('Emission sensitivity (%)')
    ax3.set_title('%s: 300 hPa release'%rt)

  
    cbar_ax = fig.add_axes([0.13, 0.09, 0.75, 0.03])#0.02-0.05

    cb1=fig.colorbar(contour_plot1, cax=cbar_ax,orientation='horizontal')
    cb1.ax.set_xlabel('EMISSION SENSITIVY (%)',**hfont)

    fig.subplots_adjust(bottom=0.2)
    fig.tight_layout()



    # Save plot
    save_loc = base_dir + 'out_figures/%s_gridded.tiff'%rt 
    print('Saving %s'%rt)
    fig.savefig(save_loc,dpi=300,format='tiff')

    plt.close(fig)
    fig.clf()
