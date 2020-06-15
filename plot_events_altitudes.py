#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 21:54:48 2020

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

# Plotting preferences:
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams.update({'font.size': 22}) 
rcParams['axes.titlepad'] = 22 
rcParams['xtick.major.pad']='8'
rcParams['ytick.major.pad']='8'

col_headers = ['time','meanLon','meanLat','meanZ','meanTopo','meanPBL','meanTropo','meanPv','rmsHBefore','rmsHAfter','rmsVBefore','rmsVAfter','pblFract','pv2Fract','tropoFract']
nclusters=5
mean_widths = [5,8,9,9,8,8,8,8,8,8,8,8,8,6,6,6]
cluster_widths = [8,8,7,6,8]
widths = mean_widths + cluster_widths * nclusters
sc = plt.rcParams['axes.prop_cycle'].by_key()['color']

# Get event files
#base_dir='/Users/heather/Desktop/Everest/everest-flexpart/'
base_dir='/nobackup/eehgu/everest-flexpart/events/'


release_times = list(set([x[-37:-23] for x in glob.glob(base_dir +'20*')]))
#end_times = list(set([x[63:77] for x in glob.glob(base_dir +'20*')]))
end_times = [dt.datetime.strftime((dt.datetime.strptime(rt,'%Y%m%d%H%M%S') - dt.timedelta(hours=72)),'%Y%m%d%H%M%S') for rt in release_times]


pressures=[650,525,400]

for i in range(0,len(release_times)):
    rt = release_times[i]
    et = end_times[i]

    tf_400 = base_dir + '%s-%s_REL_%s/'%(rt,et,400) + 'output/trajectories.txt'
    tf_525 = base_dir + '%s-%s_REL_%s/'%(rt,et,525) + 'output/trajectories.txt'
    tf_650 = base_dir + '%s-%s_REL_%s/'%(rt,et,650) + 'output/trajectories.txt'
    
    # Extract mean trajectories
    dat400 = pd.read_fwf(tf_400,skiprows=6,header=None,widths=widths)
    del dat400[0]
    df400 = dat400[dat400.columns[0:len(col_headers)]]
    df400.columns = col_headers
    
    dat525 = pd.read_fwf(tf_525,skiprows=6,header=None,widths=widths)
    del dat525[0]
    df525 = dat525[dat525.columns[0:len(col_headers)]]
    df525.columns = col_headers
    
    dat650 = pd.read_fwf(tf_650,skiprows=6,header=None,widths=widths)
    del dat650[0]
    df650 = dat650[dat650.columns[0:len(col_headers)]]
    df650.columns = col_headers
    
    
    
    
    
    # Plot
    fig = plt.figure(figsize=(40,10))
    ax1 = fig.add_subplot(111)

    #ax1.plot(-mean_df['time'].astype(int)/60/60/24,mean_df['meanTropo'],label='Tropopause',lw=1,ls=':',c='k')
    #ax1.fill_between(-mean_df['time'].astype(int)/60/60/24,0,df400['meanTopo'],label='Topography',fc='darkgreen')

    ax1.plot(-df400['time'].astype(int)/60/60/24,df400['meanZ'],label='400 hPa release',lw=6,c=sc[0])
    ax1.plot(-df525['time'].astype(int)/60/60/24,df525['meanZ'],label='525 hPa release',lw=6,c=sc[1])
    ax1.plot(-df650['time'].astype(int)/60/60/24,df650['meanZ'],label='650 hPa release',lw=6,c=sc[2])
    
    ax1.plot(-df400['time'].astype(int)/60/60/24,df400['meanTopo'],lw=2,ls='-.',c=sc[0])
    ax1.plot(-df525['time'].astype(int)/60/60/24,df525['meanTopo'],lw=2,ls='-.',c=sc[1])
    ax1.plot(-df650['time'].astype(int)/60/60/24,df650['meanTopo'],lw=2,ls='-.',c=sc[2])


    ax1.grid(True)
    ax1.set_ylabel('m asl')
    ax1.set_ylim(bottom=0)
    ax1.set_xlim(0,3)
    ax1.legend(loc='upper right',fontsize=26)
    ax1.set_xlabel('Days before release')
    fig.tight_layout()

    # Save plot
    save_loc = base_dir + 'out_figures/%s_altitude.png'%rt 
    print('Saving %s altitude plot'%rt)
    fig.savefig(save_loc)

    plt.close(fig)
    fig.clf()