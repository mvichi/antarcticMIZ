#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 09:42:49 2020
Identification of climatological MIZ variability

This script plots the std deviation of the anomaly using
NSIDC data
The anomaly is obtained from the daily value - MONTHLY climatology and the monthly
standard deviation is computed for each month and for a climatological month
@author: mvichi
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean # using this through the "cmo" variable, so must be imported
#import pandas as pd
import cartopy.crs as ccrs
import calendar

#%% ll2grid function
# this script finds the indices of the closest grid point to given coordinates
# on a geographical grid
# inputs: lon,lat the point coordinate
#         lona,lata the grid arrays  
def sphdist(lon,lat,lona,lata):
    EARTH_RADIUS = 6376.8
    x = lon*np.pi/180
    y = lat*np.pi/180
    xa = lona*np.pi/180
    ya = lata*np.pi/180
    r = 2*(1 - np.cos(y)*np.cos(ya)*np.cos(x-xa)-np.sin(y)*np.sin(ya))
    return np.sqrt(r)*EARTH_RADIUS

def ll2grid(lon,lat,lona,lata):
    dist = sphdist(lon, lat, lona, lata)
    d = dist.values
    return np.unravel_index(np.argmin(d, axis=None), d.shape)
#%% Plot the climatological stdev of the daily anomaly for each month
# This is the data projection, that is their original coordinates
data_crs = ccrs.Stereographic(-90, 0)
# The projection keyword determines how the plot will look
map_proj = ccrs.SouthPolarStereo()
# get the months name
months = calendar.month_name[1:]
#months = calendar.month_abbr[1:]

DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
df=DIR+'NSIDC_cdr_m_sigmaSIAclim.nc'
anom=xr.open_dataset(df)
data_crs = ccrs.Stereographic(-90, 0)
field = anom.sigmaSIAclim.where(anom.sigmaSIAclim>0)
cmap=plt.cm.get_cmap('cmo.matter', 10)
p = field.plot(figsize=(8.9,9.5),transform=data_crs,  # the data's projection
             col='time', col_wrap=3,  # multiplot settings
             aspect=len(field.xgrid) / len(field.ygrid),  # for a sensible figsize
             cmap=cmap, vmin=0, vmax=0.5,
             subplot_kws={'projection': map_proj},
             cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
                'label': 'Climatological $\sigma_{SIA}$ NOAA/NSIDC CDR [-]'})
for i,ax in enumerate(p.axes.flat):
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
    ax.set_title(months[i])

#%% This is the climatological mean of the monthly stdev, not pooled together
# it is much less intense (does not measure the real variance, not used in the paper) 

# DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
# df=DIR+'seaice_conc_clim_dayanom_monstd_sh_1978_2019_v03r01.nc'
# anom=xr.open_dataset(df)
# data_crs = ccrs.Stereographic(-90, 0)
# field = anom.seaice_conc_cdr.where(anom.seaice_conc_cdr>0)
# cmap=plt.cm.get_cmap('cmo.matter', 10)
# p = field.plot(figsize=(8.9,9.5),transform=data_crs,  # the data's projection
#              col='time', col_wrap=3,  # multiplot settings
#              aspect=len(field.xgrid) / len(field.ygrid),  # for a sensible figsize
#              cmap=cmap, vmin=0, vmax=0.5,
#              subplot_kws={'projection': map_proj},
#              cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
#                 'label': 'Stdev of sea ice concentration anomaly from daily data [-]'})
# for i,ax in enumerate(p.axes.flat):
#     ax.gridlines()
#     ax.coastlines()
#     ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
#     ax.set_title(months[i])

#%% Compute the median of the monthly sigmaSIA
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
df=DIR+'NSIDC_cdr_m_sigmaSIA.nc'
anom=xr.open_dataset(df)
# exclude gridpoints with no variability. This is
# preferred to identify a threshold for MIZ variability that is unbiased
# It may enhance the variance in open ocean points that are 
# occasionally ice covered
q = anom.sigmaSIA.where(anom.sigmaSIA>0.).median(dim='time')
# this computation includes all months with standard deviation = 0,
# which weights more the periods of ice-free ocean. 
# q0 = anom.seaice_conc_cdr.median(dim='time')

# compute the CDF
p_cdf = np.arange(0,1.001,0.001)
quant = q.quantile(p_cdf)
#%% Extract timeseries of monthly sigmaSIA at selected locations
# uses datasets: anom (monthly stds) and q (median)
def extract_ts(lon,lat):
    ind=ll2grid(lon, lat, q.longitude, q.latitude) # find the indices of the location
    ts = anom.sigmaSIA.isel(ygrid=ind[0],xgrid=ind[1])
    thres = q.isel(ygrid=ind[0],xgrid=ind[1]).values
    return ts.to_pandas(),thres # convert to pandas for slicing, otherwise error
def extract_ts_nonzero(lon,lat):
    ind=ll2grid(lon, lat, q.longitude, q.latitude) # find the indices of the location
    ts = anom.sigmaSIA.where(anom.sigmaSIA>0.).isel(ygrid=ind[0],xgrid=ind[1])
    thres = q.isel(ygrid=ind[0],xgrid=ind[1]).values
    return ts.to_pandas(),thres # convert to pandas for slicing, otherwise error

#ts.sel(time=slice('2010-01-16','2011-12-16')).plot() # this one does not work

#%% Figure 1
fig1 = plt.figure(figsize=(8,7),constrained_layout=True)
gs = fig1.add_gridspec(5, 2)
f1_ax2 = fig1.add_subplot(gs[3:, 0])
# locations for timeseries analyses
points=[
 [-40,[-59,-62,-70]],
 [0,[-56,-60,-68]],
 [60,[-60,-63,-67]],
 [120,[-62,-64,-66]],
 [180,[-65,-70,-76]]
 ]

q.plot.hist(ax=f1_ax2,bins=np.arange(0,0.6,0.02),
            density=True,histtype='step',color='k')
f1_ax2.set_title('')
f1_ax2.set_xlabel('Median of $\sigma_{SIA}$ NOAA/NSIDC CDR [-]')
f1_ax2.set_ylabel('PDF')
f1_ax2.set_xlim(-0.02,0.3)
ax2t = f1_ax2.twinx()
ax2t.plot(quant.values,p_cdf)
ax2t.set_ylabel('CDF')

for i in range(5):
    ax=fig1.add_subplot(gs[i, 1])
    p=points[i]
    lon=p[0]
    lats=p[1]
    data=[]
    labels=[]
    for n,l in enumerate(lats):
        tsp,thres = extract_ts_nonzero(lon,l)
        data.append(tsp.dropna())
        labels.append('{} ({:g})'.format(l,thres))
    ax.boxplot(data,labels=labels)
    ax.set_ylim(0,0.5)
    ax.set_ylabel('$\sigma_{SIA}$')
    ax.set_title('Longitude {}'.format(lon))
ax.set_xlabel('Latitude (median)')
fig1.text(0.01,0.98,'a',fontsize=14,fontweight='bold')
fig1.text(0.01,0.40,'b',fontsize=14,fontweight='bold')
fig1.text(0.97,0.98,'c',fontsize=14,fontweight='bold')
#%% Figure 1 panel a: sigmaSIA median
# This is used to define a threshold of variability that is indicative of 
# the MIZ

# This is the data projection, original coordinates of the dataset
data_crs = ccrs.Stereographic(-90, 0)
# The projection keyword determines how the plot will look
map_proj = ccrs.SouthPolarStereo()

fmed=plt.figure(figsize=(6, 6))
axmed = plt.axes(projection=map_proj)
#ax.coastlines()
x=q.xgrid
y=q.ygrid
axmed.set_extent([x.min(), x.max(), y.min()+0.7e6, y.max()],
                ccrs.Stereographic(-90, 0))
#ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
axmed.gridlines(draw_labels=True)
cmap=plt.cm.get_cmap('cmo.matter', 5)
q.plot(ax=axmed,transform=data_crs,
       cmap=cmap,vmin=0,vmax=0.25,
       cbar_kwargs={'orientation':'horizontal',
                    'shrink':0.8,'pad':0.04,
                    'label':'Median of $\sigma_{SIA}$ NOAA/NSIDC CDR [-]'})

#%% overlay selected points
# locations for timeseries analyses
points=[
 [-40,[-59,-62,-70]],
 [0,[-56,-60,-68]],
 [60,[-60,-63,-67]],
 [120,[-62,-64,-66]],
 [180,[-65,-70,-76]]
 ]
#[-80,[-66,-69,-71]] Amundsen-Bellingshausen

# execute this after the figure is produced
plt.tight_layout(pad=1.15)

for p in points:
    lon=p[0]
    lats=p[1]
    for l in lats:
        axmed.plot(lon,l,marker='o',markerfacecolor='None',
                   color='k',transform=ccrs.PlateCarree())
# axmed.text(42,-44,'a',fontsize=14,
#             fontweight='bold',transform=ccrs.PlateCarree())

#%% Figure 1 Panel b: plot the distribution of medians
f=plt.figure()
q.plot.hist()
plt.title('')
plt.xlabel('Median of $\sigma_{SIA}$ NOAA/NSIDC CDR [-]')
plt.ylabel('Counts')
plt.xlim(0,0.3)
print(q.quantile(0.99))
#%% Figure 1 panel c: Point Distributions (box-n-whisker)
f,axes = plt.subplots(5,1,figsize=(6,9))
for i,ax in enumerate(axes.flat):
    p=points[i]
    lon=p[0]
    lats=p[1]
    data=[]
    labels=[]
    for n,l in enumerate(lats):
        tsp,thres = extract_ts_nonzero(lon,l)
        data.append(tsp.dropna())
        labels.append('{} ({:g})'.format(l,thres))
    ax.boxplot(data,labels=labels)
    ax.set_ylim(0,0.5)
    ax.set_ylabel('$\sigma_{SIA}$ NOAA/NSIDC CDR [-]')
    ax.set_title('Longitude {}'.format(lon))
plt.tight_layout()

#%% Point timeseries 
# markers are plotted when stdev exists (there is sea ice)
# a missing marker indicates no sea ice during that month
f,axes = plt.subplots(5,1,figsize=(12,10))
tstart='2010-01-01'
tend='2015-01-01'
for i,ax in enumerate(axes.flat):
    p=points[i]
    lon=p[0]
    lats=p[1]
    for l in lats:
        tsp,thres = extract_ts(lon,l)
        tsp.loc[tstart:tend].plot(ax=ax,marker='o',
                                label='{} ({:g})'.format(l,thres))
    ax.set_ylim(0,0.42)
    ax.set_xlabel('')
    ax.set_ylabel('$\sigma_{SIA}$ NOAA/NSIDC CDR [-]')
    ax.legend()
    ax.grid(axis='x')
    ax.set_title('Longitude {}'.format(lon))
plt.tight_layout()

#%% Point Distributions
# timeseries of monthly stdevanom at selected locations where std>0
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

f,axes = plt.subplots(5,1,figsize=(12,10))

for i,ax in enumerate(axes.flat):
    p=points[i]
    lon=p[0]
    lats=p[1]
    for n,l in enumerate(lats):
        tsp,thres = extract_ts_nonzero(lon,l)
        tsp.hist(ax=ax,bins=np.arange(0, 0.5, 0.025), density=True,
                 label='{} ({:g})'.format(l,thres),color=colors[n])
        ax.plot([thres,thres],[0,25],color=colors[n])
    ax.set_xlim(0,0.5)
    ax.set_ylim(0,25)
    ax.set_ylabel('Density')
    ax.legend()
    ax.set_title('Longitude {}'.format(lon))
ax.set_xlabel('$\sigma_{SIA}$ NOAA/NSIDC CDR [-]')
plt.tight_layout()

#%% Goddard merged
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
df=DIR+'seaice_conc_dayanom_monstd_sh_1978_2019_v03r01.nc'
anom=xr.open_dataset(df)
q = anom.goddard_merged_seaice_conc.where(anom.goddard_merged_seaice_conc>0.).median(dim='time')
# compute the CDF
p_cdf = np.arange(0,1.001,0.001)
quant = q.quantile(p_cdf)

#%% Goddard merged

# This is the data projection, original coordinates of the dataset
data_crs = ccrs.Stereographic(-90, 0)
# The projection keyword determines how the plot will look
map_proj = ccrs.SouthPolarStereo()

fig = plt.figure(figsize=(6,8), constrained_layout=True)
gs = fig.add_gridspec(3, 1)
f_ax1 = fig.add_subplot(gs[:2,0],projection=map_proj)
f_ax2 = fig.add_subplot(gs[2, 0])
x=q.xgrid
y=q.ygrid
f_ax1.set_extent([x.min(), x.max(), y.min()+0.7e6, y.max()],
                ccrs.Stereographic(-90, 0))
f_ax1.gridlines(draw_labels=True)
cmap=plt.cm.get_cmap('cmo.matter', 5)
q.plot(ax=f_ax1,transform=data_crs,
       cmap=cmap,vmin=0,vmax=0.25,
       cbar_kwargs={'orientation':'horizontal',
                    'shrink':0.8,'pad':0.04,
                    'label':'Median of $\sigma_{SIA}$ NOAA/NSIDC Goddard merged [-]'})

q.plot.hist(ax=f_ax2,bins=np.arange(0,0.6,0.02),
            density=True,histtype='step',color='k')
f_ax2.set_title('')
f_ax2.set_xlabel('Median of $\sigma_{SIA}$ NOAA/NSIDC Goddard merged [-]')
f_ax2.set_ylabel('PDF')
f_ax2.set_xlim(-0.02,0.3)
ax2t = f_ax2.twinx()
ax2t.plot(quant.values,p_cdf)
ax2t.set_ylabel('CDF')
#%% Monthly maps of stdevanom from NSIDC for a selected year
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
df=DIR+'seaice_conc_dayanom_monstd_sh_1978_2019_v03r01.nc'
anom=xr.open_dataset(df)
data_crs = ccrs.Stereographic(-90, 0)
field = anom.seaice_conc_cdr.sel(time='2015').where(anom.seaice_conc_cdr>0)
cmap=plt.cm.get_cmap('cmo.matter', 10)
p = field.plot(figsize=(8.9,9.5),transform=data_crs,  # the data's projection
             col='time', col_wrap=3,  # multiplot settings
             aspect=len(field.xgrid) / len(field.ygrid),  # for a sensible figsize
             cmap=cmap, vmin=0, vmax=0.5,
             subplot_kws={'projection': map_proj})  # the plot's projection
for ax in p.axes.flat:
    ax.coastlines()
    ax.gridlines()
    ax.set_extent([-180, 180, -90, -55], ccrs.PlateCarree())

#%% Plot the climatological mask based on the median criterion for each month
# This is the data projection, that is their original coordinates
data_crs = ccrs.Stereographic(-90, 0)
# The projection keyword determines how the plot will look
map_proj = ccrs.SouthPolarStereo()
# get the months name
months = calendar.month_name[1:]
#months = calendar.month_abbr[1:]

DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
df=DIR+'seaice_conc_dayanomstd_sh_1978_2019_v03r01.nc'
anom=xr.open_dataset(df)
data_crs = ccrs.Stereographic(-90, 0)
field = anom.goddard_merged_seaice_conc.where(anom.goddard_merged_seaice_conc>0)
cmap=plt.cm.get_cmap('cmo.matter', 10)
p = field.plot(figsize=(8.9,9.5),transform=data_crs,  # the data's projection
             col='time', col_wrap=3,  # multiplot settings
             aspect=len(field.xgrid) / len(field.ygrid),  # for a sensible figsize
             cmap=cmap, vmin=0, vmax=0.5,
             subplot_kws={'projection': map_proj},
             cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
                'label': 'Stdev of sea ice concentration anomaly from daily data [-]'})
for i,ax in enumerate(p.axes.flat):
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
    ax.set_title(months[i])

