#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 09:42:49 2020
Identification of climatological MIZ variability

This script plots the std deviation of the anomaly using
NSIDC data
The anomaly is obtained from the daily value - climatology and the monthly
standard deviation is computed for each month and for a climatological month
@author: mvichi
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean
import pandas as pd
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
df=DIR+'seaice_conc_dayanomstd_sh_1978_2019_v03r01.nc'
anom=xr.open_dataset(df)
data_crs = ccrs.Stereographic(-90, 0)
field = anom.seaice_conc_cdr.where(anom.seaice_conc_cdr>0)
cmap=plt.cm.get_cmap('cmo.matter', 10)
p = field.plot(figsize=(8.9,9.5),transform=data_crs,  # the data's projection
             col='time', col_wrap=3,  # multiplot settings
             aspect=len(field.xgrid) / len(field.ygrid),  # for a sensible figsize
             cmap=cmap, vmin=0, vmax=0.5,
             subplot_kws={'projection': map_proj},
             cbar_kwargs={'label': 'Stdev of sea ice concentration anomaly from daily data [-]'})
for i,ax in enumerate(p.axes.flat):
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
    ax.set_title(months[i])
    
#%% Compute the median of the monthly std anomalies
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
df=DIR+'seaice_conc_dayanom_monstd_sh_1978_2019_v03r01.nc'
anom=xr.open_dataset(df)
# exclude gridpoints with no variability. This is
# preferred to identify a threshold for MIZ variability that is unbiased
# It may enhance the variance in open ocean points that are 
# occasionally ice covered
q = anom.seaice_conc_cdr.where(anom.seaice_conc_cdr>0.).median(dim='time')
# this computation includes all months with standard deviation = 0,
# which weights more the periods of ice-free ocean. 
# q0 = anom.seaice_conc_cdr.median(dim='time')
#%% Plot the median
# This is used to define a threshold of variability that is indicative of 
# the MIZ
plt.figure(figsize=(6, 6))
ax = plt.axes(projection=map_proj)
#ax.coastlines()
x=q.xgrid
y=q.ygrid
ax.set_extent([x.min(), x.max(), y.min()+0.7e6, y.max()],
                ccrs.Stereographic(-90, 0))
#ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
ax.gridlines(draw_labels=True)
cmap=plt.cm.get_cmap('cmo.matter', 6)
q.plot(ax=ax,transform=data_crs,
       cmap=cmap,vmin=0,vmax=0.3,
       cbar_kwargs={'orientation':'horizontal',
                    'shrink':0.8,'pad':0.05,
                    'label':'Median of the variability indicator [-]'})
#im=ax.pcolormesh(x, y, q, transform=data_crs, cmap=cmap,vmin=0,vmax=0.3)
#plt.colorbar(im,orientation='horizontal',shrink=0.8,pad=0.05)

#%% timeseries of monthly stdevanom at selected locations
ind=ll2grid(0,-60, q.longitude, q.latitude) # find the indices of the location
ts = anom.seaice_conc_cdr.isel(ygrid=ind[0],xgrid=ind[1])
tsp = ts.to_pandas() # convert to pandas for slicing, otherwise error
#ts.sel(time=slice('2010-01-16','2011-12-16')).plot() # this one does not work
plt.figure()
tsp.loc['2015-01-01':'2020-01-01'].plot(marker='o')
plt.grid(axis='x')

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



