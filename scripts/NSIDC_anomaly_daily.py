#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 09:42:49 2020
Climatological variability

This script plots the std deviation of the anomaly using
NSIDC data and daily anomalies
The anomaly is computed from the daily value - climatology
@author: mvichi
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean

import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

#%%

# This is the data projection, that is their original coordinates
data_crs = ccrs.Stereographic(-90, 0)
# The projection keyword determines how the plot will look
map_proj = ccrs.SouthPolarStereo()

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
             subplot_kws={'projection': map_proj})  # the plot's projection
for ax in p.axes.flat:
    ax.coastlines()
    ax.gridlines()
    ax.set_extent([-180, 180, -90, -55], ccrs.PlateCarree())
    
#%% Monthly timeseries of stdev anomaly NASA
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



