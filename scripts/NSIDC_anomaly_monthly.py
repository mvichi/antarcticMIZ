#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 09:42:49 2020
Climatological variability

This script plots the std deviation of the anomaly using
NSIDC data and monthly anomalies
The anomaly is computed from the monthly mean - climatology
@author: mvichi
"""
import xarray as xr
import matplotlib.pyplot as plt
import cmocean
import cartopy.crs as ccrs
import calendar
#%% set the directory to where your processed data are
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'

# This is the data projection 
# NSIDC data are on the stereographic polar projection
data_crs = ccrs.Stereographic(-90, 0)
# The map projection keyword determines how the plot will look
map_proj = ccrs.SouthPolarStereo()
# get the months name
# months = calendar.month_name[1:]
months = calendar.month_abbr[1:]
#%% std anomaly from monthly data
df = DIR+'seaice_conc_monanomstd_sh_1978_2019_v03r01.nc'
anom = xr.open_dataset(df)
field = anom.seaice_conc_cdr.where(anom.seaice_conc_cdr>0)
cmap = plt.cm.get_cmap('cmo.matter', 10)
# use the faceting capability of xarray
p = field.plot(figsize=(8.9,9.5),transform=data_crs,  # the data's projection
             col='time', col_wrap=3,  # multiplot settings
             aspect=len(field.xgrid) / len(field.ygrid),  # for a sensible figsize
             cmap=cmap, vmin=0, vmax=0.5,
             subplot_kws={'projection': map_proj},
             cbar_kwargs={'label': 'Stdev of sea ice concentration anomaly from monthly data [-]'}) 
for i,ax in enumerate(p.axes.flat):
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
    ax.set_title(months[i])