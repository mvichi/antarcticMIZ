#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 22 07:07:37 2021
Regional Area extent and intensity
@author: vichi
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean # using this through the "cmo" variable
import pandas as pd
import calendar
# get the months name
months = calendar.month_name[1:]
#months = calendar.month_abbr[1:]
import matplotlib as mpl
mpl.rcParams['font.size'] = 12
#%% Constants
AREA = 25.*25./1.e6 # from pixel to 10^6 km2
TH = 0.1 # sigmaSIA threshold
THP = TH*100.
#%% Load the sigmaSIA for NOAA/NSIDC CDR
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
dfile=DIR+'NSIDC_cdr_m_sigmaSIA.nc'
ds=xr.open_dataset(dfile)
miz=ds.sigmaSIA.sel(time=slice('1989','2019')).where(ds.sigmaSIA>=TH)
mizint=miz.mean(dim=['xgrid','ygrid']) # measure of the intensity of variability
miz=miz/miz
mizext=miz.sum(dim=['xgrid','ygrid'])*AREA

#%% Plot selected years
years = ['1988','2000','2015','2016','2017','2019']
x = range(1,13)
plt.figure()
for y in years:
    plt.plot(x,mizext.sel(time=y).values,label=y)
plt.xticks(x,labels=calendar.month_abbr[1:])
plt.legend()
plt.ylabel('Million km$^2$')

#%% Plot intensity timeseries for selected months
plt.figure()
for i in [6,7,8]:
    mizint.sel(time=mizint.time.dt.month.isin([i])).plot(label=months[i-1],marker='o')
plt.legend()
plt.grid(axis='x')

#%% Plot extent timeseries for selected months
plt.figure()
for i in [1,3,11,12]:
    mizext.sel(time=mizext.time.dt.month.isin([i])).plot(label=months[i-1],marker='o')
plt.legend()
plt.grid(axis='x')
plt.ylabel('Million km$^2$')

#%% REGIONAL ANALYSIS (done with GM)
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
dfile=DIR+'seaice_conc_mon_sh_1978_2019_v03r01.nc'
ds = xr.open_dataset(dfile)
miz = ds.goddard_merged_seaice_conc.where(ds.goddard_merged_seaice_conc>=TH)
#%% Compute MIZ extent and intensity for the regions (outputs a dataframe)
KH = [-15,70]   # King Haakon VII
WS = [-60,-15]  # Weddell Sea
EA = [60,170]   # East Antarctica
AB = [-110,-60] # Amundsen/Bellingshausen Seas
RS = [170,-110] # Ross Sea
SH = [-180,180] # Southern Hemisphere
regions180 = [RS,AB,WS,KH,EA,SH]
region_fn180 = ['Ross Sea',
             'Amundsen/Bellingshausen Seas',
             'Weddell Sea',
             'King Haakon VII',
             'East Antarctica',
             'Southern Hemisphere']
region_sn180 = ['RS','AB','WS','KH','EA','SH']

# create a dataframe for the regions (use same index as the dataset)
mizint_reg = miz.time.to_dataframe() # intensity (sigmaSIA)
mizext_reg = miz.time.to_dataframe() # extent (million km2)
T,LON,LAT,I,J = list(miz.coords)

for r,reg in enumerate(regions180):
    print(region_fn180[r])
    if r==0: # for RS region
        mask = ((miz[LON] >= reg[0]) | (miz[LON] < reg[1]))
    else:
        mask = ((miz[LON] >= reg[0]) & (miz[LON] < reg[1]))
    z = miz.where(mask,drop=True)
    intensity = z.mean(dim=['xgrid','ygrid'])
    mizint_reg[region_sn180[r]] = intensity.to_series()
    z = z/z # change to a binary mask
    ext = z.sum(dim=['xgrid','ygrid'])*AREA
    mizext_reg[region_sn180[r]] = ext.to_series()

#%% REGION timeseries
theReg = 'WS'
theMonths = [[12,1,2],[3,4,5],[6,7,8],[9,10,11]]
f = plt.figure(figsize=(9,9))
for p in range(4):
    ax = plt.subplot(4,1,p+1)
    for i in theMonths[p]:
        mizext_reg[theReg].loc[mizext_reg.index.month.isin([i])].plot(ax=ax,
                                                            label=months[i-1],
                                                            marker='o')
    ax.legend()
    plt.grid(axis='x')
    ax.set_ylabel('Million km$^2$')
    ax.set_xlabel('')
f.suptitle(theReg)

#%% climatological SIE NSIDC
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
df=DIR+'seaice_conc_mon_sh_1978_29_v03r.nc'
ds=xr.open_dataset(df)


sic=ds.goddard_merged_seaice_conc.where(ds.goddard_merged_seaice_conc>=0.15)
sie=sic/sic
sie=sie.sum(dim=['xgrid','ygrid'])*AREA
GMsieclim=sie.groupby('time.month').mean(skipna=True)
#GMsieclim_std=sie.groupby('time.month').std()

# Method 1
da=ds.seaice_conc_cdr.sel(time=slice('1988','2019')).groupby('time.month').mean()
sieclim=da.where(da>=0.15)
sieclim=sieclim/sieclim
sieclim=sieclim.sum(dim=['xgrid','ygrid'])*AREA
# Method 2
sic=ds.seaice_conc_cdr.sel(time=slice('1988','2019')).where(ds.seaice_conc_cdr>=0.15)
sie=sic/sic
sie=sie.sum(dim=['xgrid','ygrid'])*AREA
CDRsieclim=sie.groupby('time.month').mean(skipna=True)
#CDRsieclim_std=sie.groupby('time.month').std()

sic=ds.seaice_conc_cdr.sel(time='2017').where(ds.seaice_conc_cdr>=0.15)
sie=sic/sic
sie=sie.sum(dim=['xgrid','ygrid'])*AREA
