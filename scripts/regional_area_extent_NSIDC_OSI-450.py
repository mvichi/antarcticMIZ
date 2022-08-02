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
DIR ='/home/vichi/WORKS/SCIENCE/antarcticMIZ/data/'
dfile =DIR+'NSIDC_cdr_m_sigmaSIA.nc'
ds = xr.open_dataset(dfile)
mizsigma =ds.sigmaSIA.sel(time=slice('1987','2019')).where(ds.sigmaSIA>=THP)/100.
mizint = mizsigma.mean(dim=['x','y']) # measure of the intensity of variability

ancillary = xr.open_dataset(DIR+'G02202-cdr-ancillary-sh.nc')
longitude = ancillary.longitude
latitude = ancillary.latitude
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

#%% REGIONAL ANALYSIS 
# regional correlation between sigmaSIA and SIC-BASED EXTENT > 15% 
dfile=DIR+'seaice_conc_mon_sh_1979_2019_v04FIX.nc'
ds = xr.open_dataset(dfile)
sic = ds.cdr_seaice_conc.where(ds.cdr_seaice_conc<250) # remove landmask
sicext = sic.sel(time=slice('1987','2019')).where(sic>=15.)
sicext = sicext/sicext # change to binary mask
sicmiz = sic.sel(time=slice('1987','2019')).where((sic>=15.)&(sic<80.))
sicmiz = sicmiz/sicmiz # change to binary mask

#%% Compute MIZ extent and intensity for the regions (outputs a dataframe)
KH = [-15,70]   # King Haakon VII
WS = [-60,-15]  # Weddell Sea
EA = [60,170]   # East Antarctica
AB = [-110,-60] # Amundsen/Bellingshausen Seas
RS = [170,-110] # Ross Sea
SH = [-180,180] # Southern Hemisphere
region = [RS,AB,WS,KH,EA,SH]
region_fn = ['Ross Sea',
             'Amundsen/Bellingshausen Seas',
             'Weddell Sea',
             'King Haakon VII',
             'East Antarctica',
             'Southern Hemisphere']
region_sn = ['RS','AB','WS','KH','EA','SH']

# create dataframes for the regions (use same index as the dataset)
sicext_reg = mizint.time.to_dataframe() # sic-based extent (million km2)
sicmiz_reg = mizint.time.to_dataframe() # MIZ ic-based extent (million km2) 
mizext_reg = mizint.time.to_dataframe() # sigmaSIA-based MIZ extent 
mizint_reg = mizint.time.to_dataframe() # sigmaSIA-based intensity 

#%%
for r,reg in enumerate(region):
    print(region_fn[r])
    if r==0: # for RS region
        mask = ((longitude >= reg[0]) | (longitude < reg[1]))
    else:
        mask = ((longitude >= reg[0]) & (longitude < reg[1]))
    # sigmaSIA
    z = mizsigma.where(mask)
    intensity = z.mean(dim=['x','y'])
    mizint_reg[region_sn[r]] = intensity.to_series()
    z = z/z # change to a binary mask
    ext = z.sum(dim=['x','y'])*AREA
    mizext_reg[region_sn[r]] = ext.to_series()
    # SIC
    z = sicext.where(mask)
    ext = z.sum(dim=['x','y'])*AREA
    sicext_reg[region_sn[r]] = ext.to_series()
    z = sicmiz.where(mask)
    ext = z.sum(dim=['x','y'])*AREA
    sicmiz_reg[region_sn[r]] = ext.to_series()

#%% Annual mean
sicmiz1Y = sicmiz_reg.resample('1Y').mean()
sicext1Y = sicext_reg.resample('1Y').mean()
mizint1Y = mizint_reg.resample('1Y').mean()
mizext1Y = mizext_reg.resample('1Y').mean()
# reference lines
x=np.arange(0,5,1)
y1 = x
y05 = 0.5*x
y025 = 0.25*x
XLIM = [0,4]
YLIM = [0,2]

#%% Annual max
sicmiz1Y = sicmiz_reg.resample('1Y').max()
sicext1Y = sicext_reg.resample('1Y').max()
mizint1Y = mizint_reg.resample('1Y').mean()
mizext1Y = mizext_reg.resample('1Y').max()
# reference lines
x=np.arange(0,8,1)
y1 = x
y05 = 0.5*x
y025 = 0.25*x
XLIM = [0,7]
YLIM = [0,5]
L1='d'; L2='e'; L3='f'
TITLE='Maximum extent'
#%% Annual min
sicmiz1Y = sicmiz_reg.resample('1Y').min()
sicext1Y = sicext_reg.resample('1Y').min()
mizint1Y = mizint_reg.resample('1Y').mean()
mizext1Y = mizext_reg.resample('1Y').min()
# reference lines
x=np.arange(0,8,1)
y1 = x
y05 = 0.5*x
y025 = 0.25*x
XLIM = [0,2]
YLIM = [0,2]
L1='a'; L2='b'; L3='c'
TITLE='Minimum extent'
#%% Scatter plots
# SIC extent
f,ax=plt.subplots(3,1,figsize=(7,9))
plt.axes(ax[0])
plt.scatter(sicext1Y['RS'],sicmiz1Y['RS'])
plt.scatter(sicext1Y['AB'],sicmiz1Y['AB'])
plt.scatter(sicext1Y['WS'],sicmiz1Y['WS'])
plt.scatter(sicext1Y['KH'],sicmiz1Y['KH'])
plt.scatter(sicext1Y['EA'],sicmiz1Y['EA'])
plt.legend(region_fn[:-1],loc='upper left',fontsize=10)
plt.plot(x,y1,'k')
plt.plot(x,y05,'k--')
plt.plot(x,y025,'k-.')
plt.xlim(XLIM); plt.ylim(YLIM)
plt.ylabel('MIZ extent (SIC) [$10^6$ km$^2]$')

# Scatter plot sigma intensity
plt.axes(ax[1])
plt.scatter(sicext1Y['RS'],mizint1Y['RS'])
plt.scatter(sicext1Y['AB'],mizint1Y['AB'])
plt.scatter(sicext1Y['WS'],mizint1Y['WS'])
plt.scatter(sicext1Y['KH'],mizint1Y['KH'])
plt.scatter(sicext1Y['EA'],mizint1Y['EA'])
plt.xlim(XLIM); plt.ylim([0.18,0.24])
plt.ylabel('Mean annual $\sigma_{SIA}$ (MIZ)')

# Scatter plot sigma extent
plt.axes(ax[2])
plt.scatter(sicext1Y['RS'],mizext1Y['RS'])
plt.scatter(sicext1Y['AB'],mizext1Y['AB'])
plt.scatter(sicext1Y['WS'],mizext1Y['WS'])
plt.scatter(sicext1Y['KH'],mizext1Y['KH'])
plt.scatter(sicext1Y['EA'],mizext1Y['EA'])
plt.plot(x,y1,'k')
plt.plot(x,y05,'k--')
plt.plot(x,y025,'k-.')
plt.xlim(XLIM); plt.ylim(YLIM)
plt.xlabel('SIE [$10^6$ km$^2]$')
plt.ylabel('MIZ extent ($\sigma_{SIA}$) [$10^6$ km$^2]$')

f.text(0.01,0.9,L1,fontsize=14,fontweight='bold')
f.text(0.01,0.6,L2,fontsize=14,fontweight='bold')
f.text(0.01,0.35,L3,fontsize=14,fontweight='bold')
f.text(0.4,0.9,TITLE,fontsize=14,fontweight='bold')
#%% focus on KH region (no apparent relation)
cmap=plt.cm.get_cmap('CMRmap',31)
plt.figure()
im = plt.scatter(sicext1Y['KH'],mizext1Y['KH'],
                 c=sicext1Y.index.year,cmap=cmap,s=200)
plt.colorbar(im)

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

#%% climatological SIE NSIDC: issues with the methodology
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
df=DIR+'seaice_conc_mon_sh_1978_2019_v03r01.nc'
ds=xr.open_dataset(df)


sic=ds.goddard_merged_seaice_conc.where(ds.goddard_merged_seaice_conc>=0.15)
sie=sic/sic
sie=sie.sum(dim=['x','y'])*AREA
GMsieclim=sie.groupby('time.month').mean(skipna=True)
#GMsieclim_std=sie.groupby('time.month').std()

# Method 1
da=ds.seaice_conc_cdr.sel(time=slice('1988','2019')).groupby('time.month').mean(skipna=True)
sieclim=da.where(da>=0.15)
sieclim=sieclim/sieclim
sieclim=sieclim.sum(dim=['x','y'])*AREA
# Method 2
sic=ds.seaice_conc_cdr.sel(time=slice('1988','2019')).where(ds.seaice_conc_cdr>=0.15)
sie=sic/sic
sie=sie.sum(dim=['x','y'])*AREA
CDRsieclim=sie.groupby('time.month').mean(skipna=True)
#CDRsieclim_std=sie.groupby('time.month').std()

sieclim.plot(label='mean+extent')
CDRsieclim.plot(label='extent+mean')
plt.ylabel('Million km$^2$')
plt.legend()
plt.grid()
plt.figure()
d=sieclim-CDRsieclim
d.plot()
plt.ylabel('Million km$^2$')
#%%
sic=ds.seaice_conc_cdr.sel(time='2017').where(ds.seaice_conc_cdr>=0.15)
sie=sic/sic
sie=sie.sum(dim=['x','y'])*AREA
