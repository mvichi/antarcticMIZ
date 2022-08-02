#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 22 07:07:37 2021
Area extent
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
mpl.rcParams['font.size'] = 14
#%% Constants
AREA = 25.*25./1.e6 # from pixel to 10^6 km2
TH = 0.1 # sigmaSIA threshold
THP = TH*100.
#%% Load the sigmaSIA for NOAA/NSIDC CDR
DIR='/mnt/c/Users/marce/Documents/WORKS/SCIENCE/antarcticMIZ/data/'
dfile = DIR+'NSIDC_cdr_m_sigmaSIA.nc'
ds = xr.open_dataset(dfile)
miz = ds.sigmaSIA.sel(time=slice('1989','2019')).where(ds.sigmaSIA>=THP)

#%% Compute the mask and the MIZ areal extent from NOAA/NSIDC CDR
mizint = miz.mean(dim=['x','y']) # measure of the intensity of variability
miz = miz/miz # quick way to get a 0/1 mask
mizext = miz.sum(dim=['x','y'])*AREA
# consolidated pack
pack=ds.sigmaSIA.where(ds.sigmaSIA<THP)
pack=pack/pack
packext=pack.sum(dim=['x','y'])*AREA
#%% Compute the mask and the MIZ areal extent from OSI-450
DIR='/mnt/d/SEAICE/ESA/OSI-450-CDR/'
dfile=DIR+'ice_conc_sh_dayanom_monstd_ease2-250_cdr-v2p0_1979-2015.nc'
ds=xr.open_dataset(dfile)
OSI=ds.ice_conc.where(ds.ice_conc>=THP)
OSI=OSI/OSI
OSIext=OSI.sum(dim=['xc','yc'])*AREA
# consolidated pack
OSIpack=ds.ice_conc.where(ds.ice_conc<THP)
OSIpack=OSIpack/OSIpack
OSIpackext=OSIpack.sum(dim=['xc','yc'])*AREA
#%% Compute the mask and the MIZ areal extent from the NSIDC datasets 
DIR='/mnt/c/Users/marce/Documents/WORKS/SCIENCE/antarcticMIZ/data/'
dfile=DIR+'seaice_conc_dayanom_monstd_sh_1979_2019_v04FIX.nc'
ds=xr.open_dataset(dfile)
CDR = ds.cdr_seaice_conc.where(ds.cdr_seaice_conc>=THP)
CDR = CDR/CDR
CDRext = CDR.sum(dim=['x','y'])*AREA
BT = ds.nsidc_bt_seaice_conc.where(ds.nsidc_bt_seaice_conc>=THP)
BT = BT/BT
BText = BT.sum(dim=['x','y'])*AREA
packBT  = ds.nsidc_bt_seaice_conc.where(ds.nsidc_bt_seaice_conc<THP)
packBT = packBT/packBT
packBText=packBT.sum(dim=['x','y'])*AREA
NT  = ds.nsidc_nt_seaice_conc.where(ds.nsidc_nt_seaice_conc>=THP)
NT = NT/NT
NText=NT.sum(dim=['x','y'])*AREA
packNT  = ds.nsidc_nt_seaice_conc.where(ds.nsidc_nt_seaice_conc<THP)
packNT = packNT/packNT
packNText=packNT.sum(dim=['x','y'])*AREA


#%% MIZ based on SIC threshold from NSIDC datasets
DIR='/mnt/c/Users/marce/Documents/WORKS/SCIENCE/antarcticMIZ/data/'
dfile=DIR+'seaice_conc_mon_sh_1979_2019_v04FIX.nc'
ds = xr.open_dataset(dfile)
sic = ds.cdr_seaice_conc.where(
    (ds.cdr_seaice_conc>=15)&(ds.cdr_seaice_conc<80))
sic = sic/sic
sicext = sic.sum(dim=['x','y'])*AREA
sicBT  = ds.nsidc_bt_seaice_conc.where(
    (ds.nsidc_bt_seaice_conc>=15)&(ds.nsidc_bt_seaice_conc<80))
sicBT = sicBT/sicBT
sicBText = sicBT.sum(dim=['x','y'])*AREA
sicNT  = ds.nsidc_nt_seaice_conc.where(
    (ds.nsidc_nt_seaice_conc>=15)&(ds.nsidc_nt_seaice_conc<80))
sicNT = sicNT/sicNT
sicNText = sicNT.sum(dim=['x','y'])*AREA

#%% MIZ based on SIC threshold from OSI-450
DIR='/mnt/d/SEAICE/ESA/OSI-450-CDR/'
dfile=DIR+'ice_conc_sh_mon_ease2-250_cdr-v2p0_1979-2015.nc'
ds = xr.open_dataset(dfile,engine='netcdf4')
sicOSI = ds.ice_conc.where((ds.ice_conc>=15.)&(ds.ice_conc<80.))
sicOSI = sicOSI/sicOSI
sicOSIext = sicOSI.sum(dim=['xc','yc'])*AREA
#%% Figure 7: climatological MIZ seasonality
x = range(1,13)
f,(ax1,ax2) = plt.subplots(1,2,figsize=(10,5),sharey='row',
                           constrained_layout=True)

m=sicext.groupby('time.month').mean()
s=sicext.groupby('time.month').std()
ax1.plot(x,m,'o-',label='NOAA/NSIDC CDR')
ax1.fill_between(x,m+s,m-s,alpha=0.2,label='_')
m=sicBText.groupby('time.month').mean()
s=sicBText.groupby('time.month').std()
ax1.plot(x,m,'s-',label='NOAA/NSIDC BT')
ax1.fill_between(x,m+s,m-s,alpha=0.2,label='_')
m=sicNText.groupby('time.month').mean()
s=sicNText.groupby('time.month').std()
ax1.plot(x,m,'x-',label='NOAA/NSIDC NT')
ax1.fill_between(x,m+s,m-s,alpha=0.2,label='_')
m=sicOSIext.groupby('time.month').mean()
s=sicOSIext.groupby('time.month').std()
ax1.plot(x,m,'x-',label='OSI-450')
ax1.fill_between(x,m+s,m-s,alpha=0.2,label='_')
ax1.set_ylabel('Million km$^2$')
ax1.set_ylim(0,11)
ax1.set_xticks(x)
ax1.set_xticklabels(calendar.month_abbr[1:])
ax1.grid(axis='y')
ax1.text(12,10.5,'a',fontsize=14,fontweight='bold')
#ax1.legend()

m=CDRext.groupby('time.month').mean()
s=CDRext.groupby('time.month').std()
ax2.plot(x,m,'o-',label='NOAA/NSIDC CDR')
ax2.fill_between(x,m+s,m-s,alpha=0.2,label='_')
m=BText.groupby('time.month').mean()
s=BText.groupby('time.month').std()
ax2.plot(x,m,'s-',label='NOAA/NSIDC BT')
ax2.fill_between(x,m+s,m-s,alpha=0.2,label='_')
m=NText.groupby('time.month').mean()
s=NText.groupby('time.month').std()
plt.plot(x,m,'x-',label='NOAA/NSIDC NT')
plt.fill_between(x,m+s,m-s,alpha=0.2,label='_')
m=OSIext.groupby('time.month').mean()
s=OSIext.groupby('time.month').std()
ax2.plot(x,m,'x-',label='OSI-450')
ax2.fill_between(x,m+s,m-s,alpha=0.2,label='_')
#ax2.set_ylabel('Million km$^2$')
ax2.set_ylim(0,11)
ax2.set_xticks(x)
ax2.set_xticklabels(calendar.month_abbr[1:])
ax2.grid(axis='y')
ax2.legend(loc='upper left')
ax2.text(12,10.5,'b',fontsize=14,fontweight='bold')

#%% seasonal cycle for specific years
years = ['2000','2010','2015','2016']
f,axes = plt.subplots(2,2,figsize=(9,9)
                      #,constrained_layout=True
                      )
x = range(1,13)
for i,ax in enumerate(axes.flat):
    ax.plot(x,CDRext.sel(time=years[i]),'o-',label='NOAA/NSIDC CDR')
    ax.plot(x,BText.sel(time=years[i]),'s-',label='NOAA/NSIDC BT')
    ax.plot(x,NText.sel(time=years[i]),'x-',label='NOAA/NSIDC NT')
    if i<3: ax.plot(x,OSIext.sel(time=years[i]),'x-',label='OSI-450')
    ax.set_ylim(0,12);
    ax.set_xticks(x)
    ax.set_xticklabels(calendar.month_abbr[1:],fontsize=10)
    ax.grid(axis='y')
    ax.set_title(years[i])

axes[0,0].legend()

#%% pack ice
# Not used in the paper since the sigmaSIA method does not distinguish
# multiyear ice from open ocean (no variability)
# m=packext.groupby('time.month').mean()
# s=packext.groupby('time.month').std()
# plt.plot(x,m,'-',color='C0')
# plt.fill_between(x,m+s,m-s,color='C0',alpha=0.2,label='_')
# m=packBText.groupby('time.month').mean()
# s=packBText.groupby('time.month').std()
# plt.plot(x,m,'-',color='C1')
# plt.fill_between(x,m+s,m-s,color='C1',alpha=0.2,label='_')
# m=packNText.groupby('time.month').mean()
# s=packNText.groupby('time.month').std()
# plt.plot(x,m,'-',color='C2')
# plt.fill_between(x,m+s,m-s,color='C2',alpha=0.2,label='_')
# plt.title('Threshold = '+str(TH))

#%% Extent based on climatologies
DIR='/mnt/c/Users/marce/Documents/WORKS/SCIENCE/antarcticMIZ/data/'
ds=xr.open_dataset(DIR+'NSIDC_cdr_clim_m_maskMIZ.nc')
maskMIZ=ds.maskMIZ
mizext=maskMIZ.sum(dim=['x','y'])*AREA

dfile=DIR+'seaice_conc_monclim_sh_1979_2019_v04FIX.nc'
ds=xr.open_dataset(dfile)
sicMIZ=ds.cdr_seaice_conc.where((ds.cdr_seaice_conc>=15.)&(ds.cdr_seaice_conc<80.))
sicMIZ=sicMIZ/sicMIZ
sicext=sicMIZ.sum(dim=['x','y'])*AREA
extent = ds.cdr_seaice_conc.where((ds.cdr_seaice_conc>=15.)&(ds.cdr_seaice_conc<250.))
extent = extent/extent
extent = extent.sum(dim=['x','y'])*AREA
#%% Figure 9
x = range(1,13)
plt.figure()
plt.plot(x,mizext,'k-o',label='MIZ SIE ($\sigma_{SIA}$)')
plt.plot(x,sicext,'k-x',label='MIZ SIE (SIC range)')
plt.plot(x,extent,'k--',label='Total SIE')
plt.xticks(x,labels=calendar.month_abbr[1:])
plt.ylabel('Million km$^2$')
plt.legend()
