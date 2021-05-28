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
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
dfile=DIR+'NSIDC_cdr_m_sigmaSIA.nc'
ds=xr.open_dataset(dfile)
miz=ds.sigmaSIA.sel(time=slice('1989','2019')).where(ds.sigmaSIA>=TH)

#%% Compute the mask and the MIZ areal extent from NOAA/NSIDC CDR
mizint=miz.mean(dim=['xgrid','ygrid']) # measure of the intensity of variability
miz=miz/miz
mizext=miz.sum(dim=['xgrid','ygrid'])*AREA
# consolidated pack
pack=ds.sigmaSIA.where(ds.sigmaSIA<TH)
pack=pack/pack
packext=pack.sum(dim=['xgrid','ygrid'])*AREA
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
#%% Compute the mask and the MIZ areal extent from other NSIDC datasets 
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
dfile=DIR+'seaice_conc_dayanom_monstd_sh_1978_2019_v03r01.nc'
ds=xr.open_dataset(dfile)
GM = ds.goddard_merged_seaice_conc.where(ds.goddard_merged_seaice_conc>=TH)
GM = GM/GM
GMext=GM.sum(dim=['xgrid','ygrid'])*AREA
BT  = ds.goddard_bt_seaice_conc.where(ds.goddard_bt_seaice_conc>=TH)
BT = BT/BT
BText=BT.sum(dim=['xgrid','ygrid'])*AREA
packBT  = ds.goddard_bt_seaice_conc.where(ds.goddard_bt_seaice_conc<TH)
packBT = packBT/packBT
packBText=packBT.sum(dim=['xgrid','ygrid'])*AREA
NT  = ds.goddard_nt_seaice_conc.where(ds.goddard_nt_seaice_conc>=TH)
NT = NT/NT
NText=NT.sum(dim=['xgrid','ygrid'])*AREA
packNT  = ds.goddard_nt_seaice_conc.where(ds.goddard_nt_seaice_conc<TH)
packNT = packNT/packNT
packNText=packNT.sum(dim=['xgrid','ygrid'])*AREA


#%% MIZ based on SIC threshold from NSIDC datasets
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
dfile=DIR+'seaice_conc_mon_sh_1978_2019_v03r01.nc'
ds = xr.open_dataset(dfile)
sic = ds.seaice_conc_cdr.where(
    (ds.seaice_conc_cdr>=0.15)&(ds.seaice_conc_cdr<0.80))
sic = sic/sic
sicext = sic.sum(dim=['xgrid','ygrid'])*AREA
sicGM  = ds.goddard_merged_seaice_conc.where(
    (ds.goddard_merged_seaice_conc>=0.15)&(ds.goddard_merged_seaice_conc<0.8))
sicGM = sicGM/sicGM
sicGMext = sicGM.sum(dim=['xgrid','ygrid'])*AREA
sicBT  = ds.goddard_bt_seaice_conc.where(
    (ds.goddard_bt_seaice_conc>=0.15)&(ds.goddard_bt_seaice_conc<0.8))
sicBT = sicBT/sicBT
sicBText = sicBT.sum(dim=['xgrid','ygrid'])*AREA
sicNT  = ds.goddard_nt_seaice_conc.where(
    (ds.goddard_nt_seaice_conc>=0.15)&(ds.goddard_nt_seaice_conc<0.8))
sicNT = sicNT/sicNT
sicNText = sicNT.sum(dim=['xgrid','ygrid'])*AREA

#%% MIZ based on SIC threshold from OSI-450
DIR='/mnt/d/SEAICE/ESA/OSI-450-CDR/'
dfile=DIR+'ice_conc_sh_mon_ease2-250_cdr-v2p0_1979-2015.nc'
ds = xr.open_dataset(dfile)
sicOSI = ds.ice_conc.where((ds.ice_conc>=15.)&(ds.ice_conc<80.))
sicOSI = sicOSI/sicOSI
sicOSIext = sicOSI.sum(dim=['xc','yc'])*AREA
#%% Figure: climatological MIZ seasonality
x = range(1,13)
f,(ax1,ax2) = plt.subplots(1,2,figsize=(10,5),sharey='row',
                           constrained_layout=True)
# GM is identical to BT (not plotted)

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

m=mizext.groupby('time.month').mean()
s=mizext.groupby('time.month').std()
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

#%% pack ice
# Not used in the paper since the sigmaSIA method does not distinguish
# multiyear ice from open ocean (no variability)
m=packext.groupby('time.month').mean()
s=packext.groupby('time.month').std()
plt.plot(x,m,'-',color='C0')
plt.fill_between(x,m+s,m-s,color='C0',alpha=0.2,label='_')
m=packBText.groupby('time.month').mean()
s=packBText.groupby('time.month').std()
plt.plot(x,m,'-',color='C1')
plt.fill_between(x,m+s,m-s,color='C1',alpha=0.2,label='_')
m=packNText.groupby('time.month').mean()
s=packNText.groupby('time.month').std()
plt.plot(x,m,'-',color='C2')
plt.fill_between(x,m+s,m-s,color='C2',alpha=0.2,label='_')
plt.title('Threshold = '+str(TH))
