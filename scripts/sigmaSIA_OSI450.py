#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 09:42:49 2020
Identification of climatological MIZ variability

This script plots the std deviation of the anomaly using
OSI-450 CDR data
The anomaly is obtained from the daily value - MONTHLY climatology.
The monthly standard deviation is computed for each month and 
for a climatological month
@author: mvichi
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean # using this through the "cmo" variable
import pandas as pd
import cartopy.crs as ccrs
import calendar

# This is the data projection, that is their original coordinates
data_crs = ccrs.LambertAzimuthalEqualArea(0,-90)
# The projection keyword determines how the plot will look
map_proj = ccrs.SouthPolarStereo()
# get the months name
months = calendar.month_name[1:]
#months = calendar.month_abbr[1:]
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

#%% Load monthly sigmaSIA, the climatological monthly sigmaSIA 
#   and the climatological mask (threshold = 0.1)
# OSI-450 is in %, hence convertion to fraction is done here
DIR='/mnt/d/SEAICE/ESA/OSI-450-CDR/'
dfile=DIR+'OSI-450_cdr_clim_m_sigmaSIA.nc'
ds=xr.open_dataset(dfile)
OSIclim=ds.sigmaSIAclim.where(ds.sigmaSIAclim>0.)/100.

# Compute the median of the monthly std anomalies
DIR='/mnt/d/SEAICE/ESA/OSI-450-CDR/'
df=DIR+'OSI-450_cdr_m_sigmaSIA.nc'
anom=xr.open_dataset(df)
sigmaSIA = anom.sigmaSIA.where(anom.sigmaSIA>0.)/100.

# exclude gridpoints with no variability. This is
# preferred to identify a threshold for MIZ variability that is unbiased
# It may enhance the variance in open ocean points that are 
# occasionally ice covered
q = sigmaSIA.median(dim='time')
q['xc'] = q.xc*1.e3 # convert coordinates from km to m
q['yc'] = q.yc*1.e3 # convert coordinates from km to m

# compute the empirical CDF for the spatial distribution of the median
p_cdf = np.arange(0,1.001,0.001)
quant = q.quantile(p_cdf)
#%%
fig2 = plt.figure(figsize=(5,7),constrained_layout=True)
gs = fig2.add_gridspec(3, 1)
axmed = fig2.add_subplot(gs[:2,0], projection=map_proj)
f2_ax2 = fig2.add_subplot(gs[-1, 0])
q.plot.hist(ax=f2_ax2,bins=np.arange(0,0.6,0.02),
            density=True,histtype='step',color='k')
f2_ax2.set_title('')
f2_ax2.set_xlabel('Median of $\sigma_{SIA}$ OSI-450 CDR [-]')
f2_ax2.set_ylabel('PDF')
f2_ax2.set_xlim(-0.02,0.3)
ax2t = f2_ax2.twinx()
ax2t.plot(quant.values,p_cdf)
ax2t.set_ylabel('CDF')

# Plot the median
# This is used to define a threshold of variability that is indicative of 
# the MIZ

# fmed=plt.figure(figsize=(6, 6))
# axmed = plt.axes(projection=map_proj)
#ax.coastlines()
x=q.xc
y=q.yc
axmed.set_extent([x.min()+1.8e6, x.max()-1.6e6, y.min()+1.9e6, y.max()-1.4e6],
                data_crs)

axmed.gridlines(draw_labels=True)
cmap=plt.cm.get_cmap('cmo.matter', 5)
q.plot(ax=axmed,transform=data_crs,
       cmap=cmap,vmin=0.,vmax=0.25,
       cbar_kwargs={'orientation':'horizontal',
                    'shrink':0.8,'pad':0.08,
                    'label':'Median of $\sigma_{SIA}$ OSI-450 CDR [-]'})

#axmed.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
#plt.tight_layout(pad=1)
#%% Plot the climatological median sigmaSIA for each month
field = OSIclim/100. # convert to fraction
field['xc'] = field.xc*1.e3 # convert coordinates from km to m
field['yc'] = field.yc*1.e3 # convert coordinates from km to m

cmap=plt.cm.get_cmap('cmo.matter', 10)
p = field.plot(figsize=(8.9,9.5),transform=data_crs,  # the data's projection
             col='time', col_wrap=3,  # multiplot settings
             aspect=len(field.xc) / len(field.yc),  # for a sensible figsize
             cmap=cmap, vmin=0, vmax=0.5,
             subplot_kws={'projection': map_proj},
             cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
                'label': 'Climatological $\sigma_{SIA}$ - OSI-450 CDR [-]'})
for i,ax in enumerate(p.axes.flat):
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
    ax.set_title(months[i])




#%% overlay points
# execute this after the figure is produced
plt.tight_layout(pad=1.15)
# locations for timeseries analyses
points=[
 [0,[-56,-60,-68]],
 [60,[-60,-63,-67]],
 [120,[-62,-64,-66]],
 [180,[-65,-70,-76]],
 [-40,[-59,-62,-70]]
 ]

#[-80,[-66,-69,-71]] Amundsen-Bellingshausen

for p in points:
    lon=p[0]
    lats=p[1]
    for l in lats:
        axmed.plot(lon,l,marker='o',markerfacecolor='None',
                   color='k',transform=ccrs.PlateCarree())

# timeseries of monthly stdevanom at selected locations
def extract_ts(lon,lat):
    ind=ll2grid(lon, lat, q.lon, q.lat) # find the indices of the location
    ts = anom.ice_conc.isel(yc=ind[0],xc=ind[1])/100.
    thres = q.isel(yc=ind[0],xc=ind[1]).values
    return ts.to_pandas(),thres # convert to pandas for slicing, otherwise error
def extract_ts_nonzero(lon,lat):
    ind=ll2grid(lon, lat, q.lon, q.lat) # find the indices of the location
    ts = anom.ice_conc.where(anom.ice_conc>0.).isel(yc=ind[0],xc=ind[1])/100.
    thres = q.isel(yc=ind[0],xc=ind[1]).values
    return ts.to_pandas(),thres # convert to pandas for slicing, otherwise error

#ts.sel(time=slice('2010-01-16','2011-12-16')).plot() # this one does not work
#%% plot the distribution of medians
f=plt.figure()
q.plot.hist()
plt.title('')
plt.xlabel('Median of $\sigma_{SIA}$ OSI-450 CDR [-]')
plt.ylabel('Counts')
plt.xlim(0,0.3)
#%% plot timeseries at different meridians
# markers are plotted when stdev exists (there is sea ice)
# a missing marker indicates no sea ice during that month
f,axes = plt.subplots(5,1,figsize=(12,10))

for i,ax in enumerate(axes.flat):
    p=points[i]
    lon=p[0]
    lats=p[1]
    for l in lats:
        tsp,thres = extract_ts(lon,l)
        tsp.loc['2010-01-01':'2015-01-01'].plot(ax=ax,marker='o',
                                label='{} ({:g})'.format(l,thres))
    ax.set_ylim(0,0.42)
    ax.set_xlabel('')
    ax.set_ylabel('$\sigma_{SIA}$ OSI-450 CDR [-]')
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
ax.set_xlabel('$\sigma_{SIA}$ OSI-450 CDR [-]')
plt.tight_layout()

#%% Point Distributions (box-n-whisker)
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
    ax.set_ylabel('$\sigma_{SIA}$ OSI-450 CDR [-]')
    ax.set_title('Longitude {}'.format(lon))
plt.tight_layout()


