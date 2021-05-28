#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Exceedance probablity
@author: vichi
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean # using this through the "cmo" variable
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

# get the months name
months = calendar.month_name[1:]
#months = calendar.month_abbr[1:]

# This is the data projection, that is their original coordinates
data_crs = ccrs.Stereographic(-90, 0)
# The projection keyword determines how the plot will look
map_proj = ccrs.SouthPolarStereo()

#%%
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
dfile =DIR+'NSIDC_cdr_clim_m_maskMIZ.nc'
ds = xr.open_dataset(dfile)
maskMIZ = ds.maskMIZ.where(ds.maskMIZ>0.)

#%% monthly timeseries
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
dfile = DIR+'NSIDC_cdr_m_sigmaSIA.nc'
ds = xr.open_dataset(dfile)
# compute median (just used for finding grid locations, I guess)
q = ds.sigmaSIA.where(ds.sigmaSIA>0.).median(dim='time')

#%% exceedance
# timeseries of monthly stdevanom at selected locations where std>0
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
points=[
 [0,[-56,-60,-68]],
 [60,[-60,-63,-67]],
 [120,[-62,-64,-66]],
 [180,[-65,-70,-76]],
 [-40,[-59,-62,-70]]
 ]

def extract_exceedance(lon,lat):
    ind=ll2grid(lon, lat, q.longitude, q.latitude) # find the indices of the location
    ts = ds.sigmaSIA.where(ds.sigmaSIA>0.).isel(ygrid=ind[0],xgrid=ind[1])
    thres = q.isel(ygrid=ind[0],xgrid=ind[1]).values
    data = ts.to_pandas() # convert to pandas for convenience
    sorted_data = data.dropna()
    exceedance = 1.-np.arange(1.,len(sorted_data) + 1.)/len(sorted_data)
    return sorted_data.sort_values(),exceedance,thres 
 
f,axes = plt.subplots(5,1,figsize=(8,9))

for i,ax in enumerate(axes.flat):
    p=points[i]
    lon=p[0]
    lats=p[1]
    for n,l in enumerate(lats):
        tsp,ex,thres = extract_exceedance(lon,l)
        ax.plot(tsp,ex,label='{} ({:g})'.format(l,thres))
        ax.plot([thres,thres],[0,1],'--',color=colors[n])
    ax.set_xlim(0,0.5)
    ax.set_yticks(ticks=np.arange(0,1.1,0.25))
    ax.set_ylabel('Exceedance')
    ax.grid(axis='y')
    ax.legend()
    ax.set_title('Longitude {}'.format(lon))
ax.set_xlabel('$\sigma_{SIA}$ [-]')
plt.tight_layout()

#%% Exceedance function 
def exceedance(data,THRESHOLD):
    #print(f"data: {data.shape}")
    dn = data[np.logical_not(np.isnan(data))]
    if (len(dn)==0):
        ex = np.nan
    else:
        sorted_data = np.sort(dn)
        exceedance = 1.-np.arange(1.,len(dn) + 1.)/(len(dn))
        ind = np.argwhere(sorted_data>=THRESHOLD)
        if (len(ind)==0):
            ex = np.nan
        else:
            ex = exceedance[ind[0]]
    return ex 

#%% Compute exceedance of given threshold 
TH = 0.2
# Use whole timeseries
y = ds.sigmaSIA.where(ds.sigmaSIA>0.)

# This was supposed to work but it doesn't. 
# res = xr.apply_ufunc(exceedance,
#                      y,
#                      exclude_dims=set(('time',)),
#                      input_core_dims=[['time']],
#                      output_core_dims=[['dim0']],
#                      vectorize=True)

# Use a classical loop (quick enough!)
x=y.values
out=x[0,:,:]
I,J=out.shape
for i in range(I):
    for j in range(J):
        out[i,j] = exceedance(x[:,i,j],TH)

#%% plot the annual exceedance
cmap=plt.cm.get_cmap('viridis', 6)
plt.figure()
ax = plt.axes(projection=map_proj)
im=ax.pcolormesh(y.xgrid, y.ygrid, out, 
                 vmin=0,vmax=0.6,cmap=cmap,
                 transform=data_crs)
plt.colorbar(im,orientation='horizontal',extend='max',shrink=0.8,pad=0.05)
ax.set_extent([-180, 180, -90, -53], ccrs.PlateCarree())
ax.gridlines()
ax.coastlines()

#%% compute monthly exceedence for a given threshold
TH = 0.1
N,I,J=y.shape
mon_ex=np.zeros((12,I,J))

# loop over months
for m in range(12):
    y=ds.sigmaSIA.sel(time=ds.time.dt.month.isin([m+1])).where(ds.sigmaSIA>0.)
    x=y.values
    for i in range(I):
        for j in range(J):
            mon_ex[m,i,j] = exceedance(x[:,i,j],TH)

#%% Figure: monthly exceedance
cmap=plt.cm.get_cmap('viridis', 10)
fig=plt.figure(figsize=(7,9.5))
for m in range(12):
    ax = plt.subplot(4,3,m+1,projection=map_proj)
    im=ax.pcolormesh(ds.xgrid, ds.ygrid, mon_ex[m,:,:]*100., 
                 vmin=0,vmax=100,cmap=cmap,
                 transform=data_crs)  
    # tmiz=maskMIZ.isel(time=m).values
    # tmiz[np.isnan(tmiz)] = 0.
    # ax.contour(maskMIZ.xgrid, maskMIZ.ygrid, tmiz,
    #             levels=[1],colors='red', transform=data_crs, linewidths=0.5)
    ax.set_extent([-180, 180, -90, -53], ccrs.PlateCarree())
    ax.gridlines()
    ax.coastlines()
    ax.set_title(months[m])
cb_ax = fig.add_axes([0.2, 0.06, 0.6, 0.02])
cbar = fig.colorbar(im,cax=cb_ax,orientation='horizontal',
                    extend='max',label='Exceedance probability [%]')
