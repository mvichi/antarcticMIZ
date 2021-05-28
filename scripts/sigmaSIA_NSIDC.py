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
# get the months name
months = calendar.month_name[1:]

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
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
dfile=DIR+'NSIDC_cdr_m_sigmaSIA.nc'
anom=xr.open_dataset(dfile)
sigmaSIA = anom.sigmaSIA.where(anom.sigmaSIA>0.)

ds = xr.open_dataset(DIR+'NSIDC_cdr_clim_m_sigmaSIA.nc')
sigmaSIAclim = ds.sigmaSIAclim.where(ds.sigmaSIAclim>0.)

ds=xr.open_dataset(DIR+'NSIDC_cdr_clim_m_maskMIZ.nc')
maskMIZ=ds.maskMIZ.where(ds.maskMIZ>0.)

# Compute the overall median of the monthly sigmaSIA.
# exclude gridpoints with no variability. This is
# preferred to identify a threshold for MIZ variability that is unbiased
# It may enhance the variance in open ocean points that are 
# occasionally ice covered
q = sigmaSIA.median(dim='time')

# the computation below (not used) includes all months with standard deviation = 0,
# which weights more the periods of ice-free ocean. 
# q0 = sigmaSIA.median(dim='time')

# compute the empirical CDF for the spatial distribution of the median
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

#%% Figure 2
fig2 = plt.figure(figsize=(8,7),constrained_layout=True)
gs = fig2.add_gridspec(5, 2)
f2_ax2 = fig2.add_subplot(gs[3:, 0])
# locations for timeseries analyses
points=[
 [-40,[-59,-62,-70]],
 [0,[-56,-60,-68]],
 [60,[-60,-63,-67]],
 [120,[-62,-64,-66]],
 [180,[-65,-70,-76]]
 ]

q.plot.hist(ax=f2_ax2,bins=np.arange(0,0.6,0.02),
            density=True,histtype='step',color='k')
f2_ax2.set_title('')
f2_ax2.set_xlabel('Median of $\sigma_{SIA}$ NOAA/NSIDC CDR [-]')
f2_ax2.set_ylabel('PDF')
f2_ax2.set_xlim(-0.02,0.3)
ax2t = f2_ax2.twinx()
ax2t.plot(quant.values,p_cdf)
ax2t.set_ylabel('CDF')

for i in range(5):
    ax=fig2.add_subplot(gs[i, 1])
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
fig2.text(0.01,0.98,'a',fontsize=14,fontweight='bold')
fig2.text(0.01,0.40,'b',fontsize=14,fontweight='bold')
fig2.text(0.97,0.98,'c',fontsize=14,fontweight='bold')
#%% Figure 2 panel a: sigmaSIA median
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

# NOTE: panels have been merged graphically due to issues in axes formatting

#%% Figure supplementary: distributions of the monthly sigmaSIA
f,axs = plt.subplots(4,3,figsize=(8.5,9.5),sharex=True,sharey=True)
for m,ax in enumerate(axs.flatten()):
    sigmaSIA.sel(time=sigmaSIA.time.dt.month.isin([m+1])).plot.hist(ax=ax,density=True)
    ax.set_title(months[m]+' 1988-2019')
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_ylim(0,12.5)
f.text(0.5,0.07, '$\sigma_{SIA}$', ha="center", va="center",fontsize=14)
f.text(0.05,0.5,'Density', ha="center", va="center", rotation=90,fontsize=14)


#%% The climatological sigmaSIA for each month using GM 
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
dfile=DIR+'seaice_conc_dayanomstd_sh_1978_2019_v03r01.nc'
ds=xr.open_dataset(dfile)
GMclim = ds.goddard_merged_seaice_conc.where(ds.goddard_merged_seaice_conc>0.)

#%% Plot the climatological sigmaSIA for each month
# the plotted variable "field" is set here
field = sigmaSIAclim
label = 'Climatological $\sigma_{SIA}$ NOAA/NSIDC CDR [-]'
# field = GMclim
# label = 'Climatological $\sigma_{SIA}$ NOAA/NSIDC GM [-]'

data_crs = ccrs.Stereographic(-90, 0)
map_proj = ccrs.SouthPolarStereo()

cmap=plt.cm.get_cmap('cmo.matter', 10)
p = field.plot(figsize=(8.5,9.5),transform=data_crs,  # the data's projection
             col='time', col_wrap=3,  # multiplot settings
             aspect=len(field.xgrid) / len(field.ygrid),  # for a sensible figsize
             cmap=cmap, vmin=0, vmax=0.5,
             subplot_kws={'projection': map_proj},
             cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
                'label': label})
for i,ax in enumerate(p.axes.flat):
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
    ax.set_title(months[i])

#%% Load the conventional climatological monthly SIC
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
dfile=DIR+'seaice_conc_monclim_sh_1978_2019_v03r01.nc'
ds=xr.open_dataset(dfile)
sicMIZ=ds.seaice_conc_cdr.where(ds.seaice_conc_cdr>0.)

#%% Figure: Mask comparison: SIC (contour) vs sigmaSIA (filled)
x=q.xgrid
y=q.ygrid

plt.figure(figsize=(7,9.5))
for m in range(12):
    ax = plt.subplot(4,3,m+1,projection=map_proj)
    im=ax.pcolormesh(maskMIZ.xgrid, maskMIZ.ygrid, maskMIZ.isel(time=m), 
                     transform=data_crs,cmap=plt.cm.viridis_r)
    ax.set_extent([x.min()+0.7e6, x.max(), y.min()+0.7e6, y.max()],
                ccrs.Stereographic(-90, 0))
    ax.gridlines()
    ax.coastlines()
    ax.contour(sicMIZ.xgrid, sicMIZ.ygrid, sicMIZ.isel(time=m),
               levels=[0.15,0.80], colors='purple', transform=data_crs)
    ax.set_title(months[m])
#plt.tight_layout()

#%% Plot sigmaSIA for a given year
# This is the data projection, that is their original coordinates
data_crs = ccrs.Stereographic(-90, 0)
# The projection keyword determines how the plot will look
map_proj = ccrs.SouthPolarStereo()
# get the months name
months = calendar.month_name[1:]
#months = calendar.month_abbr[1:]
year='2019'
field=sigmaSIA.sel(time=year)
cmap=plt.cm.get_cmap('cmo.matter', 6)
p = field.plot(figsize=(8.5,9.5),transform=data_crs,  # the data's projection
             col='time', col_wrap=3,  # multiplot settings
             aspect=len(field.xgrid) / len(field.ygrid),  # for a sensible figsize
             cmap=cmap, vmin=0, vmax=0.3,
             subplot_kws={'projection': map_proj},
             cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
                'label': '$\sigma_{SIA}$ NOAA/NSIDC CDR [-]'})
for i,ax in enumerate(p.axes.flat):
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
    ax.set_title(months[i]+' '+year)
    
#%% all novembers

data_crs = ccrs.Stereographic(-90, 0)
# The projection keyword determines how the plot will look
map_proj = ccrs.SouthPolarStereo()
# get the months name
months = calendar.month_name[1:]
#months = calendar.month_abbr[1:]
year='2019'
selyear=sigmaSIA.sel(time=slice('2008-01','2019-12'))
field=selyear.sel(time=selyear.time.dt.month.isin([11]))
cmap=plt.cm.get_cmap('cmo.matter', 6)
p = field.plot(figsize=(8.5,9.5),transform=data_crs,  # the data's projection
             col='time', col_wrap=3,  # multiplot settings
             aspect=len(field.xgrid) / len(field.ygrid),  # for a sensible figsize
             cmap=cmap, vmin=0, vmax=0.3,
             subplot_kws={'projection': map_proj},
             cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
                'label': '$\sigma_{SIA}$ NOAA/NSIDC CDR [-]'})
for i,ax in enumerate(p.axes.flat):
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
    ax.set_title(str(field.time.dt.year.values[i]))