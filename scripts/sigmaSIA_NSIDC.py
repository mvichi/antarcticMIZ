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
import pandas as pd
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
#DIR='/mnt/d/SEAICE/NSIDC-G02202_V4/south/'
DIR='/mnt/c/Users/marce/Documents/WORKS/SCIENCE/antarcticMIZ/data/'
dfile=DIR+'NSIDC_cdr_m_sigmaSIA.nc'
anom=xr.open_dataset(dfile)
sigmaSIA = anom.sigmaSIA.where(anom.sigmaSIA>0.,drop=True)/100.

ancillary = xr.open_dataset(DIR+'G02202-cdr-ancillary-sh.nc')
longitude = ancillary.longitude
latitude = ancillary.latitude

ds = xr.open_dataset(DIR+'NSIDC_cdr_clim_m_sigmaSIA.nc')
sigmaSIAclim = ds.sigmaSIAclim.where(ds.sigmaSIAclim>0.,drop=True)/100.

ds = xr.open_dataset(DIR+'NSIDC_cdr_clim_m_maskMIZ.nc')
maskMIZ = ds.maskMIZ.where(ds.maskMIZ>0.)

# Load the conventional climatological monthly SIC
dfile = DIR+'seaice_conc_monclim_sh_1979_2019_v04FIX.nc'
ds = xr.open_dataset(dfile)
sicMIZ = ds.cdr_seaice_conc.where(ds.cdr_seaice_conc>0.)

# Load the monthly SIC
ds_mc = xr.open_dataset(DIR+'seaice_conc_mon_sh_1979_2019_v04FIX.nc')

# define the CRS (this is the NSIDC stereographic grid)
crs = ccrs.Stereographic(central_latitude=-90, central_longitude=0, 
                         true_scale_latitude=-70)
# The projection keyword determines how the plot will look
map_proj = ccrs.SouthPolarStereo()


#%% SigmaSIA median
# Compute the overall median of the monthly sigmaSIA and the empirical CDF 
# for the spatial distribution of the median
# I excluded the gridpoints with no variability. This is
# preferred to identify a threshold for MIZ variability that is unbiased
# It may enhance the variance in open ocean points that are 
# occasionally ice covered
# an alternative is to include all months with standard deviation == 0,
# which weights more the periods of ice-free ocean. 
q = sigmaSIA.median(dim='time')
p_cdf = np.arange(0.01,1.01,0.01)
quant = q.quantile(p_cdf)

#%% Figure 3
fig = plt.figure(figsize=(4.5,6))
gs = fig.add_gridspec(3, 1)
axmap = fig.add_subplot(gs[:2,0], projection=map_proj)
axcdf = fig.add_subplot(gs[-1, 0])

q.plot.hist(ax=axcdf,bins=np.arange(0,0.6,0.02),
            density=True,histtype='step',color='k')
axcdf.set_title('')
axcdf.set_xlabel('Median of $\sigma_{SIA}$ NOAA/NSIDC CDR [-]')
axcdf.set_ylabel('PDF')
axcdf.set_xlim(-0.02,0.3)
axcdf.set_ylim(0,17.5)
ax2t = axcdf.twinx()
ax2t.plot(quant.values,p_cdf)
ax2t.set_ylabel('CDF')
ax2t.set_yticks(np.arange(0,1.1,0.25))

# This is used to define a threshold of variability that is indicative of 
# the MIZ


#ax.coastlines()
x=q.x
y=q.y
axmap.set_extent([x.min(), x.max(), y.min()+0.7e6, y.max()],
                crs=map_proj)
#ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
axmap.gridlines(draw_labels=True)
cmap=plt.cm.get_cmap('cmo.matter', 5)
q.plot(ax=axmap,transform=crs,
       cmap=cmap,vmin=0,vmax=0.25,
       cbar_kwargs={'orientation':'horizontal',
                    'shrink':0.7,'pad':0.05,
                    'label':'Median of $\sigma_{SIA}$ NOAA/NSIDC CDR [-]'})
fig.text(0.01,0.9,'a',fontsize=14,fontweight='bold')
fig.text(0.01,0.35,'b',fontsize=14,fontweight='bold')


#%% Figure 4. Plot the climatological sigmaSIA for each month
# the plotted variable "field" is set here
field = sigmaSIAclim
label = 'Climatological $\sigma_{SIA}$ NOAA/NSIDC CDR [-]'

cmap=plt.cm.get_cmap('cmo.matter', 10)
p = field.plot(figsize=(8.5,9.5),transform=crs,  # the data's projection
             col='time', col_wrap=3,  # multiplot settings
             aspect=len(field.x) / len(field.y),  # for a sensible figsize
             cmap=cmap, vmin=0, vmax=0.5,
             subplot_kws={'projection': map_proj},
             cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
                'label': label})
for i,ax in enumerate(p.axes.flat):
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
    ax.set_title(months[i])

#%% Evaluate uncertainty
# load the mean monthly standard deviation (this is the spatial std)
std = ds_mc.stdev_of_cdr_seaice_conc.where(ds_mc.stdev_of_cdr_seaice_conc>0.,
                                           drop=True)
# get the indexes of each month
month_idxs=std.groupby('time.month').groups

# Figure S
cmap=plt.cm.get_cmap('cmo.amp', 10)
fig=plt.figure(figsize=(7,9.5))
for m in range(12):
    # compute the ratio between clim sigmaSIA and the median spatial std_dev
    std_max = std.isel(time=month_idxs[m+1]).median(dim='time')
    s = sigmaSIAclim.isel(time=m)
    ratio = std_max/s
    ax = plt.subplot(4,3,m+1,projection=map_proj)
    im = plt.pcolormesh(ratio.x,ratio.y,ratio,transform=crs,
                    vmin=0,vmax=1,cmap=cmap)
    ax.gridlines()
    ax.set_title(months[m])
    ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
ax=fig.add_axes([0.27,0.96,0.5,0.02])
plt.colorbar(im,cax=ax,orientation='horizontal',
             shrink=0.8,label='Ratio spatial vs temporal variability')  
#%% Figure 5
f = plt.figure(figsize=(6,6))
cmap1=plt.cm.get_cmap('cmo.amp', 10)
cmap0=plt.cm.get_cmap('gist_stern_r', 8)
# plot the standard deviation for January and August
stdm = std.isel(time=month_idxs[12]).median(dim='time')
s = sigmaSIAclim.isel(time=11)
ratio = stdm/s
ax = plt.subplot(2,2,1,projection=map_proj)
stdm.plot(ax=ax,transform=crs,
          vmin=0,vmax=0.2,cmap=cmap0,
          cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
                       'label':''})
ax.set_extent([x.min()+0.7e6, x.max(), y.min()+0.7e6, y.max()],
            crs=crs)
ax.gridlines()
ax.set_title(months[11])
ax = plt.subplot(2,2,3,projection=map_proj)
ratio.plot(ax=ax,transform=crs,
          vmin=0,vmax=1,cmap=cmap1,
          cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
                       'label':''})
ax.set_extent([x.min()+0.7e6, x.max(), y.min()+0.7e6, y.max()],
            crs=crs)
ax.gridlines()
ax.set_title('')

stdm = std.isel(time=month_idxs[8]).median(dim='time')
s = sigmaSIAclim.isel(time=7)
ratio = stdm/s
ax = plt.subplot(2,2,2,projection=map_proj)
stdm.plot(ax=ax,transform=crs,
          vmin=0,vmax=0.2,cmap=cmap0,
          cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
                       'label':''})
ax.set_extent([x.min()+0.7e6, x.max(), y.min()+0.7e6, y.max()],
            crs=crs)
ax.gridlines()
ax.set_title(months[7])
ax = plt.subplot(2,2,4,projection=map_proj)
ratio.plot(ax=ax,transform=crs,
          vmin=0,vmax=1,cmap=cmap1,
          cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
                       'label':''})
ax.set_extent([x.min()+0.7e6, x.max(), y.min()+0.7e6, y.max()],
            crs=crs)
ax.gridlines()
ax.set_title('')
f.text(0.13,0.88,'a',fontsize=14,fontweight='bold')
f.text(0.55,0.88,'b',fontsize=14,fontweight='bold')
f.text(0.13,0.46,'c',fontsize=14,fontweight='bold')
f.text(0.55,0.46,'d',fontsize=14,fontweight='bold')




#%% Figure 6: all Novembers from 2008 to 2019

# Panel b (sigmaSIA)
# get the months name
months = calendar.month_name[1:]
#months = calendar.month_abbr[1:]
M=11 # November
M=7 # July
selyear=sigmaSIA.sel(time=slice('2014-01','2019-12'))
field=selyear.sel(time=selyear.time.dt.month.isin([M]))
cmap=plt.cm.get_cmap('cmo.matter', 6)
p = field.plot(figsize=(8.5,4.5),transform=crs,  # the data's projection
             col='time', col_wrap=3,  # multiplot settings
             aspect=len(field.x) / len(field.y),  # for a sensible figsize
             cmap=cmap, vmin=0, vmax=0.3,
             subplot_kws={'projection': map_proj},
             cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
                'label': '$\sigma_{SIA}$ NOAA/NSIDC CDR [-]'})
for i,ax in enumerate(p.axes.flat):
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
    ax.set_title('November '+str(field.time.dt.year.values[i]))
p.fig.text(0.1,0.9,'b',fontsize=14,fontweight='bold')     

# Panel a (MIZ SIC)
sic = ds_mc.cdr_seaice_conc.where(ds_mc.cdr_seaice_conc<=80.)/100.
sic = sic.where(sic > 0.15)
selyear = sic.sel(time=slice('2014-01','2019-12'))
field = selyear.sel(time=selyear.time.dt.month.isin([11]))
cmap = plt.cm.get_cmap('viridis_r',13)
p = field.plot(figsize=(8.5,4.5),transform=crs,  # the data's projection
             col='time', col_wrap=3,  # multiplot settings
             aspect=len(field.x) / len(field.y),  # for a sensible figsize
             cmap=cmap, vmin=0.15, vmax=0.80,
             subplot_kws={'projection': map_proj},
             cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
                'label': 'SIC NOAA/NSIDC CDR [-]','extend':'neither'})
for i,ax in enumerate(p.axes.flat):
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
    ax.set_title('November '+str(field.time.dt.year.values[i]))
p.fig.text(0.1,0.9,'a',fontsize=14,fontweight='bold') 


#%% Figure 8: Mask comparison: SIC (contour) vs sigmaSIA (filled)
x=maskMIZ.x
y=maskMIZ.y

plt.figure(figsize=(7,9.5))
for m in range(12):
    ax = plt.subplot(4,3,m+1,projection=map_proj)
    im=ax.pcolormesh(x, y, maskMIZ.isel(time=m), 
                     transform=crs,cmap=plt.cm.viridis_r)
    ax.set_extent([x.min()+0.7e6, x.max(), y.min()+0.7e6, y.max()],
                ccrs.Stereographic(-90, 0))
    ax.gridlines()
    ax.coastlines()
    ax.contour(sicMIZ.x, sicMIZ.y, sicMIZ.isel(time=m),
               levels=[15.,80.], colors='purple', transform=crs)
    ax.set_title(months[m])
#plt.tight_layout()




#%%  July 2017 study case with drifter
BDIR = '/mnt/c/Users/marce/Documents/WORKS/SCIENCE/antarctic-buoys/data/' # directory of the drifter data

#Read in the buoy drifter data
theBuoy = 'Trident_U1'
drifter = pd.read_csv(BDIR+theBuoy+'.csv',index_col='time',parse_dates=True)
lats = drifter['latitude (deg)']
lons = drifter['longitude (deg)']

# Panel b (sigmaSIA)
# get the months name
map_proj = ccrs.SouthPolarStereo(central_longitude=50)
months = calendar.month_name[1:]
#months = calendar.month_abbr[1:]
M=7 # July
field = sigmaSIA.sel(time=slice('2017-07','2017-12'))
miz_sic = ds_mc.cdr_seaice_conc.sel(time=slice('2017-07','2017-12'))
cmap=plt.cm.get_cmap('cmo.matter', 10)
p = field.plot(figsize=(7,4.5),transform=crs,  # the data's projection
         col='time', col_wrap=3,  # multiplot settings
         sharex=False, sharey=False,
         #aspect=len(field.x) / len(field.y),  # for a sensible figsize
         cmap=cmap, vmin=0, vmax=0.5,
         subplot_kws={'projection': map_proj},
         cbar_kwargs={'shrink':0.4,'pad':0.1,'fraction':0.1,'aspect':25,
            'label': '$\sigma_{SIA}$ NOAA/NSIDC CDR [-]',
            'orientation':'horizontal'})
for i,ax in enumerate(p.axes.flat):
    gl=ax.gridlines(draw_labels=['bottom','left'])
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    ax.coastlines()
    ax.set_extent([30, 60, -65, -58], ccrs.PlateCarree())
    ax.set_title(months[M+i-1]+' 2017')
    ax.plot(lons,lats,color='black',linewidth=1,transform=ccrs.Geodetic())
    start = '2017-'+str(M+i)
    x=lons.loc[start:start]
    y=lats.loc[start:start]
    ax.plot(x,y,color='magenta',linewidth=2,transform=ccrs.Geodetic())
    msic = miz_sic.isel(time=i)
    ax.contour(msic.x,msic.y,msic,transform=crs,
               levels=[15,80],colors='green',linewidths=1)
    
#%% ADDITIONAL PLOTS
# not used in the manuscript
# Plot sigmaSIA for a given year
# This is the data projection, that is their original coordinates
data_crs = ccrs.Stereographic(-90, 0)
# The projection keyword determines how the plot will look
map_proj = ccrs.SouthPolarStereo()
# get the months name
months = calendar.month_name[1:]
#months = calendar.month_abbr[1:]
year='2014'
field=sigmaSIA.sel(time=year)
cmap=plt.cm.get_cmap('cmo.matter', 10)
p = field.plot(figsize=(8.5,9.5),transform=data_crs,  # the data's projection
             col='time', col_wrap=3,  # multiplot settings
             aspect=len(field.x) / len(field.y),  # for a sensible figsize
             cmap=cmap, vmin=0, vmax=0.5,
             subplot_kws={'projection': map_proj},
             cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
                'label': '$\sigma_{SIA}$ NOAA/NSIDC CDR [-]'})
for i,ax in enumerate(p.axes.flat):
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([-180, 180, -90, -52], ccrs.PlateCarree())
    ax.set_title(months[i]+' '+year)
   
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

