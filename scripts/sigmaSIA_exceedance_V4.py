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
import cartopy.crs as ccrs
import calendar

# get the months name
months = calendar.month_name[1:]
#months = calendar.month_abbr[1:]

# define the CRS (this is the NSIDC stereographic grid)
data_crs = ccrs.Stereographic(central_latitude=-90, central_longitude=0, 
                         true_scale_latitude=-70)
# The projection keyword determines how the plot will look
map_proj = ccrs.SouthPolarStereo()

#%% Load sigmaSIA indicator 
DIR='/home/vichi/WORKS/SCIENCE/antarcticMIZ/data/'
dfile = DIR+'NSIDC_cdr_m_sigmaSIA.nc'
ds = xr.open_dataset(dfile)
sigmaSIA = ds.sigmaSIA.where(ds.sigmaSIA>0.)/100.
q = sigmaSIA.median(dim='time')
ancillary = xr.open_dataset(DIR+'G02202-cdr-ancillary-sh.nc')
longitude = ancillary.longitude
latitude = ancillary.latitude
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
    ind=ll2grid(lon, lat, longitude, latitude) # find the indices of the location
    ts = ds.sigmaSIA.where(ds.sigmaSIA>0.).isel(y=ind[0],x=ind[1])
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
y = sigmaSIA

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
plt.figure(figsize=(8,8))
ax = plt.axes(projection=map_proj)
im=ax.pcolormesh(y.x, y.y, out, 
                 vmin=0,vmax=0.6,cmap=cmap,
                 transform=data_crs)
plt.colorbar(im,orientation='horizontal',extend='max',shrink=0.6,pad=0.05)
ax.set_extent([-180, 180, -90, -53], ccrs.PlateCarree())
ax.gridlines(draw_labels=True)
ax.coastlines()

#%% compute monthly exceedence for a given threshold
TH = 0.1
N,I,J=sigmaSIA.shape
mon_ex=np.zeros((12,I,J))

# loop over months
for m in range(12):
    y=sigmaSIA.sel(time=ds.time.dt.month.isin([m+1]))
    x=y.values
    for i in range(I):
        for j in range(J):
            mon_ex[m,i,j] = exceedance(x[:,i,j],TH)

#%% Figure: monthly exceedance
cmap=plt.cm.get_cmap('viridis', 10)
fig=plt.figure(figsize=(7,9.5))
for m in range(12):
    ax = plt.subplot(4,3,m+1,projection=map_proj)
    im=ax.pcolormesh(ds.x, ds.y, mon_ex[m,:,:]*100., 
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
                    extend='max',
                    label='Exceedance probability $\sigma_{SIA}\geq0.2$ [%]')


#%% Add the exceedance fields to to the dataset
df = xr.open_dataset('NSIDC_cdr_clim_m_sigmaSIA.nc')
time = df.time.values
x = df.x.values
y = df.y.values
da = xr.DataArray(mon_ex*100,coords=[('time',time),('y',y),('x',x)])
da.attrs['units']='%'
da.attrs['long_name']='Probability of exceeding sigmaSIA>=0.1'
da.attrs['standard_name']='exceedance_probability'
da
df['exceedance01'] = da
df.to_netcdf('NSIDC_cdr_clim_m_sigmaSIA_exceedance.nc')

#%% Figure 10
field=df.exceedance01
p = field.plot(figsize=(7,9),transform=data_crs,  # the data's projection
             col='time', col_wrap=3,  # multiplot settings
             aspect=len(field.x) / len(field.y),  # for a sensible figsize
             cmap=cmap, vmin=0, vmax=100,
             subplot_kws={'projection': map_proj},
             cbar_kwargs={'shrink':0.8,'pad':0.05,'fraction':0.1,'aspect':25,
                          'orientation':'horizontal',
                          'extend':'max',
                          'label': 'Exceedance probability $\sigma_{SIA}\geq0.1$ [%]'})
for i,ax in enumerate(p.axes.flat):
    ax.gridlines()
    ax.coastlines()
    ax.set_extent([field.x.min()+0.7e6, field.x.max(), 
                   field.y.min()+0.7e6, field.y.max()], crs=data_crs)
    #ax.set_extent([-180, 180, -90, -55], ccrs.PlateCarree())
    ax.set_title(months[i])
    