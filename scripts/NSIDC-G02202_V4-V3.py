#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 18:15:55 2022

This script process the NSIDC-G02202_V4 data to remove the filters

@author: vichi
"""

import xarray as xr

DIRIN = '/mnt/d/SEAICE/NSIDC-G02202_V4/north/PROC/'
DIROUT = '/mnt/d/SEAICE/NSIDC-G02202_V4/north/V3/'
#DIRIN = '/mnt/d/SEAICE/NSIDC-G02202_V4/south/PROC/'
#DIROUT = '/mnt/d/SEAICE/NSIDC-G02202_V4/south/V3/'
for y in range(1979,2020):
    filein = 'seaice_conc_daily_nh_'+str(y)+'_v04FIX.nc'
    print(filein)
    fileout = filein # overwrite!
    ds = xr.open_dataset(DIRIN+filein)
    cdr = ds.cdr_seaice_conc
    cdr_std = ds.stdev_of_cdr_seaice_conc
    cdr_nt = ds.nsidc_nt_seaice_conc
    cdr_bt = ds.nsidc_bt_seaice_conc
    flagt = ds.temporal_interpolation_flag
    flags = ds.spatial_interpolation_flag
    cond = ~(flagt>0)&~(flags>0)
    ds['cdr_seaice_conc']=cdr.where(cond,drop=True)
    ds['stdev_of_cdr_seaice_conc']=cdr_std.where(cond,drop=True)
    ds['nsidc_nt_seaice_conc']=cdr_nt.where(cond,drop=True)
    ds['nsidc_bt_seaice_conc']=cdr_bt.where(cond,drop=True)
    ds.to_netcdf(DIROUT+fileout)
