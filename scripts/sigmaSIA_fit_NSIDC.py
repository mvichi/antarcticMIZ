#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 21 07:30:06 2021
Analyse the sigmaSIA distribution and fitting distro
@author: vichi
"""

import pandas as pd
import xarray as xr
import numpy as np
import scipy
from sklearn.preprocessing import StandardScaler
import scipy.stats
import matplotlib.pyplot as plt

#%% Load monthly sigmaSIA
DIR='/mnt/d/SEAICE/NSIDC-G02202_V3/south/'
dfile=DIR+'NSIDC_cdr_m_sigmaSIA.nc'
ds = xr.open_dataset(dfile)
sigmaSIA=ds.sigmaSIA.where(ds.sigmaSIA>0,drop=True)
#%% Distribution fitting
# sample too big for running K-S (see below)
df=pd.DataFrame(sigmaSIA.values.flatten(),columns=['sSIA'])
df=df.dropna()
df.describe()
x = np.arange(0,0.51,0.01)
distribution = 'pareto'
dist = getattr(scipy.stats, distribution)
param = dist.fit(df) # returns b, location and scale
cdf_fit=dist.cdf(x,param[0], loc=param[1], scale=param[2])

# distribution = 'expon'
# dist_exp = getattr(scipy.stats, distribution)
# param_exp = dist_exp.fit(y) # returns b, location and scale
# cdf_fit_exp=dist.cdf(x,*param_exp[:-2], loc=param[-2], scale=param_exp[-1])
#%% Figure 1: sigmaSIA distribution
# use larger bins to eliminate gaps
f,ax = plt.subplots()
h1 = sigmaSIA.plot.hist(ax=ax,bins=np.arange(0,0.5,0.02),
            density=True,histtype='step',color='k')
p_cdf = np.arange(0,1.001,0.001)
quant=sigmaSIA.quantile(p_cdf)
ax.set_title('')
ax.set_xlim([-0.01,0.5])
ax.set_xlabel('$\sigma_{SIA}$ NOAA/NSIDC CDR [-]')
ax.set_ylabel('Density')
ax.grid(axis='x')
ax2=ax.twinx()
h2 = ax2.plot(quant.values,p_cdf,label='CDF')
h3 = ax2.plot(x,cdf_fit,label='CDF fitted')
ax2.set_ylabel('Probability')
ax2.legend(loc='center right')
# Solution for joining two legends (does not work with xarray plot)
# objects = h1 + h2
# labs = [l.get_label() for l in objects]
# ax.legend(objects, labs)

#%% Compute the Kolmogorov-Smirnov test
import random
# normalize and subsample for K-S
y = random.sample(list(df.values), 5000)
sc=StandardScaler()
sc.fit(y) # first fit the estimator
y_std=sc.transform(y)

distribution = 'pareto'
dist = getattr(scipy.stats, distribution)
param = dist.fit(y_std)
p = scipy.stats.kstest(y_std, distribution, args=param)
print(p)
