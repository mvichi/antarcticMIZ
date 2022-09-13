# antarcticMIZ

This repository contains code to compute a statistical indicator of sea ice variability 
computed from daily maps of sea ice concentration (SIC) from satellites (Vichi, 2022). 
The $\sigma_{SIA}$ indicator is based on the standard deviation of SIC daily anomalies computed over
the monthly time scales, using a monthly climatology as the baseline.
The code has been applied to both the Southern and the Northern Hemispheres, but it has been 
specifically designed to capture the year-to-year variability of the Antarctic marginal ice zone (MIZ), 
and to expand beyond the classical threshold-based definition (15%<SIC<80%).

This code has been applied to the following data sets:
* NOAA/NSIDC Climate Data Record of Passive Microwave Sea Ice Concentration, Version 4 <https://doi.org/10.7265/efmz-2t65>
* EUMETSAT OSI-450 CDR <https://doi.org/10.15770/EUM_SAF_OSI_0008>

However, as explained in Vichi (2022), MIZ variability is better captured by using satellite data that have not been 
gap-filled. For this reason, the NSIDC-G02202\_V4 product has been re-processed to remove the spatial and temporal 
filters. In addition, there are a few issues in V4 that will likely be solved in a future version.
The NOAA/NSIDC support team made available a shell script to fix these issues.
To produce the same results shown in the paper, the user should follow this workflow:

* Download the annual files from NSIDC-G02202\_V4
* Run fixcdr\_yearly.sh to fix issues
* Remove the time and space interpolation with python script NSIDC-G02202\_V4-V3.py
* Run process\_sigmaSIA\_NSIDC\_sh(nh).sh to compute the sigmaSIA

The following processed climatologies are directly available on Zenodo <https://doi.org/10.5281/zenodo.7077676>:
1. NSIDC_cdr_clim_m_sigmaSIA.nc: Climatological monthly $\sigma_{SIA}$ (climatological 'pooled' monthly std of daily anomaly) from NOAA/NSIDC CDR Southern Hemisphere Sea Ice Concentration
2. NSIDC_cdr_clim_m_maskMIZ.nc: Climatological mask of the MIZ extent using a $\sigma_{SIA}>0.1$ criterion
3. NSIDC_cdr_clim_m_sigmaSIA_exceedance.nc: Climatological probability maps of exceeding $\sigma_{SIA}>0.1$ and $\sigma_{SIA}>0.2$

# References

Vichi, M.: An indicator of sea ice variability for the Antarctic marginal ice zone, The Cryosphere, in press, 2022
