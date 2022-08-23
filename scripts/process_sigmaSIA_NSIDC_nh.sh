#!/bin/bash
# chain of processing
#cdo -P 8 -monmean seaice_conc_daily_nh_1979_2019_v04FIX.nc seaice_conc_mon_nh_1979_2019_v04FIX.nc
#cdo -P 8 -ymonmean seaice_conc_mon_nh_1979_2019_v04FIX.nc seaice_conc_monclim_nh_1979_2019_v04FIX.nc
# needs 8 bits for a small value
# -b 8 may be needed here
cdo -b 8 -P 8 -ymonsub seaice_conc_daily_nh_1979_2019_v04FIX.nc seaice_conc_monclim_nh_1979_2019_v04FIX.nc seaice_conc_dayanom_nh_1979_2019_v04FIX.nc

cdo -P 4 -monstd seaice_conc_dayanom_nh_1979_2019_v04FIX.nc seaice_conc_dayanom_monstd_nh_1979_2019_v04FIX.nc
cdo -P 4 -ymonstd seaice_conc_dayanom_nh_1979_2019_v04FIX.nc seaice_conc_dayanomstd_nh_1979_2019_v04FIX.nc

# set attributes
cdo -setattribute,sigmaSIA@standard_name=sea_ice_area_fraction_anomaly,sigmaSIA@long_name="Monthly sigmaSIA (monthly std of daily anomaly) from NOAA/NSIDC CDR Northern Hemisphere Sea Ice Concentration" -chname,cdr_seaice_conc,sigmaSIA -selvar,cdr_seaice_conc -selyear,1987/2019 seaice_conc_dayanom_monstd_nh_1979_2019_v04FIX.nc tmp.nc
cdo -setattribute,comment="This file contains a processed version of the original data. All attributes are maintained for consistency. Author: Marcello Vichi (UCT, MARIS)" tmp.nc NSIDC_cdr_m_sigmaSIA_nh.nc
rm tmp.nc

cdo -setattribute,sigmaSIAclim@standard_name=sea_ice_area_fraction_anomaly,sigmaSIAclim@long_name="Climatological monthly sigmaSIA (climatological 'pooled' monthly std of daily anomaly) from NOAA/NSIDC CDR Northern Hemisphere Sea Ice Concentration" -chname,cdr_seaice_conc,sigmaSIAclim -selvar,cdr_seaice_conc seaice_conc_dayanomstd_nh_1979_2019_v04FIX.nc tmp.nc
cdo -setattribute,comment="This file contains a processed version of the original data. All attributes are maintained for consistency. Author: Marcello Vichi (UCT, MARIS)" tmp.nc NSIDC_cdr_clim_m_sigmaSIA_nh.nc
rm tmp.nc

# process diagnostics
# SIE
cdo -setctomiss,0 -selvar,cdr_seaice_conc seaice_conc_monclim_nh_1979_2019_v04FIX.nc tmp.nc
ncap2 -s "where(cdr_seaice_conc<15 && cdr_seaice_conc>=250) cdr_seaice_conc=0." -s "where(cdr_seaice_conc>=15) cdr_seaice_conc=1." tmp.nc NSIDC_cdr_clim_m_extent_nh.nc
cdo chname,cdr_seaice_conc,extent NSIDC_cdr_clim_m_extent_nh.nc tmp.nc
mv tmp.nc NSIDC_cdr_clim_m_extent_nh.nc
