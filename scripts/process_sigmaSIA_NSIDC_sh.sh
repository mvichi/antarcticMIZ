#!/bin/bash
# chain of processing
cdo -P 8 -monmean seaice_conc_daily_sh_1979_2019_v04FIX.nc seaice_conc_mon_sh_1979_2019_v04FIX.nc
cdo -P 8 -ymonmean seaice_conc_mon_sh_1979_2019_v04FIX.nc seaice_conc_monclim_sh_1979_2019_v04FIX.nc
# needs -b 8 bits for a small value
cdo -b 8 -P 8 -ymonsub seaice_conc_daily_sh_1979_2019_v04FIX.nc seaice_conc_monclim_sh_1979_2019_v04FIX.nc seaice_conc_dayanom_sh_1979_2019_v04FIX.nc

cdo -P 8 -monstd seaice_conc_dayanom_sh_1979_2019_v04FIX.nc seaice_conc_dayanom_monstd_sh_1979_2019_v04FIX.nc
cdo -P 4 -ymonstd seaice_conc_dayanom_sh_1979_2019_v04FIX.nc seaice_conc_dayanomstd_sh_1979_2019_v04FIX.nc

# set attributes
cdo -setattribute,sigmaSIA@standard_name=sea_ice_area_fraction_anomaly,sigmaSIA@long_name="Monthly sigmaSIA (monthly std of daily anomaly) from NOAA/NSIDC CDR Southern Hemisphere Sea Ice Concentration" -chname,cdr_seaice_conc,sigmaSIA -selvar,latitude,longitude,cdr_seaice_conc -selyear,1987/2019 seaice_conc_dayanom_monstd_sh_1979_2019_v04FIX.nc tmp.nc
cdo -setattribute,comment="This file contains a processed version of the original data. All attributes are maintained for consistency. Author: Marcello Vichi (UCT, MARIS)" tmp.nc NSIDC_cdr_m_sigmaSIA.nc
rm tmp.nc

cdo -setattribute,sigmaSIAclim@standard_name=sea_ice_area_fraction_anomaly,sigmaSIAclim@long_name="Climatological monthly sigmaSIA (climatological 'pooled' monthly std of daily anomaly) from NOAA/NSIDC CDR Southern Hemisphere Sea Ice Concentration" -chname,cdr_seaice_conc,sigmaSIAclim -selvar,latitude,longitude,cdr_seaice_conc seaice_conc_dayanomstd_sh_1979_2019_v04FIX.nc tmp.nc
cdo -setattribute,comment="This file contains a processed version of the original data. All attributes are maintained for consistency. Author: Marcello Vichi (UCT, MARIS)" tmp.nc NSIDC_cdr_m_sigmaSIAclim.nc
rm tmp.nc

# process diagnostics
# SIE
cdo -setctomiss,0 -selvar,cdr_seaice_conc seaice_conc_monclim_sh_1979_2019_v04FIX.nc tmp.nc
ncap2 -s "where(cdr_seaice_conc<15 && cdr_seaice_conc>=250) cdr_seaice_conc=0." -s "where(cdr_seaice_conc>=15) cdr_seaice_conc=1." tmp.nc NSIDC_cdr_clim_m_extent.nc
cdo chname,cdr_seaice_conc,extent NSIDC_cdr_clim_m_extent.nc tmp.nc
mv tmp.nc NSIDC_cdr_clim_m_extent.nc
# maskMIZ
ncap2 -s "where(sigmaSIAclim<=0.1) sigmaSIAclim=0." -s "where(sigmaSIAclim>0.1) sigmaSIAclim=1." NSIDC_cdr_clim_m_sigmaSIA.nc tmp.nc
cdo -setattribute,maskMIZ@standard_name=mask,maskMIZ@long_name="MIZ climatological mask" -chname,sigmaSIAclim,maskMIZ tmp.nc NSIDC_cdr_clim_m_maskMIZ.nc
rm tmp.nc
# maskPACK 
ncap2 -s "where(sigmaSIAclim>0.1) sigmaSIAclim=2." -s "where(sigmaSIAclim<=0.1 && sigmaSIAclim>0.) sigmaSIAclim=1." NSIDC_cdr_m_sigmaSIAclim.nc tmp.nc
cdo -setattribute,maskPACK@standard_name=mask,maskPACK@long_name="Pack ice climatological mask" -chname,sigmaSIAclim,maskPACK tmp.nc NSIDC_cdr_clim_m_maskPACK.nc
rm tmp.nc
