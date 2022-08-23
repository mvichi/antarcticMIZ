#!/bin/bash

# fixcdr_yearly.sh

# Corrects known issues in CDR v4 files
# This version is for daily aggregate files
# usage: fixcdr_yearly.sh input_filename output_filename

# CF compliance checker:
#   https://pumatest.nerc.ac.uk/cgi-bin/cf-checker.pl

ifn="$1"
if [[ "$ifn" == "" ]]; then
  echo "No input name provided"
  exit
fi

if [[ ! -f $ifn ]]; then
  echo "Input file does not exist: $ifn"
  exit
fi

ofn="$2"
if [[ "$ofn" == "" ]]; then
  echo "No output name provided"
  exit
fi

echo " "
echo "Input filename  (${0}): $ifn"
echo "Output filename (${0}): $ofn"
echo " "

# Checking variable: cdr_seaice_conc
# ERROR: (8.1): must be of type byte, short or int
# Correction: change the datatype from ubyte to byte
# Note: even though the _Unsigned attribute is set, ncatted still converts
#       these flag values to the signed-byte equivalent of -5 through -1
ncap2 -O -s 'cdr_seaice_conc=cdr_seaice_conc.convert(NC_BYTE)' $ifn $ofn
ncatted -O -a flag_values,cdr_seaice_conc,m,b,"251,252,253,254,255" $ofn

# Checking variable: nsidc_bt_seaice_conc
# ERROR: (8.1): must be of type byte, short or int
# Correction: change the datatype from ubyte to byte
# Note: even though the _Unsigned attribute is set, ncatted still converts
#       these flag values to the signed-byte equivalent of -5 through -1
ncap2 -O -s 'nsidc_bt_seaice_conc=nsidc_bt_seaice_conc.convert(NC_BYTE)' $ofn $ofn
ncatted -O -a flag_values,nsidc_bt_seaice_conc,m,b,"251,252,253,254,255" $ofn

# Checking variable: nsidc_nt_seaice_conc
# ERROR: (8.1): must be of type byte, short or int
# Correction: change the datatype from ubyte to byte
# Note: even though the _Unsigned attribute is set, ncatted still converts
#       these flag values to the signed-byte equivalent of -5 through -1
ncap2 -O -s 'nsidc_nt_seaice_conc=nsidc_nt_seaice_conc.convert(NC_BYTE)' $ofn $ofn
ncatted -O -a flag_values,nsidc_nt_seaice_conc,m,b,"251,252,253,254,255" $ofn

# Checking variable: spatial_interpolation_flag
# ERROR: (3.3): Invalid standard_name: interpolation_flag
# Correction: rename the standard_name to status_flag
ncatted -O -a standard_name,spatial_interpolation_flag,o,c,status_flag $ofn

# Checking variable: stdev_of_cdr_seaice_conc
# ERROR: Attribute flag_values of incorrect type (expecting 'Data Variable' type, got 'Numeric' type)
# Corrections: (Note: flag values are not included in this field)
#  remove the flag_values attribute
#  remove the flag_meanings attribute
#  add units attribute
ncatted -O -a flag_values,stdev_of_cdr_seaice_conc,d,, $ofn
ncatted -O -a flag_meanings,stdev_of_cdr_seaice_conc,d,, $ofn
ncatted -O -a units,stdev_of_cdr_seaice_conc,a,c,1 $ofn

# Checking variable: temporal_interpolation_flag
# ERROR: (3.3): Invalid standard_name: interpolation_flag
# Correction: rename the standard_name to status_flag
ncatted -O -a standard_name,temporal_interpolation_flag,o,c,status_flag $ofn

# Note: this line becomes obsolete if the time dimension is renamed
#       from tdim to time as below
# Correct the coordinates attributes
# ncatted -O -a coordinates,,m,c,'tdim y x' $ofn

# Rename the tdim dimension from 'tdim' to 'time'
# Note: renaming the dimension causes the value to be lost,
#       so need to copy the time value, then reapply it after the
#       dimension name has been changed
ncap2 -O -s 'timeval=time' $ofn $ofn
ncrename -O -d tdim,time $ofn
ncap2 -O -s 'time=timeval' $ofn $ofn
ncks -O -x -v timeval $ofn $ofn

# Rename the xgrid variable from 'xgrid' to 'x'
# Note: This needs to be done by copying xgrid to x, then deleting xgrid
#       though the deletion can only occur when there are no refs to xgrid
ncap2 -O -s 'x=xgrid' $ofn $ofn

# Rename the ygrid variable from 'ygrid' to 'y'
# Note: This needs to be done by copying ygrid to y, then deleting ygrid
#       though the deletion can only occur when there are no refs to ygrid
ncap2 -O -s 'y=ygrid' $ofn $ofn

# Correct the coordinates attributes for variables that have this attr
# Note: This needs to occur before xgrid and ygrid can be removed
ncatted -O -a coordinates,,m,c,'time y x' $ofn

# Now that coordinates have been redefined, xgrid and ygrid can be removed
ncks -O -C -x -v xgrid $ofn $ofn
ncks -O -C -x -v ygrid $ofn $ofn

# Remove the cell_methods attribute
ncatted -O -a cell_methods,,d,, $ofn
