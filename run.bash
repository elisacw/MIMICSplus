#!/usr/bin/env bash

#Review these variables before running: 
description='test_case'
namelist_file='options.nml'
spinup_years=10
only_spinup='True' #set false to also get historical files
CLM_version='old' #"old" or "new", "old" uses old variable names.
make
mkdir -p  Spinup_values
mkdir -p ./results/${description}
mkdir -p ./results/${description}/nml_options
mkdir -p ./results/${description}/Spinup_values

cp ./src/paramMod.f90 ./results/${description}/parameters_${description}_${site}.txt
cp ./src/mycmimMod.f90 ./results/${description}/main_${description}_${site}.txt

for site in langtjern
do
  echo ${site}
  cp options.nml ./results/${description}/nml_options/options_${description}_${site}.txt
  #clm_data_file='./test/'${site}'_historical/lnd/hist/'${site}'_historical.clm2.'
  #clm_mortality_file='./test/31464_Hurdal_mortality/mort_'${site}'_'
  #clm_surface_file='./test/surfdata_'${site}'_simyr2000.nc'

  clm_data_file='./fake_data/all/'${site}'_hist_all.'
  clm_mortality_file='./fake_data/mortality/'mort'_'${site}'_'
  clm_surface_file='./fake_data/surfdata_0.9x1.25_hist_16pfts_Irrig_CMIP6_simyr2000_langtjern_c230705.nc'
 ./run_script_dev $site $description $clm_data_file $clm_mortality_file $clm_surface_file $namelist_file "./results/"${description}"/" $spinup_years $only_spinup $CLM_version
done
mv Spinup_values/* ./results/${description}/Spinup_values/
rm -r Spinup_values

