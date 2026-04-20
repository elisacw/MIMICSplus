#!/usr/bin/env bash

# =============================================================================
# run.bash — MIMICSplus run script for Bygland ERA5 data (1850–2025)
# =============================================================================

description='Bygland_ERA5_2025'
namelist_file='options.nml'
spinup_years=10
only_spinup='True'    # set to 'False' once spinup is verified to run historical
CLM_version='new'     # ERA5L output uses new CLM variable names

make
mkdir -p Spinup_values
mkdir -p ./results/${description}
mkdir -p ./results/${description}/nml_options
mkdir -p ./results/${description}/Spinup_values

cp ./src/paramMod.f90    ./results/${description}/parameters_${description}.txt
cp ./src/mycmimMod.f90   ./results/${description}/main_${description}.txt

# Root directory where preprocess_clm_input.bash wrote its output
MIMICS_INPUT="/home/elisacw/mimicsplus_input"

for site in Bygland
do
  echo "Running site: ${site}"
  cp ${namelist_file} ./results/${description}/nml_options/options_${description}_${site}.txt

  # clm_data_file: path prefix up to (not including) the year.
  # preprocess_clm_input.bash wrote files named: Bygland_hist_all.YYYY.nc
  clm_data_file="${MIMICS_INPUT}/${site}/clm_h1/${site}_hist_all."

  # clm_mortality_file: path prefix up to (not including) the year.
  # process_mortality.ipynb wrote files named: mort_Bygland_YYYY.nc
  clm_mortality_file="${MIMICS_INPUT}/${site}/mortality/mort_${site}_"

  # clm_surface_file: single surface data file for this site
  clm_surface_file="${MIMICS_INPUT}/${site}/surfdata_${site}_simyr2000.nc"

  ./run_script_dev \
    $site \
    $description \
    $clm_data_file \
    $clm_mortality_file \
    $clm_surface_file \
    $namelist_file \
    "./results/${description}/" \
    $spinup_years \
    $only_spinup \
    $CLM_version

done

mv Spinup_values/* ./results/${description}/Spinup_values/
rm -r Spinup_values
