#!/bin/bash
# =============================================================================
# preprocess_clm_input.bash
#
# Prepares CLM h1 output files for use with MIMICSplus.
# Handles sites split across two simulation periods (e.g. HIST1 and HIST2),
# copies WATSAT from h0 into h1, renames all files to a single unified prefix
# so MIMICSplus can read them with one clm_data_file path, and concatenates
# a spinup file (1850-1869).
#
# Usage:
#   bash preprocess_clm_input.bash
#
# Edit the CONFIGURATION section below before running.
# =============================================================================

# =============================================================================
# CONFIGURATION — edit these for your setup
# =============================================================================

SITE="Bygland"

# Two simulation periods. Set HIST2_DIR="" if you only have one period.
HIST1_CASE="ERA5L_HIST1_Bygland_century"
HIST2_CASE="ERA5L_HIST2_Bygland_century"

HIST1_DIR="/home/elisacw/mimics_yasso_results/${HIST1_CASE}/lnd/hist"
HIST2_DIR="/home/elisacw/mimics_yasso_results/${HIST2_CASE}/lnd/hist"

# Where to write the processed h1 files (will be created if absent)
OUT_DIR="/home/elisacw/mimicsplus_input/${SITE}/clm_h1"

# Unified prefix for all output h1 files.
# MIMICSplus builds filenames as: <prefix><YEAR>.nc  (e.g. Bygland_hist_all.1901.nc)
# This must match clm_data_file in run.bash (everything up to the year).
UNIFIED_PREFIX="${SITE}_hist_all."

# Spinup years (model uses 1850-1869 by default)
SPINUP_START=1850
SPINUP_END=1869
SPINUP_OUT="${OUT_DIR}/${UNIFIED_PREFIX}for_spinup.${SPINUP_START}-${SPINUP_END}.nc"

# =============================================================================
# SETUP
# =============================================================================

source ~/anaconda3/etc/profile.d/conda.sh
conda activate nco

# Strict mode on after conda activation — conda's own scripts have unbound variables
set -euo pipefail

mkdir -p "${OUT_DIR}"
mkdir -p "${OUT_DIR}/raw"   # archive originals before any changes

echo "=== Preprocessing CLM input for site: ${SITE} ==="
echo "    Output prefix: ${UNIFIED_PREFIX}"

# =============================================================================
# FUNCTION: process one simulation directory
# Copies each h1 file to OUT_DIR with a unified prefix name,
# injects WATSAT from the matching h0, and archives the original.
# Arguments:
#   $1  source directory (lnd/hist/)
#   $2  case name prefix (e.g. ERA5L_HIST1_Bygland_century)
# =============================================================================
process_dir() {
    local src_dir="$1"
    local case_name="$2"

    if [[ ! -d "${src_dir}" ]]; then
        echo "  WARNING: directory not found, skipping: ${src_dir}"
        return
    fi

    echo ""
    echo "--- Processing: ${src_dir} ---"

    for h1_file in "${src_dir}"/${case_name}.clm2.h1.*.nc; do
        [[ -e "${h1_file}" ]] || { echo "  No h1 files found."; return; }

        filename=$(basename "${h1_file}")
        datestamp="${filename#${case_name}.clm2.h1.}"  # → 1850-02-01-00000.nc
        datestamp="${datestamp%.nc}"                    # → 1850-02-01-00000
        year="${datestamp%%-*}"                         # → 1850

        # Final output file uses the unified prefix so all years share one prefix
        out_file="${OUT_DIR}/${UNIFIED_PREFIX}${year}.nc"

        if [[ -f "${out_file}" ]]; then
            echo "  Already exists, skipping: ${UNIFIED_PREFIX}${year}.nc"
            continue
        fi

        echo "  Processing year ${year} ..."
        cp "${h1_file}" "${out_file}"

        # Inject WATSAT from matching h0 file
        h0_file="${src_dir}/${case_name}.clm2.h0.${datestamp}.nc"
        if [[ -f "${h0_file}" ]]; then
            if ncdump -h "${h0_file}" 2>/dev/null | grep -q "WATSAT"; then
                echo "    Adding WATSAT from h0 ..."
                ncks -A -v WATSAT "${h0_file}" "${out_file}"
            else
                echo "    WATSAT not found in h0, skipping injection."
            fi
        else
            echo "    h0 file not found, skipping WATSAT injection."
        fi

        # Archive original
        cp "${h1_file}" "${OUT_DIR}/raw/${filename}"
    done
}

# =============================================================================
# PROCESS BOTH PERIODS
# =============================================================================

process_dir "${HIST1_DIR}" "${HIST1_CASE}"

if [[ -n "${HIST2_DIR}" ]]; then
    process_dir "${HIST2_DIR}" "${HIST2_CASE}"
fi

# =============================================================================
# CONCATENATE SPINUP FILE (1850-1869)
# =============================================================================

echo ""
echo "--- Creating spinup file (${SPINUP_START}-${SPINUP_END}) ---"

spinup_files=()
for year in $(seq "${SPINUP_START}" "${SPINUP_END}"); do
    candidate="${OUT_DIR}/${UNIFIED_PREFIX}${year}.nc"
    if [[ -f "${candidate}" ]]; then
        spinup_files+=("${candidate}")
    else
        echo "  WARNING: missing spinup year ${year}: ${candidate}"
    fi
done

if [[ ${#spinup_files[@]} -eq 0 ]]; then
    echo "  WARNING: No spinup files found. Skipping."
else
    echo "  Concatenating ${#spinup_files[@]} files into ${SPINUP_OUT} ..."
    ncrcat "${spinup_files[@]}" "${SPINUP_OUT}"
    echo "  Spinup file written: ${SPINUP_OUT}"
fi

echo ""
echo "=== Done. Output in: ${OUT_DIR} ==="
echo ""
echo "Set in run.bash:"
echo "  clm_data_file='${OUT_DIR}/${UNIFIED_PREFIX}'"
