#!/usr/bin/env bash
# Compute the mzQC files for both datasets and draw the figure (see ../README.md).
#   ./run.sh [input_dir] [output_dir]    (default: ../data/input ../data/output)
printf 'computing mzQC and preparing plots for the imaging datasets\n'
set -euox pipefail
here="$(dirname "$0")"
input_dir="${1:-../data/input}"
output_dir="${2:-../data/output}"

uv run --project "$here" "$here/compute_mzqc.py" \
    "$input_dir/20220425_MouseBrain_male_212G_382x305_Att30_25um.imzml" \
    --ref-mz 601.03790 --window-ppm 20 --spectrum-type centroid \
    -o "$output_dir/20220425_MouseBrain_male_212G_382x305_Att30_25um.mzQC"
uv run --project "$here" "$here/compute_mzqc.py" \
    "$input_dir/controlspheroids-nonormalization.imzML" \
    --ref-mz 1296.6848 --window-ppm 100 --spectrum-type profile \
    -o "$output_dir/controlspheroids-nonormalization.mzQC"
uv run --project "$here" "$here/plot_figures.py" \
    "$output_dir/20220425_MouseBrain_male_212G_382x305_Att30_25um.mzQC" \
    "$output_dir/controlspheroids-nonormalization.mzQC" \
    --out "$output_dir/imaging_qc"
