#!/usr/bin/env bash
# Compute the mzQC files for both datasets and draw the figure (see README.md).
printf 'computing mzQC and preparing plots for the imaging datasets\n'
set -euox pipefail
cd "$(dirname "$0")"

uv run compute_mzqc.py data/20220425_MouseBrain_male_212G_382x305_Att30_25um.imzml \
    --ref-mz 601.03790 --window-ppm 20 --spectrum-type centroid
uv run compute_mzqc.py data/controlspheroids-nonormalization.imzML \
    --ref-mz 1296.6848 --window-ppm 100 --spectrum-type profile
uv run plot_figures.py out/20220425_MouseBrain_male_212G_382x305_Att30_25um.mzQC \
    out/controlspheroids-nonormalization.mzQC
