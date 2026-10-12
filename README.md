# mzqc-usecases-manuscript

These scripts reproduce the mzQC files and figures for the use cases in the manuscript "mzQC: a versatile way to communicate quality information for biological mass spectrometry"

- **Code:** this github repository (https://github.com/MS-Quality-Hub/mzqc-use-cases)
- **Data:** zenodo, DOI [10.5281/zenodo.22959481](https://doi.org/10.5281/zenodo.22959481) and each use case has two folders on zenodo:
  - `input/` - files the script read
  - `output/` - mzQC files and figures the script produce

running the code of each use case on its `input/` folder will produce files of its `output/` folder

## setup

download the zenodo data and put it with the code so each use case will be like:

```
<use-case>/
├── code/          #from github
└── data/
    ├── input/     #from zenodo
    └── output/    #created by running related code but also available from zenodo
```

requirements:

- **python 3** with `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`, `openpyxl`, `plotly`
- **uv** for the imaging use case (installs its own pinned environment)
- **R** with `jsonlite`, `ggplot2`, `MASS`


## 1.proteogenomics outliers (figure 1)

outlier detection (LoOP and robust PCA) on 120 DDA runs of M. tuberculosis dataset (MSV000081205 / PXD006843)

**input**
- `Mtb-120-outlier-metrics.xlsx` - QuaMeter IDFree QC metrics for the 120 raw files
- `report_v1.1.6_txt_*_heatmap_*.txt`, `report_v1.1.6_txt_*_filename_sort_*.txt`- PTXQC metrics from MaxQuant results
- `report_v1.1.6_*.html` - PTXQC/QC reports 

**run**
```bash
python proteogenomics-outliers.py ../data/input ../data/output
```

**output**
- `Mtb-120-outlier-metrics.mzQC` - QuaMeter metrics as mzQC
- `Mtb-120-outlier-loop-quameter.mzQC`, `Mtb-120-outlier-loop-ptxqc-noMBR.mzQC`, `Mtb-120-outlier-loop-ptxqc-MBR.mzQC` - LoOP outlier scores 
- `outlier_scores.png`, `attributions_quameter.png`, `timestamps.png` - figure 1



## 2. longitudinal system suitability (figure 2)

MSstatsQC river plot for the CPTAC study 9.1, site 65 system suitability runs

**input**
- `Study_9_1_Site_65.csv` - Skyline derived peptide metrics per run (acquisition time, precursor, retention time, peak area, FWHM, peak start/end)
- `Study9.1_Site65SSS_runs_QCED_.../` - Skyline document the metrics were exported from

**run**
```bash
Rscript longitudinal-system-suitability.R ../data/input/Study_9_1_Site_65.csv ../data/output/longitudinal-system-suitability.png
```

**output**
- `longitudinal-system-suitability.png` - figure 2
- `Study9_1_Site_65.mzqc`, `longitudinal-system-suitability.mzqc` - the peptide level QC metrics


## 3. DIA method evaluation (supplementary figure 1)

window level TIC quartile retention times for 24 DIA runs (MSV000081828 / PASS00589) computed with DIAuditor

**input**
- `DIAMetric-byRun.tsv` - DIAuditor run level summary
- `DIAMetric-byIsolationWindow.tsv` - DIAuditor summary per isolation window

**run**
```bash
Rscript reproduce_TIC_DIAMetric.R ../data/input ../data/output
```

**output**
- `DIAMetric.mzqc` - DIA metrics as mzQC
- `diametric.png` - supplementary figure 1

## 4. repository scale metabolomics (figure 3)

QC metrics for public metabolomics files from GNPS/MassIVE, MetaboLights and Metabolomics workbench

**input**
- `metabolomics_repository_qc.csv` - one row per ms file with acquisition level metrics: MS1/MS2 scan counts per polarity, profile/centroid, RT range, MS1 m/z range, MS2 isolation properties, number of unique precursor m/z

**run**
```bash
#mzQC files -> (one per dataset): --n_rows 0 processes all rows
python mzQC_producer_repository-metabolomics.py \
    --input_file ../data/input/metabolomics_repository_qc.csv \
    --output_dir ../data/output/metabolomics_repository_qc_mzQC \
    --n_rows 0

#figure 3
python mzqc_repository-metabolomics.py
```

**output**
- `metabolomics_repository_qc_mzQC/` - one mzQC file per dataset (example: `MSV000079622.mzqc`, `MTBLS136.mzqc` & `ST002008.mzqc`)
- `figure3.png` - figure 3

## 5. mass spectrometry imaging

QC metrics for two imzML runs, computed with pyOpenMS: per pixel maps of total ion current, base peak fraction and reference m/z error with their quantiles, plus acquired pixel fraction, duplicate pixel coordinate count and pixel grid dimensions per run.

python environment is pinned in `pyproject.toml`/`uv.lock` (pyOpenMS 3.6.0, pymzqc 1.0.1), [uv](https://docs.astral.sh/uv/) installs it on first run.

**input**
- `20220425_MouseBrain_male_212G_382x305_Att30_25um.imzml` / `.ibd` - mouse brain, Orbitrap, centroided, METASPACE [`2025-04-27_13h26m30s`](https://metaspace2020.org/dataset/2025-04-27_13h26m30s)
- `controlspheroids-nonormalization.imzML` / `.ibd` - control spheroids, TOF, profile, PRIDE [PXD013069](https://www.ebi.ac.uk/pride/archive/projects/PXD013069)

**run**
```bash
./run.sh ../data/input ../data/output
```
which runs `compute_mzqc.py` per dataset and `plot_figures.py` to produce png + pdf of the figure.

the reference ion and search window are set per dataset in `run.sh`:

| dataset | reference ion | window | spectrum type |
|---|---|---|---|
| mouse brain | [4DHB−3H₂O+K]⁺ 601.03790 | ±20 ppm | centroid |
| control spheroids | Angiotensin I [M+H]⁺ 1296.6848 | ±100 ppm | profile |

**output**
- `20220425_MouseBrain_male_212G_382x305_Att30_25um.mzQC`, `controlspheroids-nonormalization.mzQC` - imaging QC metrics as mzQC
- `imaging_qc.png`, `imaging_qc.pdf` - Figure 6 of the paper.
