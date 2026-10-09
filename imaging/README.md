# Mass spectrometry imaging QC

 `compute_mzqc.py`  imzML → one mzQC file per run (pyOpenMS) 
 `plot_figures.py`  takes the previously produced mzQC files to reproduce the paper's plot. 

## Data

Download both datasets into `data/`:

**Mouse brain, Orbitrap, centroided:** METASPACE dataset
[`2025-04-27_13h26m30s`](https://metaspace2020.org/dataset/2025-04-27_13h26m30s)

```
data/20220425_MouseBrain_male_212G_382x305_Att30_25um.imzml
data/20220425_MouseBrain_male_212G_382x305_Att30_25um.ibd
```

**Control spheroids, TOF, profile:** PRIDE
[PXD013069](https://www.ebi.ac.uk/pride/archive/projects/PXD013069)

```
data/controlspheroids-nonormalization.imzML
data/controlspheroids-nonormalization.ibd
```

## Running

The environment is pinned in `pyproject.toml`/`uv.lock` (pyOpenMS 3.6.0, pymzqc 1.0.1);
[uv](https://docs.astral.sh/uv/) installs it on first use. `./run.sh` runs everything;
it is equivalent to:

```bash
uv run compute_mzqc.py data/20220425_MouseBrain_male_212G_382x305_Att30_25um.imzml \
    --ref-mz 601.03790 --window-ppm 20 --spectrum-type centroid
uv run compute_mzqc.py data/controlspheroids-nonormalization.imzML \
    --ref-mz 1296.6848 --window-ppm 100 --spectrum-type profile
uv run plot_figures.py out/20220425_MouseBrain_male_212G_382x305_Att30_25um.mzQC \
    out/controlspheroids-nonormalization.mzQC
```

This writes `out/<run>.mzQC` and `figures/imaging_qc.{png,pdf}`.

| dataset | reference ion | window | spectrum type |
|---|---|---|---|
| mouse brain | [4DHB−3H₂O+K]⁺ 601.03790 | ±20 ppm | centroid |
| control spheroids | Angiotensin I [M+H]⁺ 1296.6848 | ±100 ppm | profile |