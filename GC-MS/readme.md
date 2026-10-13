# GC-MS Quality Control

Steps to reproduce the paper results:

1) Download the raw data (MTBLS4827) from https://www.ebi.ac.uk/metabolights  
2) convert to mzML using ProteoWizard
3) run `python py_scripts/GC_MS_pickPeaks_and_match.py -c compounds.csv -o mzQC-metrics *.mzML` to create mzQC metrics for each mzML input file (use `requirements.txt` to install dependencies)
4) run the R script `plot_GC-MS.R` to produce the plots as shown in the paper

`compounds.csv` lists the expected (target) compounds of the QC mix used in the paper, with the columns `compound_name`, `external_id` (a database identifier, here a ChEBI accession) and `rt_seconds` (or `rt_minutes`). For other compound lists, `--write-compounds-template` writes an example to adapt.

The mzQC files use the PSI-MS CV terms MS:4000238 to MS:4000247 (targeted compound metrics, PSI-MS CV 4.2.5 or later).
