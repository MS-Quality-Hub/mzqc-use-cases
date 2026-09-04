# GC-MS Quality Control

Steps to reproduce the paper results:

1) Download the raw data (MTBLS4827) from https://www.ebi.ac.uk/metabolights  
2) convert to mzML using ProteoWizard
3) run `python py_scripts/GC_MS_pickPeaks_and_match.py -o mzQC-metrics *.mzML` to create mzML metrics for each mzML input file (use `requirements.txt` to install dependencies)
4) run the R script `plot_GC-MS.R` to produce the plots as shown in the paper
