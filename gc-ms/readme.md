# GC-MS Quality Control

The expected (target) compounds of the QC mix used in the paper are provided in `compounds.csv`, with the columns `compound_name`, `external_id` (a database identifier, here a ChEBI accession) and `rt_seconds` (or `rt_minutes`).
For other custom compound lists, use `--write-compounds-template` to write an example and adapt for your needs.

The resulting mzQC files contain the PSI-MS CV terms MS:4000238 to MS:4000247 (targeted compound metrics, PSI-MS CV 4.2.5 or later).

Steps to reproduce the paper results:

1) Download the raw data (MTBLS4827) from https://www.ebi.ac.uk/metabolights  
2) convert to mzML using ProteoWizard
3) install Python (we used with Python 3.12, but any later version should be fine as well)
4) install Python dependencies: `pip install -r requirements.txt`
5) run `python py_scripts/GC_MS_pickPeaks_and_match.py -c compounds.csv -o mzQC-metrics *.mzML` to create mzQC metrics for each mzML input file 
6) install R 
7) run the R script `plot_GC-MS.R` (any missing dependency packages are installed automatically) to produce the GC-MS QC plots as shown in the paper (Figure 5)


