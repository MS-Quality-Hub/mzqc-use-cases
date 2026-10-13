## Install missing packages from CRAN, then load all of them. They go into the
## personal library, which a fresh R installation does not have yet (and the
## system-wide one is usually not writable).
REQUIRED_PACKAGES = c("ggplot2", "rmzqc", "reshape2")
missing_packages = REQUIRED_PACKAGES[!sapply(REQUIRED_PACKAGES, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  user_lib = strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep)[[1]][1]
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
  ## use the configured CRAN mirror, if any (Rscript has none by default)
  repos = getOption("repos", c(CRAN = "@CRAN@"))
  if (is.na(repos["CRAN"]) || repos["CRAN"] == "@CRAN@") repos["CRAN"] = "https://cloud.r-project.org"
  cat("Installing missing R package(s) into", user_lib, ":", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, lib = user_lib, repos = repos)
}
for (pkg in REQUIRED_PACKAGES) library(pkg, character.only = TRUE)

## The per-compound values are stored as one table metric per run (one column of
## compound names, one of external database identifiers, one of values), see the
## PSI-MS CV terms below.
ACC_DELTA_RT_TABLE  = "MS:4000245"  # targeted compound retention time deviations (metric)
ACC_INTENSITY_TABLE = "MS:4000246"  # targeted compound intensities (metric)
ACC_COMPOUND_NAME   = "MS:4000238"  # targeted compound name (column)
ACC_DELTA_RT        = "MS:4000239"  # targeted compound retention time deviation (column)
ACC_INTENSITY       = "MS:4000240"  # targeted compound intensity (column)

MZQC_DIR = "./mzQC-metrics/"

## Pool the runs of every mzQC file in MZQC_DIR (each file holds one or more runs).
## Files are read in alphabetical order and their runs keep the order they have within
## the file; together that is the acquisition order the plots use as x-axis.
mzqc_files = list.files(MZQC_DIR, pattern = "\\.mzQC$", full.names = TRUE, ignore.case = TRUE)
if (length(mzqc_files) == 0) stop("No .mzQC files found in '", normalizePath(MZQC_DIR, mustWork = FALSE), "'")
runQualities = do.call(c, lapply(mzqc_files, function(f) readMZQC(f)$runQualities))
if (length(runQualities) == 0) stop("The mzQC files in '", MZQC_DIR, "' contain no runQualities")
cat("Read", length(runQualities), "run(s) from", length(mzqc_files), "mzQC file(s) in", MZQC_DIR, "\n")

run_names = sapply(runQualities, function(rq) rq$metadata$inputFiles[[1]]$name)
if (anyDuplicated(run_names)) {
  warning("run(s) present in more than one mzQC file: ",
          paste(unique(run_names[duplicated(run_names)]), collapse = ", "))
}

getMetric = function(runQuality, accession)
{
  hit = Filter(function(qm) qm$accession == accession, runQuality$qualityMetrics)
  if (length(hit) == 0) NULL else hit[[1]]
}

## Combine one table metric of all runs into a data.frame: one row per run (in
## acquisition order, plus an 'acqOrder' column), one '<prefix><compound>' column per
## compound. A compound not detected in a run is missing from that run's table -> NA.
tableMetricToDF = function(runQualities, table_accession, value_accession, prefix)
{
  per_run = lapply(runQualities, function(rq) {
    qm = getMetric(rq, table_accession)
    if (is.null(qm)) return(NULL)
    setNames(unlist(qm$value[[value_accession]]), unlist(qm$value[[ACC_COMPOUND_NAME]]))
  })
  names(per_run) = sapply(runQualities, function(rq) rq$metadata$inputFiles[[1]]$name)
  per_run = per_run[!sapply(per_run, is.null)]
  if (length(per_run) == 0) stop("No run contains the metric ", table_accession)

  compounds = unique(unlist(lapply(per_run, names)))
  df = data.frame(filename = names(per_run),
                  do.call(rbind, lapply(per_run, function(values) values[compounds])),
                  check.names = FALSE, row.names = NULL)
  colnames(df) = c("filename", paste0(prefix, compounds))
  df$acqOrder = 1:nrow(df)
  df
}

## One line per compound over the acquisition order, plus the mean across all
## compounds as a thick black line.
plotOverTime = function(d, prefix, title, y_label, outfile)
{
  idx = grep(prefix, colnames(d), value = TRUE)
  d_long = melt(d[, c(idx, "acqOrder")], id.vars = "acqOrder", measure.vars = idx)
  d_long$analyte = gsub(prefix, "", d_long$variable)
  print(head(d_long))
  all_mean = data.frame("acqOrder" = d$acqOrder, "mean" = apply(d[, idx], 1, function(row) mean(row)))
  p = ggplot(d_long) + geom_line(aes(acqOrder, value, col = analyte)) + ggtitle(title) +
    scale_x_continuous("acquisition order", labels = d$acqOrder, breaks = d$acqOrder) + ylab(y_label) +
    geom_line(data = all_mean, aes(acqOrder, mean), colour = "black", linewidth = 1.2, alpha = 0.5)
  ggsave(outfile, plot = p, scale = 0.5)
  p
}

d_dRT = tableMetricToDF(runQualities, ACC_DELTA_RT_TABLE, ACC_DELTA_RT, "deltaRT_")
d_dRT
plotOverTime(d_dRT, "deltaRT_", "GC-MS RT deviation over time", "ΔRT [sec]", "GC-MS_deltaRT.png")

d_int = tableMetricToDF(runQualities, ACC_INTENSITY_TABLE, ACC_INTENSITY, "intensity_")
d_int
plotOverTime(d_int, "intensity_", "GC-MS intensity over time", "intensity", "GC-MS_intensity.png")
