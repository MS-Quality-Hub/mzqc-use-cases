suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
})

#metadata that is not in the TSV reports
MZML_LOCATION    <- "C:\\Research\\MSV000081828-PASS00589-DIA-QC\\mzMLs"
DIAMETRIC_VERSION <- "20250723 beta"

#dia_metric TSV  ->  mzQC

metric <- function(accession, name, value)
  list(accession = accession, name = name, value = value)

run_quality <- function(run, win) {
  num <- function(col) as.numeric(run[[col]])
  int <- function(col) as.integer(run[[col]])
  target <- (win$LoMZ + win$HiMZ) / 2
  offset <- (win$HiMZ - win$LoMZ) / 2

  list(
    metadata = list(
      analysisSoftware = list(list(
        accession = "MS:4000189", name = "DIAMetric",
        value = "Data-Independent Acquisition QC Metric Generator",
        uri = "https://github.com/dtabb73/DIAMetric",
        version = DIAMETRIC_VERSION)),
      inputFiles = list(list(
        fileFormat = list(accession = "MS:1000584", name = "mzML format"),
        fileProperties = list(
          list(accession = "MS:1000031", name = "instrument model",
               value = run$Instrument),
          list(accession = "MS:1000747", name = "completion time",
               value = run$StartTimeStamp)),
        location = paste0(MZML_LOCATION, "\\", run$SourceFile, ".mzML"),
        name = run$SourceFile)),
      label = run$SourceFile),
    qualityMetrics = list(
      metric("MS:1000529", "instrument serial number", run$SerialNumber),
      metric("MS:4000053", "chromatography duration", num("RTDuration")),
      metric("MS:4000059", "number of MS1 spectra", int("mzMLMS1Count")),
      metric("MS:4000060", "number of MS2 spectra", int("mzMLMSnCount")),
      metric("MS:1000800", "MS1 mass resolving power", int("MS1Resolution")),
      metric("MS:4000190",
             "The retention times at which 25%, 50%, and 75% of MS1 TIC have been accumulated",
             c(num("MS1TIC25ileRT"), num("MS1TIC50ileRT"), num("MS1TIC75ileRT"))),
      metric("MS:1000285", "Summed Total Ion Current from MS1",
             num("MS1TotalTIC")),
      metric("MS:4000192",
             "The median cycle time between MS1 measurements",
             num("MS1CycleTime")),
      metric("MS:4000061", "MS1 peak count quartiles",
             c(int("MS1PkCountMin"), int("MS1PkCount25ile"),
               int("MS1PkCount50ile"), int("MS1PkCount75ile"),
               int("MS1PkCountMax"))),
      metric("MS:4000194", "Count of isolation windows in this DIA method",
             int("IsolationWindowCount")),
      metric("MS:4000196",
             "Range of MS/MS measurements per isolation window in this DIA method",
             c(int("CyclesMin"), int("CyclesMax"))),
      metric("MS:1003159",
             "Range of lowest m/z in lowest isolation window to highest m/z of highest isolation window in this DIA method",
             c(num("MZRangeMin"), num("MZRangeMax"))),
      metric("MS:4000195",
             "Range of isolation window widths in this DIA method",
             c(num("IsolationWindowWidthMin"), num("IsolationWidowWidthMax"))),
      metric("MS:4000197",
             "The range of retention times at half TIC among isolation windows",
             c(num("TICMedianRTMin"), num("TICMedianRTMax"))),
      metric("MS:4000198",
             "The range of TIC sums among isolation windows",
             c(num("TotalTICMin"), num("TotalTICMax"))),
      metric("MS:4000199",
             "The range of Peak Count medians among isolation windows",
             c(int("PkCountMedianMin"), int("PkCountMedianMax"))),
      metric("MS:1000827",
             "The target m/z values for each isolation window", target),
      metric("MS:1000828",
             "The lower offset m/z values for each isolation window", offset),
      metric("MS:1000829",
             "The upper offset m/z values for each isolation window", offset),
      metric("MS:1001581",
             "FAIMS compensation voltage for each isolation window",
             as.numeric(win$IonMobility)),
      metric("MS:1000234",
             "Mass Resolving Power for each isolation window",
             as.integer(win$MassResolvingPower)),
      metric("MS:4000060",
             "Number of MS2 Spectra for each isolation window",
             as.integer(win$MSMSCount)),
      metric("MS:4000193",
             "The median cycle time between MS2 measurements for each isolation window",
             as.numeric(win$CycleTimeMedian)),
      metric("MS:4000191",
             "The retention times at which 25%, 50%, and 75% of MS2 TIC have been accumulated for each isolation window",
             unname(as.matrix(win[, c("TIC25ileRT", "TIC50ileRT",
                                      "TIC75ileRT")]))),
      metric("MS:1000285",
             "Summed Total Ion Current from MS2 for each isolation window",
             as.numeric(win$TotalTIC)),
      metric("MS:4000062",
             "The quartiles for peak counts per MS/MS for each isolation window",
             unname(as.matrix(win[, c("PkCountMin", "PkCount25ile",
                                      "PkCount50ile", "PkCount75ile",
                                      "PkCountMax")])))
    )
  )
}

writeMzQC <- function(by_run, by_win, path) {
  runs <- lapply(seq_len(nrow(by_run)), function(i) {
    run <- by_run[i, ]
    win <- by_win[by_win$SourceFile == run$SourceFile, ]
    win <- win[order(win$LoMZ), ]
    run_quality(run, win)
  })

  doc <- list(mzQC = list(
    controlledVocabularies = list(list(
      name = "Proteomics Standards Initiative Mass Spectrometry Ontology",
      uri = "https://github.com/HUPO-PSI/psi-ms-CV/releases/download/v4.1.186/psi-ms.obo",
      version = "4.1.186")),
    creationDate = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
    description = "DIAMetric QC report",
    runQualities = runs,
    version = "1.0.0"))

  writeLines(toJSON(doc, auto_unbox = TRUE, pretty = TRUE, digits = NA),
             path)
  message("Wrote ", path)
  invisible(path)
}

#mzQC  ->  data frame for plotting

get_metric <- function(run, accession) {
  for (qm in run$qualityMetrics)
    if (qm$accession == accession) return(qm$value)
  NULL
}

readMzQCWindows <- function(path) {
  doc <- fromJSON(path, simplifyVector = FALSE)
  do.call(rbind, lapply(doc$mzQC$runQualities, function(run) {
    target <- unlist(get_metric(run, "MS:1000827"))
    tic50  <- vapply(get_metric(run, "MS:4000191"),
                     function(q) as.numeric(q[[2]]), numeric(1))
    data.frame(SourceFile = run$metadata$inputFiles[[1]]$name,
               TargetMZ = target, TIC50ileRT = tic50,
               stringsAsFactors = FALSE)
  }))
}

plotTIC <- function(SWATHs) {
  ggplot(SWATHs, aes(x = TargetMZ, y = TIC50ileRT, group = SourceFile)) +
    geom_line(aes(color = SourceFile), linewidth = 0.6) +
    labs(x = "target m/z of isolation window",
         y = "retention time at which half of TIC was measured") +
    theme_minimal(base_family = "sans") +
    theme(
      legend.position  = "none",
      plot.title       = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.20),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.08),
      panel.border     = element_rect(color = "#000000", fill = NA,
                                      linewidth = 0.8),
      axis.title.x = element_text(family = "sans", color = "#000000"),
      axis.title.y = element_text(family = "sans", color = "#000000"),
      axis.text.x  = element_text(family = "sans", color = "#000000"),
      axis.text.y  = element_text(family = "sans", color = "#000000"),
      axis.ticks   = element_line(color = "#000000", linewidth = 0.5)
    )
}

#main

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  input_dir  <- if (length(args) >= 1) args[1] else "../data/input"
  output_dir <- if (length(args) >= 2) args[2] else "../data/output"

  read_tsv <- function(f) {
    path <- file.path(input_dir, f)
    if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
    header <- strsplit(readLines(path, n = 1L), "\t", fixed = TRUE)[[1]]
    tab <- read.delim(path, header = FALSE, skip = 1L, sep = "\t",
                      quote = "", stringsAsFactors = FALSE)
    tab <- tab[, seq_along(header)]
    names(tab) <- header
    tab
  }
  by_run <- read_tsv("DIAMetric-byRun.tsv")
  by_win <- read_tsv("DIAMetric-byIsolationWindow.tsv")
  message("Loaded ", nrow(by_run), " runs, ", nrow(by_win),
          " isolation windows")

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  mzqc_file <- file.path(output_dir, "DIAMetric.mzqc")
  writeMzQC(by_run, by_win, mzqc_file)

  p <- plotTIC(readMzQCWindows(mzqc_file))
  png_file <- file.path(output_dir, "diametric.png")
  ggsave(png_file, p, width = 6, height = 6, dpi = 300, bg = "white")
  message("Wrote ", png_file)
}

if (sys.nframe() == 0L) main()
