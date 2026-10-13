USE_PACKAGE <- FALSE

X_AXIS_LABEL       <- "MS run acquisition order"
Y_AXIS_LABEL       <- "proportion of out of control \npeptides"
Y_LIMITS           <- c(-1.4, 1.4)   
Y_LABEL_UPPER      <- "Increase"     
Y_LABEL_LOWER      <- "Decrease"     
Y_LABEL_SIDE       <- "right"        
Y_LABEL_OFFSET     <- 0.75           
Y_LABEL_ANGLE      <- 0              
Y_LABEL_SIZE       <- 14             
Y_LABEL_PAD        <- 0.9            
Y_LABEL_WIDTH      <- 1.4            
PANEL_LABEL_UPPER  <- "Mean"         
PANEL_LABEL_LOWER  <- "Variability"  
PANEL_LABEL_Y      <- 1.25           
PANEL_LABEL_HJUST  <- 0              
PANEL_LABEL_SIZE   <- 4.8


SIGN_BY_DIRECTION <- FALSE


PALETTE <- c(
  "Mean increase"        = "#eb861e", 
  "Mean decrease"        = "#6da7de", 
  "Variability increase" = "#5ea15d",  
  "Variability decrease" = "#d82222",  
  "Change point"         = "black"
)

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
  library(grid)
})

mzQCToMSstatsQC <- function(mzQCfile) {
  if (!file.exists(mzQCfile))
    stop("file not found: ", mzQCfile, "\n  (working directory is ", getwd(),
         ")", call. = FALSE)
  doc <- jsonlite::fromJSON(mzQCfile, simplifyVector = FALSE)
  runs <- doc$mzQC$runQualities
  if (is.null(runs) || length(runs) == 0L)
    stop("no runQualities found in ", mzQCfile)

  out <- list()
  for (i in seq_along(runs)) {
    infile <- runs[[i]]$metadata$inputFiles[[1]]

    acquired <- NA_character_
    for (fp in infile$fileProperties) {
      if (!is.null(fp$accession) && fp$accession == "MS:1000747")
        acquired <- as.character(fp$value)
    }
    if (is.na(acquired) && length(infile$fileProperties) > 0)
      acquired <- as.character(infile$fileProperties[[1]]$value)

    tab <- NULL
    for (qm in runs[[i]]$qualityMetrics) {
      v <- qm$value
      if (is.list(v) && !is.null(names(v)) && length(v) > 1L) { tab <- v; break }
    }
    if (is.null(tab)) {
      message("  ! run ", i, " (", infile$name,
              ") has no table-valued metric - skipped")
      next
    }

    tab <- lapply(tab, function(col) unlist(col, use.names = FALSE))
    df <- as.data.frame(tab, stringsAsFactors = FALSE)

    ren <- c(Peptide = "Precursor", IsotopeLabelType = "Annotations",
             BestRetentionTime = "RT", MaxFwhm = "FWHM", MaxFWHM = "FWHM")
    hit <- names(ren)[names(ren) %in% names(df)]
    names(df)[match(hit, names(df))] <- ren[hit]

    df$AcquiredTime <- acquired
    df$RunName      <- infile$name
    df$QCno         <- i
    out[[length(out) + 1L]] <- df
  }

  data <- do.call(rbind, out)
  if (!"Annotations" %in% names(data)) data$Annotations <- NA
  data$Annotations <- as.character(data$Annotations)   

  metrics <- intersect(c("RT", "FWHM", "TotalArea", "MinStartTime",
                         "MaxEndTime"), names(data))
  data <- data[, c("AcquiredTime", "RunName", "QCno", "Precursor",
                   "Annotations", metrics)]
  data$Precursor <- factor(data$Precursor)
  data
}

readMSstatsQCcsv <- function(path) {
  if (!file.exists(path))
    stop("file not found: ", path, "\n  (working directory is ", getwd(), ")",
         call. = FALSE)
  data <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  ren <- c("Acquired Time" = "AcquiredTime", "Best RT" = "RT",
           "BestRetentionTime" = "RT", "Max FWHM" = "FWHM",
           "MaxFWHM" = "FWHM", "MaxFwhm" = "FWHM",
           "Total Area" = "TotalArea", "Min Start Time" = "MinStartTime",
           "Max End Time" = "MaxEndTime")
  hit <- names(ren)[names(ren) %in% names(data)]
  names(data)[match(hit, names(data))] <- ren[hit]
  if (!"Annotations" %in% names(data)) data$Annotations <- NA
  data$Annotations <- as.character(data$Annotations)

  parsed <- as.POSIXct(data$AcquiredTime, tz = "UTC",
                       tryFormats = c("%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M",
                                      "%m/%d/%y %H:%M:%S", "%m/%d/%y %H:%M"))
  key <- unique(data.frame(AcquiredTime = data$AcquiredTime, t = parsed,
                           stringsAsFactors = FALSE))
  key <- key[order(key$t), ]
  key$QCno <- seq_len(nrow(key))
  data <- merge(data, key[, c("AcquiredTime", "QCno")], by = "AcquiredTime")
  data <- data[order(data$QCno, data$Precursor), ]
  data$Precursor <- factor(data$Precursor)
  rownames(data) <- NULL
  data
}



find_custom_metrics <- function(data) {
  data <- data[, which(colnames(data) == "Annotations"):ncol(data), drop = FALSE]
  nums <- vapply(data, is.numeric, logical(1))
  cols <- colnames(data[, nums, drop = FALSE])
  cols <- cols[seq_len(min(length(cols), 10L))]
  setdiff(cols, "QCno")
}

getMetricData <- function(data, peptide, L, U, metric, normalization,
                          selectMean, selectSD) {
  precursor.data <- data[data$Precursor == peptide, ]
  if (is.null(metric)) return(NULL)
  metricData <- precursor.data[, metric]

  if (normalization) {
    if (is.null(selectMean) && is.null(selectSD)) {
      mu <- mean(metricData[L:U])   
      sd <- sd(metricData[L:U])    
    } else {
      mu <- selectMean
      sd <- selectSD
    }
    if (sd == 0) sd <- 0.0001
    return(scale(metricData[seq_along(metricData)], mu, sd))
  }
  metricData
}

XmR.data.prepare <- function(metricData, L, U, type, selectMean, selectSD) {
  t <- numeric(length(metricData) - 1)
  UCL <- 0; LCL <- 0
  InRangeOutRange <- rep(0, length(metricData))

  for (i in 2:length(metricData)) t[i] <- abs(metricData[i] - metricData[i - 1])
  QCno <- seq_along(metricData)

  if (type == "mean") {
    if (is.null(selectMean) && is.null(selectSD)) {
      UCL <- mean(metricData[L:U]) + 2.66 * sd(t[L:U])
      LCL <- mean(metricData[L:U]) - 2.66 * sd(t[L:U])
    } else {
      UCL <- selectMean + 2.66 * selectSD
      LCL <- selectMean - 2.66 * selectSD
    }
    t <- metricData
  } else if (type == "variability") {
    if (is.null(selectMean) && is.null(selectSD)) {
      UCL <- 3.267 * sd(t[1:L - U])
    } else {
      UCL <- 3.267 * selectSD
    }
    LCL <- 0
  }

  for (i in seq_along(metricData))
    InRangeOutRange[i] <- if (t[i] > LCL && t[i] < UCL) "InRange" else "OutRange"

  data.frame(QCno, IndividualValue = metricData, mR = t, UCL, LCL,
             InRangeOutRange)
}

XmR.River.prepare <- function(data, metric, L, U, type, selectMean, selectSD) {
  y.poz <- rep(0, nrow(data)); y.neg <- rep(0, nrow(data))
  counter <- rep(0, nrow(data))
  precursors <- levels(data$Precursor)

  for (j in seq_along(precursors)) {
    metricData <- getMetricData(data, precursors[j], L = L, U = U,
                               metric = metric, normalization = TRUE,
                               selectMean, selectSD)
    counter[seq_along(metricData)] <- counter[seq_along(metricData)] + 1
    plot.data <- XmR.data.prepare(metricData, L, U, type, selectMean, selectSD)

    sub.poz <- plot.data[plot.data$IndividualValue >= plot.data$UCL, ]
    sub.neg <- plot.data[plot.data$IndividualValue <= plot.data$LCL, ]
    y.poz[sub.poz$QCno] <- y.poz[sub.poz$QCno] + 1
    y.neg[sub.neg$QCno] <- y.neg[sub.neg$QCno] + 1
  }

  max_QCno <- max(which(counter != 0))
  pr.y.poz <- y.poz[1:max_QCno] / counter[1:max_QCno]
  pr.y.neg <- y.neg[1:max_QCno] / counter[1:max_QCno]

  data.frame(
    QCno = rep(1:max_QCno, 2),
    pr.y = c(pr.y.poz, pr.y.neg),
    group = if (type == "mean")
      c(rep("Mean increase", max_QCno), rep("Mean decrease", max_QCno))
    else
      c(rep("Variability increase", max_QCno),
        rep("Variability decrease", max_QCno)),
    metric = rep(metric, max_QCno * 2)
  )
}

XmR.River.DataFrame <- function(data, data.metrics, L, U, listMean, listSD) {
  dat <- NULL
  for (metric in data.metrics) {
    data.1 <- XmR.River.prepare(data, metric, L, U, type = "mean",
                                selectMean = listMean[[metric]],
                                selectSD = listSD[[metric]])
    data.2 <- XmR.River.prepare(data, metric, L, U, type = "variability",
                                selectMean = listMean[[metric]],
                                selectSD = listSD[[metric]])
    if (SIGN_BY_DIRECTION) {
      d1 <- grepl("decrease$", data.1$group)
      d2 <- grepl("decrease$", data.2$group)
      data.1$pr.y[d1] <- -(data.1$pr.y[d1])
      data.2$pr.y[d2] <- -(data.2$pr.y[d2])
    } else {
      data.2$pr.y <- -(data.2$pr.y)
    }
    dat <- rbind(dat, data.1, data.2)
  }
  dat
}

CP.data.prepare <- function(metricData, type) {
  n <- length(metricData)
  Et <- numeric(n - 1); SS <- numeric(n - 1); SST <- numeric(n - 1)

  if (type == "mean") {
    for (i in seq_len(n - 1))
      Et[i] <- (n - i) * (((1 / (n - i)) * sum(metricData[(i + 1):n])) - 0)^2
    QCno <- seq_len(n - 1)
  } else if (type == "variability") {
    for (i in seq_along(metricData)) SS[i] <- metricData[i]^2
    for (i in seq_along(metricData)) {
      SST[i] <- sum(SS[seq_along(metricData)])
      Et[i] <- (SST[i] / 2) - ((n - i + 1) / 2) * log(SST[i] / (n - i + 1)) -
        (n - i + 1) / 2
    }
    QCno <- seq_along(metricData)
  }

  data.frame(QCno, Et, tho.hat = which(Et == max(Et)))
}

get_CP_tho.hat <- function(data, L, U, data.metrics, listMean, listSD) {
  tho.hat <- NULL
  precursors <- levels(data$Precursor)
  for (metric in data.metrics) {
    for (j in seq_along(precursors)) {
      metricData <- getMetricData(data, precursors[j], L, U, metric = metric,
                                 normalization = TRUE,
                                 selectMean = listMean[[metric]],
                                 selectSD = listSD[[metric]])
      mix <- rbind(
        data.frame(tho.hat = CP.data.prepare(metricData, "mean")$tho.hat[1],
                   metric = metric, group = "Individual Value", y = 1.1),
        data.frame(tho.hat = CP.data.prepare(metricData,
                                             "variability")$tho.hat[1],
                   metric = metric, group = "Moving Range", y = -1.1)
      )
      tho.hat <- rbind(tho.hat, mix)
    }
  }
  tho.hat
}

SummaryPlot <- function(data = NULL, L = 1, U = 5, method = "XmR",
                        listMean = NULL, listSD = NULL) {
  if (is.null(data)) return()
  if (!is.data.frame(data)) stop(data)

  data.metrics <- find_custom_metrics(data)
  data.metrics <- data.metrics[!data.metrics %in% c("MinStartTime",
                                                    "MaxEndTime")]
  dat <- XmR.River.DataFrame(data, data.metrics, L, U, listMean, listSD)
  tho.hat.df <- get_CP_tho.hat(data, L, U, data.metrics, listMean, listSD)

  gg <- ggplot(dat)
  gg <- gg + geom_hline(yintercept = 0, alpha = 0.5)
  gg <- gg + geom_smooth(method = "loess", formula = y ~ x,
                         aes(x = QCno, y = pr.y, colour = group, group = group))
  gg <- gg + geom_point(data = tho.hat.df,
                        aes(x = tho.hat, y = y, colour = "Change point"))
  gg <- gg + scale_color_manual(
    breaks = c("Mean increase", "Mean decrease", "Variability increase",
               "Variability decrease", "Change point"),
    values = PALETTE,
    guide = "legend")
  gg <- gg + guides(colour = guide_legend(

    override.aes = list(linetype = c(1, 1, 1, 1, 0),
                        shape = c(NA, NA, NA, NA, 16), fill = NA),
    nrow = 1))
  gg <- gg + facet_wrap(~metric, nrow = 1)

  if (!SIGN_BY_DIRECTION) {
    x_label <- min(dat$QCno, na.rm = TRUE)
    gg <- gg + annotate("text", x = x_label, y = PANEL_LABEL_Y,
                        label = PANEL_LABEL_UPPER, hjust = PANEL_LABEL_HJUST,
                        size = PANEL_LABEL_SIZE)
    gg <- gg + annotate("text", x = x_label, y = -PANEL_LABEL_Y,
                        label = PANEL_LABEL_LOWER, hjust = PANEL_LABEL_HJUST,
                        size = PANEL_LABEL_SIZE)
  }

  gg <- gg + scale_y_continuous(
    expand = c(0, 0), limits = Y_LIMITS,
    breaks = c(1, 0.5, 0, -0.5, -1),
    labels = c("1", "0.5", "0", "0.5", "1"),
    minor_breaks = c(-0.75, -0.25, 0.25, 0.75))
  gg <- gg + labs(x = X_AXIS_LABEL, y = Y_AXIS_LABEL)
  gg <- gg + theme_minimal(base_family = "sans")
  gg <- gg + theme(
    plot.title  = element_text(size = 17, face = "bold",
                               margin = margin(10, 0, 10, 0)),
    axis.text.x = element_text(size = 14, vjust = 0.5, family = "sans",
                               color = "#000000"),
    axis.text.y = element_text(size = 14, hjust = 0.5, family = "sans",
                               color = "#000000"),
    axis.title.y = element_text(size = 15.5, family = "sans",
                                color = "#000000"),
    axis.title.x = element_text(size = 15.5, family = "sans",
                                color = "#000000", hjust = 0.5, vjust = 0, margin = margin(t = 10)),
    legend.text = element_text(size = 15.5),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.spacing.x = unit(0.8, "cm"),
    plot.margin = unit(c(1, 3, 1, 1), "lines"),

    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),

    panel.grid.major = element_line(color = "grey80", linewidth = 0.20),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.08),

    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.8),
    panel.spacing = unit(0.6, "lines"),
    strip.background = element_rect(fill = "white", color = "#000000",
                                    linewidth = 0.8),
    strip.text = element_text(size = 15, family = "sans", color = "#000000"),

    axis.ticks = element_line(color = "#000000", linewidth = 0.5)
  )
  gg
}

add_y_direction_labels <- function(gg) {
  g <- ggplotGrob(gg)
  panels <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
  if (nrow(panels) == 0L) return(g)

  right <- identical(Y_LABEL_SIDE, "right")
  horiz <- (Y_LABEL_ANGLE %% 180) == 0
  gp <- grid::gpar(fontsize = Y_LABEL_SIZE, fontfamily = "sans",
                   col = "#000000")

  if (right) {
    pos <- max(panels$r)
  } else {

    axis_l <- g$layout$l[grepl("^axis-l", g$layout$name)]
    if (length(axis_l) == 0L) return(g)
    pos <- min(axis_l) - 1L
  }

  width <- if (horiz)
    grid::grobWidth(grid::textGrob(c(Y_LABEL_UPPER, Y_LABEL_LOWER), gp = gp)) +
      unit(Y_LABEL_PAD, "lines")
  else
    unit(Y_LABEL_WIDTH, "lines")

  g <- gtable::gtable_add_cols(g, width, pos = pos)
  col <- pos + 1L

  rel <- function(y) (y - Y_LIMITS[1]) / diff(Y_LIMITS)
  labs <- grid::textGrob(
    label = c(Y_LABEL_UPPER, Y_LABEL_LOWER),
    x     = if (!horiz) unit(0.5, "npc")
            else if (right) unit(Y_LABEL_PAD / 2, "lines")
            else unit(1, "npc") - unit(Y_LABEL_PAD / 2, "lines"),
    y     = unit(c(rel(Y_LABEL_OFFSET), rel(-Y_LABEL_OFFSET)), "npc"),
    hjust = if (!horiz) 0.5 else if (right) 0 else 1,
    rot   = Y_LABEL_ANGLE,
    gp    = gp)

  gtable::gtable_add_grob(g, labs, t = min(panels$t), b = max(panels$b),
                          l = col, clip = "off", name = "y-direction-labels")
}

RiverPlot <- function(data = NULL, L = 1, U = 5, method = "XmR",
                      listMean = NULL, listSD = NULL) {
  if (USE_PACKAGE && requireNamespace("MSstatsQC", quietly = TRUE))
    return(MSstatsQC::RiverPlot(data = MSstatsQC::DataProcess(data),
                                L = L, U = U, method = method,
                                listMean = listMean, listSD = listSD))
  SummaryPlot(data, L, U, method = method, listMean = listMean,
              listSD = listSD)
}


main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1)
    stop("usage: Rscript mzqc_riverplot.R <input.mzqc|input.csv> [output.pdf] ",
         "[L] [U]")

  infile  <- args[1]
  outfile <- if (length(args) >= 2) args[2] else "riverplot.pdf"
  L <- if (length(args) >= 3) as.integer(args[3]) else 1L
  U <- if (length(args) >= 4) as.integer(args[4]) else 10L

  data <- if (grepl("\\.(mzqc|json)$", infile, ignore.case = TRUE))
    mzQCToMSstatsQC(infile) else readMSstatsQCcsv(infile)

  message("Loaded ", infile, ": ", length(unique(data$QCno)),
          " acquisitions x ", nlevels(data$Precursor), " peptides")
  message("Metrics: ", paste(sort(setdiff(find_custom_metrics(data),
                                          c("MinStartTime", "MaxEndTime"))),
                             collapse = ", "),
          "   guide set: runs ", L, "-", U)

  gg <- RiverPlot(data = data, L = L, U = U, method = "XmR")
  gg <- add_y_direction_labels(gg)

  png_out <- sub("\\.[^.]+$", "", outfile)
  png_out <- paste0(png_out, ".png")
  outputs <- unique(c(outfile, png_out))
  for (f in outputs) {
    ggsave(f, gg, width = 12.5, height = 4.4, dpi = 300, bg = "white")
    message("Wrote ", f)
  }
}

if (sys.nframe() == 0L) main()