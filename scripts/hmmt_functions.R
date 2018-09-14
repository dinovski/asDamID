suppressPackageStartupMessages(library(HMMt))
suppressPackageStartupMessages(library(GenomicRanges))

## HMM functions
HMM <- function(normalized, na_solution) {
  # Run a basic HMM for each data column of a damid data frame.
  if (!na_solution %in% c("NA", "keep", "-")) {
    stop("Unknown na_solution")
  }
  
  # Add index to convert to the original order
  normalized.tmp <- normalized
  normalized.tmp$idx <- 1:nrow(normalized)
  
  normalized.tmp <- normalized.tmp[order(normalized.tmp$chr, normalized.tmp$start), ]
  
  ## HMM
  br <- bridge(normalized.tmp[, 1:4])
  
  ## Flush undesired output to the sink
  sink(file="/dev/null")
  fit <- BaumWelchT(x=br$x, series.length=br$series.length)
  sink()
  
  boundstate <- which.max(fit$mu)
  model <- 0 + (fit$ViterbiPath[br$nonvirtuals] == boundstate)
  model <- ifelse(model == 1, "AD", "iAD")
  
  df.hmm <- cbind(normalized.tmp[, c(1:3, 5)], model)
  
  # assume bins with NA cannot be called LAD
  if (na_solution == "NA") {
    df.hmm$model[is.na(normalized.tmp$score)] <- NA
  } else if (na_solution == "-") {
    df.hmm$model[is.na(normalized.tmp$score)] <- "-"
  }
  
  # Convert into the original order
  df.hmm <- df.hmm[order(df.hmm$idx), ]
  
  df.hmm[, c(1:3, 5)]
}


getHMMRanges <- function(gr, i = 1, score = 1) {
  # Get GRanges bins of a GRanges object with a binary score column
  # remove sequences with NA
  gr <- gr[! is.na(mcols(gr)[, i])]
  
  gr <- gr[mcols(gr)[, i] == score]
  gr <- reduce(gr)
  
  gr
}
