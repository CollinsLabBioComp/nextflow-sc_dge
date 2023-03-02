#!/usr/bin/env Rscript

######################## Required Packages #####################################
suppressPackageStartupMessages(library(MAST))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
################################################################################

################################ Functions #####################################

DE_get_input_type <- function() {
  return("logcp10")
}

DE_allow_random_effect <- function() {
  return(TRUE)
}

DE_calculate_dge <- function(
    input_matrix,
    feature_metadata,
    sample_metadata,
    testing_var,
    coef_value,
    coef_barcodes,
    formula,
    method = "singlecell::bayesglm",
    n_cores = 1,
    verbose = TRUE
) {

  resolution <- strsplit(method, split="::", fixed=T)[[1]][1]
  if (resolution == "pseudobulk") {
    stop(paste("You are trying to run MAST with pseudobulk. MAST was designed",
               "for single-cell data"))
  }
  de_method <- strsplit(method, split="::", fixed=T)[[1]][2]

  if (any(grepl("|", attr(terms(formula), "term.labels"), fixed=T))) {
    if (de_method != "glmer") {
      stop(sprintf(paste(
        "You are trying to run MAST with a random effect using method",
        "`%s`. MAST only runs with a random effect using glmer."
      )))
    }
  }

  if (n_cores > 1) {
    # Re-set options to allow multicore
    old <- options(stringsAsFactors = FALSE, mc.cores=n_cores)
    on.exit(options(old), add = TRUE)
  }

  # First, format data into sca object and fit model
  sca <- MAST::FromMatrix(
    exprsArray = input_matrix,
    cData = sample_metadata,
    fData = feature_metadata
  )

  # https://github.com/RGLab/MAST/issues/107
  ebayes <- TRUE
  if (de_method == "glmer") {
    ebayes <- FALSE
  }

  # MAST::zlm is the hurdle model
  # For more information: https://rdrr.io/bioc/MAST/man/zlm.html
  if (verbose) {
    cat("Fitting the model using MAST:zlm...\n")
  }
  zlm_fit <- MAST::zlm(
    formula,
    sca,
    method = de_method,
    ebayes = ebayes,
    silent = FALSE,
    parallel = TRUE
  )
  if (verbose) {
    cat("Done fitting the model.\n")
  }

  # Check for testing var in model
  if (!(coef_value %in% colnames(zlm_fit@coefD))) {
    cat(sprintf(paste("Returning empty dataframe for term %s",
                      "because it does't exist in the model. It was most",
                      "likely removed because it is linearly associated",
                      "with another term."),
                  coef_value))
    rez <- data.frame()
    return(rez)
  }

  if (verbose) {
    cat(sprintf("Performing Likelihood Ratio Test for the variable %s...\n",
                coef_value))
  }

  # There are two ways to calculate LRT:
  # results <- MAST::lrTest(zlm_fit, coef_value) and
  # summary(zlm_fit, doLRT = coef_value, logFC = TRUE)
  results <- summary(zlm_fit, doLRT = coef_value, logFC = TRUE)
  results_dt <- results$datatable

  # Structure results into interpretable format
  # Here we select for:
  # 1. P-values from the hurdle model... there are three different p-values
  #   from the different models (a) discrete, (b) continuous, and (c) hurdle.
  #   We want the hurdle p-values (which combine the conditionally independent
  #   discrete and continuous models).
  # 2. logFC values and related information
  # 3. Z-scores from Stouffer method, combining the conditionally independent
  #    discrete and continuous models
  fcHurdle <- merge(
    results_dt[
      contrast == coef_value & component == "H",
      .(primerid, `Pr(>Chisq)`)
    ], # hurdle P values
    results_dt[
      contrast == coef_value & component == "logFC",
      .(primerid, coef, ci.hi, ci.lo)
    ], # .(primerid, coef, ci.hi, ci.lo, z) select for logFC z scores
    by = "primerid"
  )
  fcHurdle <- merge(
      fcHurdle,
      results_dt[
        contrast == coef_value & component == "S",
        .(primerid, z)
      ], # select for Stouffer z scores
      by = "primerid"
    )
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, method = "BH")]
  fcHurdle <- merge(fcHurdle, as.data.frame(mcols(sca)), by="primerid")

  # Order by pvalue
  # fcHurdle <- fcHurdle[order(fcHurdle[[`Pr(>Chisq)`]], decreasing = F)]
  if (verbose) {
    cat("Done performing Likelihood Ratio Test.\n")
  }

  # Rename MAST results for downstream analysis
  fcHurdle <- fcHurdle %>%
    dplyr::rename(
      gene = primerid,
      pvalue = `Pr(>Chisq)`,
      log2fc = coef,
      test_statistic = z ## Use test_statistic for downstream anaylsis
    )

  # If model doesn't converge, MAST puts NA for log2FC.
  # See: https://github.com/RGLab/MAST/issues/148
  # Add NA for pvalue so it isn't included in pvalue correction, but leave
  fcHurdle[which(is.na(fcHurdle$log2fc)), "pvalue"] <- NA

  fcHurdle <- fcHurdle %>%
    dplyr::arrange(pvalue) %>%
    as.data.frame()

  # Add other necessary items
  fcHurdle$std_err <- (fcHurdle$ci.hi - fcHurdle$ci.lo) / 3.92 # SE from 95% CI
                                                            # 95% CI is default
  fcHurdle$test_statistic_type <- "z_score"
  
  # Gathering distances
  fcHurdle$cooks_d_max <- NA
  fcHurdle$cooks_d_pval_min <- NA

  return(fcHurdle)
}

################################################################################
