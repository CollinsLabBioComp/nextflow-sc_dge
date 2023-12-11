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
  #    from the different models (a) discrete, (b) continuous, and (c) hurdle.
  #    We want the hurdle p-values (which combine the conditionally independent
  #    discrete and continuous models). This corresponds to component == H.
  #    NOTE: there are no test statistics or coefficients for component H,
  #    just p-values (components D and C have test stats and coef).
  # 2. logFC values and related information for coef.
  #    NOTE: log scale is in the same scale as the input. In this case
  #    NOTE: from the MAST manual...
  #      The log fold change can be small, but the Hurdle p-value small and
  #      significant when the sign of the discrete and continuous model
  #      components are discordant so that the marginal log fold change
  #      cancels out. The large sample sizes present in many single cell
  #      experiments also means that there is substantial power to detect
  #      even small changes.
  # 3. Z-scores from Stouffer method, combining the conditionally independent
  #    discrete and continuous models. This correction to component S.
  #    There is nothing else in component S apart from test statistics (stored
  #    in z column).
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
  # Edit: 12/7/23 - Same goes for std_errs and test_stat
  na_ix <- which(
    is.na(fcHurdle$log2fc) |
      is.na(fcHurdle$test_statistic) |
      is.na(fcHurdle$ci.hi) |
      is.na(fcHurdle$ci.lo)
  )
  fcHurdle[na_ix, "pvalue"] <- NA
  fcHurdle[na_ix, "log2fc"] <- NA
  
  # NOTE: the log scale of MAST is whatever the user passed in
  # https://github.com/theislab/single-cell-tutorial/issues/20
  # https://github.com/RGLab/MAST/issues/112
  # In this case, we prep by default in 011-run_differential_expression.R
  # log1p which is uses exp(1) base. Therefore to get log2FC we need
  # to divide by log(2, base = exp(1))
  for (i in c("ci.hi", "ci.lo", "log2fc")) {
    fcHurdle[[i]] <- fcHurdle[[i]] / log(2, base = exp(1))
  }

  fcHurdle <- fcHurdle %>%
    dplyr::arrange(pvalue) %>%
    as.data.frame()

  # Add other necessary items
  fcHurdle$std_err <- (fcHurdle$ci.hi - fcHurdle$ci.lo) / 3.92 # SE from 95% CI
                                                            # 95% CI is default
  fcHurdle$test_statistic_type <- "z_score"

  return(fcHurdle)
}

################################################################################
