#!/usr/bin/env Rscript

######################## Required Packages #####################################
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
################################################################################

################################ Functions #####################################

DE_get_input_type <- function() {
  return("counts")
}

DE_allow_random_effect <- function() {
  return(FALSE)
}

## Code is lifted from
## https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
DE_calculate_dge <- function(
    input_matrix,
    feature_metadata,
    sample_metadata,
    testing_var,
    coef_value,
    coef_barcodes,
    formula,
    method = "singlecell::glmQLFit",
    n_cores = 1,
    verbose=TRUE
) {

  resolution <- strsplit(method, split="::", fixed=T)[[1]][1]
  de_method <- strsplit(method, split="::", fixed=T)[[1]][2]

  if (n_cores > 1) {
    # Re-set options to allow multicore
    old <- options(stringsAsFactors = FALSE,
                   mc.cores=n_cores)
    on.exit(options(old), add = TRUE)
  }

  # First get model - Format data into dds object and fit model
  if (verbose) {
    cat("Fitting the model using EdgeR...\n")
  }
  dge <- edgeR::DGEList(
    counts = input_matrix,
    samples = sample_metadata,
    group = sample_metadata[[testing_var]]
  )

  dge <- edgeR::calcNormFactors(dge)
  design <- model.matrix(formula, sample_metadata)

  if (!limma::is.fullrank(design)) {
    non_est_terms <- limma::nonEstimable(design)
    cat(sprintf(paste("Can not model formula. The following variables are",
                      "linearly correlated with a previous term: %s\n."),
                paste(non_est_terms, collapse = ", ")))
    cat(paste("Ending job silently. Returning an empty dataframe for",
              "this test."))
    return(data.frame())
  }

  dge <- edgeR::estimateDisp(dge, design)
  if (de_method == "glmQLFit") {
    fit <- edgeR::glmQLFit(dge, design)
  } else if (de_method == "glmLRT") {
    fit <- edgeR::glmFit(dge, design)
  }

  if (verbose) {
    cat("Done fitting the model.\n")
    cat(sprintf("Performing Likelihood Ratio Test for the variable %s...\n",
                coef_value))
  }

  if (de_method == "glmQLFit") {
    de_genes <- edgeR::glmQLFTest(fit, coef = coef_value)
  } else if (de_method == "glmLRT") {
    de_genes <- edgeR::glmLRT(fit, coef = coef_value)
  }
  # Sort
  de_genes <- as.data.frame(edgeR::topTags(de_genes, n = Inf, p.value = 1))
  if (verbose) {
    cat("Done performing Likelihood Ratio Test.\n")
  }

  ## Rename edgeR results for downstream analysis
  de_genes <- de_genes %>%
    dplyr::rename(
      log2fc = logFC,
      pvalue = PValue,
      fdr = FDR,
      test_statistic = `F` ## Use test_statistic for downstream anaylsis
    )

  # Add other columns to dataframe
  de_genes$std_err <- NA
  de_genes$gene <- rownames(de_genes)
  de_genes$gene_symbol <- feature_metadata[de_genes$gene, ]$gene_symbol
  de_genes$test_statistic_type <- "f_score"
  
  # Gathering distances
  de_genes$cooks_d_max <- NA
  de_genes$cooks_d_min <- NA
  
  return(de_genes)
}

################################################################################
