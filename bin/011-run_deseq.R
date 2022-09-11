#!/usr/bin/env Rscript

######################## Required Packages #####################################
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(glmGamPoi))
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

## Not required function
is_interaction <- function(val) {
  return(length(strsplit(val, split=":", fixed=T)[[1]]) > 1)
}

DE_calculate_dge <- function(
    input_matrix,
    feature_metadata,
    sample_metadata,
    testing_var,
    coef_value,
    formula,
    method = "singlecell::glmGamPoi",
    n_cores = 1,
    verbose=TRUE
) {

  resolution <- strsplit(method, split="::", fixed=T)[[1]][1]
  de_method <- strsplit(method, split="::", fixed=T)[[1]][2]

  # First get model - Format data into dds object and fit model
  if (verbose) {
    cat("Fitting the model using DESeq2...\n")
  }

  mm <- model.matrix(
    object = formula,
    data = sample_metadata
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = input_matrix,
    colData = sample_metadata,
    design = mm
  )

  if (verbose) {
    cat("Done fitting the model.\n")
    cat(sprintf("Performing Likelihood Ratio Test for the variable %s...\n",
                coef_value))
  }

  ## Grab design and remove the `test_var`
  des <- design(dds)
  

  reduced_des <- des[,-which(colnames(des) == coef_value), drop=F]

  ## Set parallel cores
  par <- F
  if (n_cores > 1) {
    par <- T
  }

  ## Set default DESeq2 parameters to bulk (or pseudobulk).
  # These are the DESeq2 defaults
  sfType <- "poscounts" # "ratio", "poscounts", or "iterate"
  minmu <- 0.5 # lower bound on the estimated count for gene-wise dispersion
  minReplicatesForReplace <- 7 # default is 7
  if (resolution == "singlecell") {
    ## Change default DESeq2 parameters if single cell
    # https://github.com/mikelove/DESeq2/pull/24
    # Defaults are taken from DESeq2 SC recommendations:
    # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis
    sfType <- "poscounts"
    minmu <- 1e-06
    minReplicatesForReplace <- Inf
  }

  ## DESeq reformats interactions to have a '.' instead of ':'. Fix if need
  if (is_interaction(coef_value)) {
    coef_value <- gsub(x=coef_value, pattern=":", replacement=".", fixed=T)
  }

  ## Defaults are taken from DESeq2 SC recommendations:
  # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis
  if (TRUE) {
	  de_results <- DESeq2::DESeq(
	    dds,
	    test = "LRT",
	    fitType = de_method,
	    sfType = sfType,
	    reduced = reduced_des,
	    quiet = !verbose,
	    minReplicatesForReplace = minReplicatesForReplace,
	    useT = T, ## Shouldn't matter. Only used for Wald
	    minmu = minmu,
	    parallel = par,
	    BPPARAM = BiocParallel::MulticoreParam(n_cores)
	  )

  } else {
	# Added by Caleb, Sept 7 2022
	#   This is an alternative curvefitting method in case 
	#   the dispersion estimates are < 2 orders of magnitude from
	#   the minimum. (i.e. the target phenotype has a small range)
	#   Turns out there was an error in my calculation and this isn't
	#   needed anymore, but I left it in as an option.
	dds <- estimateSizeFactors(dds)
	dds <- estimateDispersionsGeneEst(dds)
	dispersions(dds) <- mcols(dds)$dispGeneEst

  	de_results <- DESeq2::nbinomLRT(
	    dds,
	    reduced = reduced_des,
	    quiet = !verbose
	    )
  }
  # Get results
  rez <- DESeq2::results(
    de_results,
    name = coef_value,
    cooksCutoff = F,
    independentFiltering = T,
    alpha = 0.1,
    pAdjustMethod = "BH",
    format = "DataFrame",
    addMLE = F,
    minmu = minmu,
    parallel = par,
    BPPARAM = BiocParallel::MulticoreParam(n_cores)
  )

  ## Order by pval
  rez <- rez[order(rez$pvalue, decreasing = F), ]
  if (verbose) {
    cat("Done performing Likelihood Ratio Test.\n")
  }

  ## Rename DESeq results for downstream analysis
  rez <- rez %>%
    as.data.frame(.) %>%
    dplyr::rename(
      log2fc = log2FoldChange,
      std_err = lfcSE,
      fdr = padj,
      test_statistic = stat ## Use test_statistic for downstream anaylsis
    )

  # Add other columns to dataframe
  rez$gene <- rownames(rez)
  rez$gene_symbol <- feature_metadata[rez$gene, ]$gene_symbol
  rez$test_statistic_type <- "LRTStatistic"

  return(rez)
}

################################################################################
