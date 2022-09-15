#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(0)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
  optparse::make_option(c("-a", "--input_dir"),
                        type = "character",
                        default = "matrix_dir",
                        help = "Directory containing input files for MAST"
  ),

  optparse::make_option(c("-b", "--cell_label_column"),
                        type = "character",
                        default = "cluster",
                        help = "Column to use for cell label."
  ),

  optparse::make_option(c("-c", "--cell_label"),
                        type = "character",
                        default = "",
                        help = "Cell label."
  ),

  optparse::make_option(c("-d", "--variable_target"),
                        type = "character",
                        default = "condition",
                        help = "Column to test."
  ),

  optparse::make_option(c("-e", "--discrete_variables"),
                        type = "character",
                        default = "",
                        help = "Discrete covariates to include in the model."
  ),

  optparse::make_option(c("-f", "--continuous_variables"),
                        type = "character",
                        default = "",
                        help = "Continuous covariates to include in the model."
  ),

  optparse::make_option(c("-g", "--discrete_levels"),
                        type = "character",
                        default = "",
                        help = "Levels of discrete covariates to include in
                            the model. Format should be as follows:
                            var1::reference,alt1,alt2;;var2::reference,alt1.
                            For example:
                            sex::M,F;;disease::healthy,unwell,sick"
  ),

  optparse::make_option(c("-i", "--method"),
                        type = "character",
                        default = "",
                        help = "Method to use. Needs to correspond with
                        `method_script`. Three part method:
                        1) Package
                        2) single cell or pseudobulk
                        2) DE method
                        Ex for DESeq2:
                        deseq::singlecell::glmGamPoi or
                        deseq::pseudobulk::glmGamPoi"
  ),

  optparse::make_option(c("-j", "--method_script"),
                        type = "character",
                        default = "",
                        help = "Script to use to run DE."
  ),

  optparse::make_option(c("-k", "--experiment_key"),
                        type = "character",
                        default = "sanger_sample_id",
                        help = "Key to use to determine sample source of cells."
  ),

  optparse::make_option(c("-l", "--filter_threshold"),
                        type = "double",
                        default = 1,
                        help = "Filter to remove genes with fewer cp10k/counts
                        averages. See option `--filter_type`"
  ),
  optparse::make_option(c("-x", "--filter_type"),
                        type = "character",
                        default = 1,
                        help = "Either 'counts' or 'cp10k'"
  ),

  optparse::make_option(c("-m", "--pre_filter_genes"),
                        action = "store_true",
                        default = FALSE,
                        help = "Filter genes before differential expression
                        analysis."
  ),

  optparse::make_option(c("-n", "--out_file"),
                        type = "character",
                        default = "",
                        help = "Base output name."
  ),

  optparse::make_option(c("-o", "--cores_available"),
                        type = "integer",
                        default = 1,
                        help = "Number of cores to use."
  ),

  optparse::make_option(c("-p", "--formula"),
                        type = "character",
                        default = "",
                        help = "Formula to model (if none, will automatically
                            be built). Useful to set random effects such as
                            (1|experiment_id). Example:
                             ~ sex + age + (1|experiment_id). Note: for random
                             effects one must use glmer method."
  ),

  optparse::make_option(c("-q", "--include_proportion_covariates"),
                        action = "store_true",
                        default = FALSE,
                        help = "If TRUE, include proportion covarates. These
                        metadata fields are found under the
                        `proportion_covariate_keys__autogen` column and end in
                        `_proportion__autogen`."
  ),

  optparse::make_option(c("-r", "--include_cluster_identity"),
                        action = "store_true",
                        default = FALSE,
                        help = "If FALSE, do not include own cell identity as
                        a covariate. If TRUE, do include own cell identity as
                        a covariate."
  ),
  
  optparse::make_option(c("-s", "--run_ruvseq"),
                        action = "store_true",
                        default = FALSE,
                        help = "If TRUE, perform RUVseq based on empirically
                        defined control genes. If FALSE, do nothing."
  ),
  
  optparse::make_option(c("-t", "--ruvseq_n_empirical_genes"),
                        type = "double",
                        default = 0.5,
                        help = "Number of empirical genes to use for RUVseq. If
                        `value`<1, we use (# of total genes * `value`). If 
                        `value`>1, we use `value`."
  ),
  
  optparse::make_option(c("-u", "--ruvseq_min_pvalue"),
                        type = "double",
                        default = 0.25,
                        help = "Minimum pvalue threshold for RUVseq empirical
                        genes. Any gene with pvalue<`value` will be removed from
                        control set."
  ),
  
  optparse::make_option(c("-w", "--ruvseq_k_factors"),
                        type = "integer",
                        default = 2,
                        help = "Number of factors to calculate with RUVseq."
  ),

  optparse::make_option(c("-v", "--verbose"),
                        action = "store_true",
                        default = TRUE,
                        help = ""
  )
)

parser <- optparse::OptionParser(
  usage = "%prog",
  option_list = optionList,
  description = paste0(
    "Calculates differentially expressed genes using method passed by ",
    "`method_script`"
  )
)

# a hack to fix a bug in optparse that won't let you use positional args
# if you also have non-boolean optional args:
getOptionStrings <- function(parserObj) {
  optionStrings <- character()
  for (item in parserObj@options) {
    optionStrings <- append(optionStrings,
                            c(item@short_flag, item@long_flag))
  }
  optionStrings
}

optStrings <- getOptionStrings(parser)
arguments <- optparse::parse_args(parser, positional_arguments = TRUE)
################################################################################

######################## Required Packages #####################################
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RUVSeq))
# Source the method script
suppressPackageStartupMessages(source(arguments$options$method_script))
################################################################################

################################ Functions #####################################

## THIS IS A DEMO FOR THE TWO METHODS THAT NEED TO BE IMPLEMENTED IN METHOD
## SCRIPT
# DE_get_input_type <- function() {
#   return("counts")
# }
#
# DE_allow_random_effect <- function() {
#   return(FALSE)
# }
#
# DE_calculate_dge <- function(input_matrix,
#                              feature_metadata,
#                              sample_metadata,
#                              testing_var,
#                              coef_value,
#                              formula,
#                              method = "glmGamPoi",
#                              n_cores = 1,
#                              verbose=TRUE) {
#   return(data.frame())
# }


read_matrix <- function(dir) {
  matrix <- as(Matrix::readMM(sprintf("%s/matrix.mtx.gz", dir)), "matrix")
  features <- read.csv(sprintf("%s/features.tsv.gz", dir),
                       sep ="\t",
                       header=F,
                       col.names = c("gene_id", "gene_symbol", "type"),
                       row.names = 1)
  barcodes <- read.csv(sprintf("%s/barcodes.tsv.gz", dir),
                       sep ="\t",
                       header=F)$V1
  rownames(matrix) <- rownames(features)
  colnames(matrix) <- barcodes
  return(matrix)
}

get_pseudobulk_mtx <- function(
    matrix,
    metadata,
    sample_key
) {
  new_mtx <- lapply(unique(metadata[[sample_key]]), function(key) {
      barcodes <- rownames(metadata)[which(metadata[[sample_key]] == key)]
      sample_mtx <- matrix[, barcodes, drop=F]
      df <- data.frame(key = Matrix::rowSums(sample_mtx))
      names(df) <- key
      return(df)
    }) %>%
    do.call(cbind, .) %>%
    as.matrix(.)
  return(new_mtx)
}

normalize_counts <- function(
    counts_matrix,
    target= 1e6,
    log_transform=T
) {
  col_sums <- Matrix::colSums(counts_matrix)
  # Using median to normalize the counts -- it is what scanpy does
  # https://github.com/theislab/scanpy/blob/5533b644e796379fd146bf8e659fd49f92f718cd/scanpy/preprocessing/_normalization.py#L13
  size_factors <- col_sums / target
  new_mtx <- apply(counts_matrix, 1, function(x) {
      new_val <- x / size_factors
      if (log_transform) {
        new_val <- log(new_val + 1)
      }
      return(new_val)
    }) %>%
    t(.)
  return(new_mtx)
}

get_pseudobulk_metadata <- function(metadata, sample_key) {
  new_metadata <- metadata %>%
    dplyr::group_by_at(.vars = sample_key) %>%
    dplyr::group_split() %>%
    lapply(., function(x) {
      df <- data.frame(row.names = unique(x[[sample_key]]))
      for (col in colnames(x)) {
        if (length(unique(x[[col]])) == 1) {
          df[[col]] <- unique(x[[col]])
        } else {
          df[[col]] <- mean(x[[col]])
        }
      }

      # Add import metadata
      df$n_cells_pseudobulk <- nrow(x)
      return(df)
    }) %>%
    do.call(rbind, .)
  return(new_metadata)
}

pct_exprs_n <- function(matrix, threshold = 1) {
  total <- ncol(matrix)
  pcts <- apply(matrix, 1, function(x) {
    return((sum(x >= threshold, na.rm = T) / total) * 100)
  })
  return(pcts)
}

is_interaction <- function(val) {
  return(length(strsplit(val, split=":", fixed=T)[[1]]) > 1)
}

is_random_effect <- function(val) {
  # return(
  #   grepl("|", val, fixed = TRUE) &
  #     startsWith(val, "(") &
  #     endsWith(val, ")")
  # )
  return(grepl("|", val, fixed = TRUE))
}

should_remove_term <- function(term, metadata, discrete_levels, verbose=T) {
  if (term %in% colnames(metadata) & length(unique(metadata[[term]])) <= 1) {
    if (verbose) {
      print(sprintf(
        "Covariate `%s` only has one value, removing from DE list.",
        term
      ))
    }
    return(T)
  } else if (is_interaction(term)) {
    var_terms <- strsplit(term, split=":", fixed=T)[[1]]
    unique_terms <- metadata %>%
      dplyr::group_by_at(.vars = var_terms) %>%
      dplyr::count() %>%
      nrow()
    if (unique_terms <= 1) {
      if (verbose) {
        print(sprintf(
          "Covariate `%s` only has one value, removing from DE list.",
          term
        ))
      }
      return(T)
    }
  } else if (is_random_effect(term)) {
    if (!DE_allow_random_effect()) {
      print(sprintf(
        "The specified model does not support random effects. Removing `%s` from DE list.",
        term
      ))
      return(T)
    }

    var_terms <- strsplit(
      strsplit(term, split="|", fixed=T)[[1]][2],
      split=")",
      fixed=T
    )[[1]]
    if (var_terms %in% colnames(metadata) &
        length(unique(metadata[[var_terms]])) <= 1) {
      print(sprintf(
        "Covariate `%s` only has one value, removing from DE list.",
        term
      ))
      return(T)
    }
  }

  ## Check for correct values for discrete covariates
  if (term %in% names(discrete_levels)) {
    data_vals <- unique(metadata[[term]])
    reference <- discrete_levels[[term]][1]
    # Check to make sure reference is a value
    if (!(reference %in% data_vals)) {
      if (verbose) {
        print(sprintf(paste("Metadata do not contain reference value `%s`",
                            "for covariate `%s`. Removing covariate."),
                      reference, term))
      }
      return(T)
    }
    # Check to make sure data doesn't contain a value not specified by levels.
    non_specified_values <- setdiff(data_vals, discrete_levels[[term]])
    if (length(non_specified_values) > 0) {
      if (verbose) {
        print(sprintf(paste("Removing covariate `%s`. The following values",
                            "appear in the metadata, but not the specified",
                            "levels: %s"),
                      term,
                      paste(non_specified_values, collapse = ", ")))
      }
      return(T)
    }
  }
  return(F)
}

cast_covariates <- function(
    df,
    cols,
    cast_func,
    cast_func_description,
    verbose = T
) {
  if (verbose) {
    print(sprintf("Casting columns to be %s...", cast_func_description))
  }
  for (col in cols) {
    if (!(col %in% colnames(df))) {
      print(sprintf("Column `%s` not in dataframe. Skipping...", col))
      next()
    }
    df[col] <- cast_func(df[[col]])
  }
  return(df)
}

check_empty_string <- function(df, cols) {
  for (col in cols) {
    filt <- (df[[col]] == "")
    if (any(filt)) {
      print(sprintf("Empty values in %s.", col))
      stop(sprintf("Empty values in %s.", col))
    }
  }
}

mean_impute_nan_numeric <- function(df, cols, verbose = T) {
  for (col in cols) {
    # is.finite catches NA, NaN, Inf, -Inf
    filt <- !is.finite(df[[col]])
    if (any(filt)) {
      if (verbose) {
        print(sprintf("Mean imputing non finite values in %s.", col))
        warning(sprintf("Mean imputing non finite values in %s.", col))
      }
      df[[col]][filt] <- mean(df[[col]], na.rm = TRUE)
    }
  }
  return(df)
}

get_testing_data <- function(test_var, metadata, formula) {
  if (is_interaction(test_var)) {
    terms <- strsplit(test_var, split=":", fixed=T)[[1]]
    test_term <- terms[1]
    for (term in terms[-1]) {
      if (is.character(metadata[[term]])) {
        stop(sprintf(paste("%s is a discrete covariate. Interactions",
                           "with discrete covariates are not supported."),
                     term))
      }
    }
    # Get information for testing var
    df <- get_testing_data(test_term, metadata)

    # Add back interactions
    form_terms <- attr(terms(formula), "term.labels")
    form_terms <- unlist(lapply(attr(terms(formula), "term.labels"),
                                function(x) {
                                  return(strsplit(x, split=":", fixed=T)[[1]])
                                }))
    terms <- terms[order(match(terms, form_terms))]
    format_str <- paste(replace(terms, which(terms == test_term), "%s"),
                        collapse = ":")
    df$alt_var <- sprintf(format_str, df$alt_var)
    if (all(!is.na(df$ref_var))) {
      df$ref_var <- sprintf(format_str, df$ref_var)
    }
    return(df)
  } else {
    if (is.numeric(metadata[[test_var]])) {
      df <- data.frame(
        "alt_var" = test_var,
        "ref_var" = NA,
        "barcodes" = paste(rownames(metadata), collapse = "$$")
      )
      return(df)
    } else {
      reference_val <- levels(metadata[[test_var]])[1]
      alt_vals <- setdiff(unique(metadata[[test_var]]), reference_val)
      barcodes <- lapply(alt_vals, function(x) {
        barcode <- rownames(metadata)[metadata[[test_var]] %in%
                                        c(x, reference_val)]
        return(paste(barcode, collapse = "$$"))
      })
      df <- data.frame("alt_var" = paste(test_var, alt_vals, sep=""),
                       "ref_var" = paste(test_var, reference_val, sep=""),
                       "barcodes" = unlist(barcodes))
      return(df)
    }
  }
}

reconstruct_formula <- function(original_terms, terms_to_remove) {
  form_terms <- setdiff(original_terms, terms_to_remove)

  # Need to add parentheses back to random effects
  form_terms <- sapply(form_terms, function(term) {
    if (is_random_effect(term)) {
      return(sprintf("(%s)", term))
    }
    return(term)
  })
  form <- formula(sprintf("~ %s", paste(form_terms, collapse = " + ")))
  return(form)
}

run_ruvseq <- function(
  nominal_results,
  n_empirical_genes,
  min_pvalue,
  k,
  counts_mtx
) {
  # If `n_empirical_genes` < 1, assume proportional
  n_emp_genes <- if (n_empirical_genes < 1) {
      round(nrow(nominal_results) * n_empirical_genes) 
    } else {
      n_empirical_genes
    }
  
  # Remove values not converged or <= min_pvalue
  emp_genes_df <- nominal_results[!is.na(nominal_results$pvalue), ]
  if (is.null(min_pvalue)) {
    emp_genes_df <- emp_genes_df[emp_genes_df$pvalue > min_pvalue, ]  
  }
  
  # Now order by p-value and get control genes
  emp_genes <- emp_genes_df$gene[order(emp_genes_df$pvalue, decreasing = T)]
  emp_genes <- emp_genes[1:min(n_emp_genes, length(emp_genes))]
  print(paste('Number of empirical genes:', length(emp_genes)))
  
  # Now get RUVseq factors
  ruv_covs <- RUVSeq::RUVg(counts_mtx, emp_genes, k=k, isLog=F)
  return(as.data.frame(ruv_covs$W))
}

get_empty_df <- function() {
  return(data.frame(
    gene=character(),
    gene_symbol=character(),
    log2fc=character(),
    pvalue=character(),
    test_statistic=character(),
    de_method=character(),
    test_statistic_type=character()
  ))
}

################################################################################

######################## Read Data & Manipulate ################################
verbose <- arguments$options$verbose
output_file_base <- arguments$options$out_file
mean_filter <- arguments$options$filter_threshold
filter_type <- arguments$options$filter_type
formula_str <- arguments$options$formula

# Read all data in
mtx_file_dir <- arguments$options$input_dir
if (verbose) {
  print("Reading in the data...")
}

# Read all data matrices in
counts_matrix <- read_matrix(sprintf("%s/counts", mtx_file_dir))
logcp10_matrix <- read_matrix(sprintf("%s/log1p_cp10k", mtx_file_dir))
cp10_matrix <- read_matrix(sprintf("%s/cp10k", mtx_file_dir))

if (arguments$options$pre_filter_genes) {
  print(sprintf(paste("Filtering genes with fewer than %s counts per 10k",
                      "BEFORE finding differentially expressed genes."),
                mean_filter))
  n_genes_pre_filter <- nrow(cp10_matrix)

  # Filter and update all matrices
  cp10_matrix <- cp10_matrix[Matrix::rowMeans(cp10_matrix) > mean_filter, ]
  counts_matrix <- counts_matrix[rownames(cp10_matrix), ]
  logcp10_matrix <- logcp10_matrix[rownames(cp10_matrix), ]

  print(
    sprintf("Retaining %s / %s genes.", nrow(cp10_matrix), n_genes_pre_filter)
  )
}

# May also need feature metadata
feature_metadata <- read.csv(
  sprintf("%s/counts/features.tsv.gz", mtx_file_dir),
  sep ="\t",
  header=F,
  col.names = c("gene_id", "gene_symbol", "type"),
  row.names = 1
)
feature_metadata <- feature_metadata[rownames(counts_matrix), , drop=F]

## Assure metadata is in correct order
metadata <- read.csv(
  sprintf("%s/cell_metadata.tsv.gz", mtx_file_dir),
  sep ="\t",
  header=T,
  row.names=1
)
metadata <- metadata[colnames(counts_matrix), , drop=F]

# Include proportion covariates if needed
include_cell_props <- arguments$options$include_proportion_covariates
if (include_cell_props) {

  # First make sure only 1 value across cell type. If heirarchical, different
  # than cell label, but still should be only one value.
  cell_type_prop <- unique(metadata$proportion_covariate_keys__autogen)
  if (length(cell_type_prop) > 1) {
    stop(sprintf(
      "ERROR: More than one value for cell type proportions. Keys: %s",
      paste(cell_type_prop, collapse=", ")
    ))
  }

  prop_cols <- colnames(metadata)[
    endsWith(colnames(metadata), "_proportion__autogen")
  ]

  if (!arguments$options$include_cluster_identity) {
    prop_cols <- setdiff(prop_cols, cell_type_prop)
  }

  # Cast to be numeric
  metadata <- cast_covariates(metadata, prop_cols, as.numeric, "numeric")

  # Add to formula
  formula_str <- paste(c(formula_str, prop_cols), collapse=" + ")
}

# Variable casting
discrete_covs <- strsplit(
  x=arguments$options$discrete_variables,
  split=",",
  fixed=TRUE
)[[1]]

metadata <- cast_covariates(
  metadata,
  discrete_covs,
  as.character,
  "characters",
  verbose
)
check_empty_string(metadata, discrete_covs)

continuous_covs <- strsplit(
  x=arguments$options$continuous_variables,
  split=",",
  fixed=TRUE
)[[1]]
metadata <- cast_covariates(
  metadata,
  continuous_covs,
  as.numeric,
  "numeric",
  verbose
)

# Get discrete variable levels
discrete_levels_str <- arguments$options$discrete_levels
discrete_levels <- sapply(
  strsplit(x=discrete_levels_str, split=";;", fixed=T)[[1]],
  function(x) {
    items <- strsplit(x, split="::", fixed=T)[[1]]
    ret <- strsplit(items[2], split=",",fixed=T)
    names(ret) <- items[1]
    return(ret)
  },
  USE.NAMES=F
)


# Check for pseudobulk
## If pseudo:
## For expression: average across samples
## For cont. metadata: average across samples
## For discrete metadata: just get unique data point
pseudobulk <- strsplit(arguments$options$method, split="::", fixed=T)[[1]][2]
if (pseudobulk == "pseudobulk") {
  if (verbose) {
    cat("Creating a pseudobulk dataset...\n")
  }
  # Merge counts matrices
  counts_matrix <- get_pseudobulk_mtx(
    counts_matrix,
    metadata,
    arguments$options$experiment_key
  )

  # Recalc logtpm and cp10k using the per sample matrix
  cp10_matrix <- normalize_counts(counts_matrix, target=1e4, log_transform=F)
  logcp10_matrix <- normalize_counts(counts_matrix, target=1e4, log_transform=T)

  # Merge metadata
  metadata <- get_pseudobulk_metadata(
    metadata,
    arguments$options$experiment_key
  )
  metadata <- metadata[colnames(counts_matrix), ]
  # Mean impute continuous_covs after the metadata merged
  metadata <- mean_impute_nan_numeric(metadata, continuous_covs)
  if (verbose) {
    cat("Done creating a pseudobulk dataset.\n")
  }
} else {
  # Mean impute continuous_covs per cell
  metadata <- mean_impute_nan_numeric(metadata, continuous_covs)
}


# Get formula
## For each term:
## 1) Check to make sure it exists in metadata
## 2) Check to make sure all levels exist if levels specified
##    If discrete and no specification, pick first one
## 3) Check for a single value -- if only one value exists, remove from vars
##      MAST throws an error if not
formula <- formula(formula_str)
testing_var <- arguments$options$variable_target
formula_variables_passed <- attr(terms(formula), "term.labels")
vars_removed <- c()
# First evaluate any terms
for (var in formula_variables_passed) {
  if (grepl("I(", var, fixed = TRUE)) {
    cat(paste0("Evaluating:\t", var, "\n"))
    metadata[[var]] <- eval(parse(text = var), metadata)
  }
}

for (var in formula_variables_passed) {
  # Sanitize formula
  remove_term <- should_remove_term(var, metadata, discrete_levels)
  if (remove_term) {
    vars_removed <- c(vars_removed, var)
    next
  }

  var_terms <- c(var)
  if (is_interaction(var)) {
    var_terms <- strsplit(var, split=":", fixed=T)[[1]]
  } else if (is_random_effect(var)) {
    var_terms <- strsplit(
      strsplit(var, split="|", fixed=T)[[1]][2],
      split=")",
      fixed=T
    )[[1]]
  }

  ## Deal with specific types
  for (var_term in var_terms) {
    if (var_term %in% names(discrete_levels)) {
      metadata[[var_term]] <- factor(
        metadata[[var_term]],
        levels = discrete_levels[[var_term]],
        labels = discrete_levels[[var_term]]
      )
    } else if (is.character(metadata[[var_term]])) {
      factors <- factor(metadata[[var_term]])
      metadata[var_term] <- factors
      if (verbose) {
        print(sprintf(
          paste("Selected `%s` as the reference value for the `%s`",
            "column in fitting the model."),
          levels(factors)[1],
          var_term
        ))
      }
    } else if (is.numeric(metadata[[var_term]]) &
        !(paste0(var_term, "_unscaled") %in% colnames(metadata))) {
      print(sprintf("Scaling the covariate %s...", var_term))
      # Save the unscaled variable as _unscaled
      metadata[[paste0(var_term, "_unscaled")]] <- metadata[[var_term]]
      metadata[[var_term]] <- scale(
        metadata[[paste0(var_term, "_unscaled")]],
        center=T,
        scale=T
      )
    }
  }
}

if (testing_var %in% vars_removed) {
  stop("Testing variable was removed from the formula. See output to debug.")
}
formula <- reconstruct_formula(formula_variables_passed, vars_removed)

if (verbose) {
  print(sprintf("The final formula: `%s`", deparse(formula)))
}

## Now get every alt hypothesis and corresponding barcodes
test_data <- get_testing_data(testing_var, metadata, formula)

if (nrow(test_data) == 0) {
  ## if test data is empty--issue with the rank. Returning a null dataframe so
  ## the pipeline doesn't fail completely.
  de_results <- get_empty_df()
} else {
  ## Run LRT and merge
  ldfs <- apply(test_data, 1, function(i) {

    de_method <- paste(
      strsplit(x=arguments$options$method, split="::", fixed=T)[[1]][-1],
      collapse = "::"
    )

    if (verbose) {
      cat(sprintf(
        paste(
          "Calculating DGE using the following parameters:\n",
          "Script: %s\n",
          "Method: %s\n",
          "Variable to test: %s\n"
        ),
        arguments$options$method_script,
        de_method,
        i["alt_var"]
      ))
    }

    ## IMPORTANT:
    ## Every method scripts needs: DE_get_input_type() and DE_calculate_dge()
    ## The inputs need to be the same for each script.
    input_type <- DE_get_input_type()
    if (input_type == "counts") {
      input <- counts_matrix
    } else if (input_type == "logcp10") {
      input <- logcp10_matrix
    } else if (input_type == "cp10") {
      input <- cp10_matrix
    }

    ## This method should return a DF with the following columns:
    # log2fc, pvalue, test_statistic_type, test_statistic, de_method
    # gene (Ensembl IDs), gene_symbol
    rez <- try(
      DE_calculate_dge(
        input_matrix = input,
        feature_metadata = feature_metadata,
        sample_metadata = metadata,
        testing_var = testing_var,
        coef_value = i["alt_var"],
        formula = formula,
        n_cores = arguments$options$cores_available,
        method = de_method
      )
    )

    if (class(rez) == "try-error") {
      print(paste(
        'The function call `DE_calculate_dge` returned the following error:',
        geterrmessage()
      ))
      print('Creating an empty dataframe to continue the pipeline...')
      return(get_empty_df())
    }
    
    # If RUVseq, perform
    if (arguments$options$run_ruvseq) {
      ruvseq_n_emp <- arguments$options$ruvseq_n_empirical_genes
      ruvseq_min_pvalue <- arguments$options$ruvseq_min_pvalue
      ruvseq_k <- arguments$options$ruvseq_k_factors
      
      print('Performing RUVSeq with the following settings:')
      print(paste('# empirical genes:', ruvseq_n_emp))
      print(paste('Minimum p-value:', ruvseq_min_pvalue))
      print(paste('# factors:', ruvseq_k))
      
      ruvseq_factors <- run_ruvseq(
        rez,
        ruvseq_n_emp,
        ruvseq_min_pvalue,
        ruvseq_k,
        counts_matrix
      )
      
      ruvseq_factors__formatted <- ruvseq_factors
      rownames(ruvseq_factors__formatted) <- rownames(metadata)
      ruv_file <- gzfile(
        sprintf("%s_ruvseq_factors.tsv.gz", output_file_base),
        "w",
        compression = 9
      )
      write.table(
        x=ruvseq_factors__formatted,
        file=ruv_file,
        sep="\t",
        col.names=T,
        row.names=T,
        quote=F
      )
      close(ruv_file)
      
      # Now update data and re-run
      metadata <- cbind(metadata, ruvseq_factors)
      formula <- update(
        formula,
        sprintf('~ . + %s', paste(colnames(ruvseq_factors), collapse = '+'))
      )
      print(paste('Performing DGE with the updated formula:', deparse(formula)))
      
      rez <- try(
        DE_calculate_dge(
          input_matrix = input,
          feature_metadata = feature_metadata,
          sample_metadata = metadata,
          testing_var = testing_var,
          coef_value = i["alt_var"],
          formula = formula,
          n_cores = arguments$options$cores_available,
          method = de_method
        )
      )
      
      # Since passed once, shouldn't hit but keep in case
      if (class(rez) == "try-error") {
        print(paste(
          'The function call `DE_calculate_dge` returned the following error:',
          geterrmessage()
        ))
        print('Creating an empty dataframe to continue the pipeline...')
        return(get_empty_df())
      }
    } else {
      ruvseq_n_emp <- NULL
      ruvseq_min_pvalue <- NULL
      ruvseq_k <- NULL
    }

    ## Get only the columns we want
    cols_retain <- c("gene", "gene_symbol", "log2fc", "std_err", "pvalue",
                     "test_statistic", "test_statistic_type")
    rez <- rez[, cols_retain] # NOTE: must be data.frame

    ## Now assign necessary columns back to dataframe
    rez$target_variable <- testing_var
    rez$reference_value <- i["ref_var"]
    rez$coef_value <- i["alt_var"]
    rez$formula_passed <- gsub(
      arguments$options$formula,
      pattern=" ",
      replacement="",
      fixed=T
    )
    rez$formula <- gsub(
      paste(deparse(formula), collapse=""),
      pattern=" ",
      replacement="",
      fixed=T
    )
    rez$cell_label_column <- arguments$options$cell_label_column
    rez$cell_label <- arguments$options$cell_label
    rez$de_method <- arguments$options$method
    rez$pre_filtered <- arguments$options$pre_filter_genes
    rez$include_cell_proportions <- include_cell_props
    rez$ruvseq_n_empirical_genes <- ruvseq_n_emp
    rez$ruvseq_min_pvalue <- ruvseq_min_pvalue
    rez$ruvseq_k_factors <- ruvseq_k

    ## Get count averages
    barcodes <- strsplit(i["barcodes"], split="$$", fixed=T)[[1]]
    rez$mean_counts <- Matrix::rowMeans(counts_matrix[rez$gene, barcodes])
    rez$mean_cp10k <- Matrix::rowMeans(cp10_matrix[rez$gene, barcodes])
    rez[[paste0("pct_", mean_filter, "cp10k")]] <- pct_exprs_n(
        cp10_matrix[rez$gene, barcodes], mean_filter
    )
    ## Add number of cells
    rez$n_cells <- ncol(input)

    return(rez)
  })
  de_results <- do.call(rbind, ldfs)
}

if (nrow(de_results) > 0) {

  # Save unfiltered result
  # p.adjust ignores NA pvalues unless `n` argument is specified
  de_results <- de_results %>%
    dplyr::group_by(coef_value) %>%
    dplyr::mutate(
      qvalue_bh_percelltype = p.adjust(pvalue, method = "BH")
    )

  if (verbose) {
    print("Writing DE results...")
  }
  gz_file <- gzfile(
    sprintf("%s_unfiltered-de_results.tsv.gz", output_file_base),
    "w",
    compression = 9
  )
  write.table(
    x=de_results,
    file=gz_file,
    sep="\t",
    col.names=T,
    row.names=F,
    quote=F
  )
  close(gz_file)

  # Filter
  if (verbose) {
    cat(sprintf("Filtering out genes with mean %s < %s...\n",
		filter_type,
                mean_filter))
  }
  n_genes_before <- nrow(de_results)
  print(filter_type)
  # Removed because it wasn't working right
  if (filter_type == "cp10k") {
  	de_results <- de_results[which(de_results$mean_cp10k > mean_filter), ]
  } else if (filter_type == "counts") {
  	de_results <- de_results[which(de_results$mean_counts > mean_filter), ]
  } else { 
	stop(sprintf("Error: invalid argument for `filter_type`: expected 'cp10k' or 'counts', got %s", filter_type))
  }
  if (verbose) {
    cat(sprintf("Done. Filtered %s genes.\n",
                n_genes_before - nrow(de_results)))
  }

  # p.adjust ignores NA pvalues unless `n` argument is specified
  de_results <- de_results %>%
    dplyr::group_by(coef_value) %>%
    dplyr::mutate(
      qvalue_bh_percelltype = p.adjust(pvalue, method = "BH")
    )

  # Save filtered
  gz_file <- gzfile(
    sprintf("%s_filtered-de_results.tsv.gz", output_file_base),
    "w",
    compression = 9
  )
  write.table(
    x=de_results,
    file=gz_file,
    sep="\t",
    col.names=T,
    row.names=F,
    quote=F
  )
  close(gz_file)

}

if (verbose) {
  print("Done.")
}


################################################################################
