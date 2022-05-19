#!/usr/bin/env Rscript

SCRIPT_NAME <- "summarize_fgsea.R"

# disable strings as factors, but re-enable upon exit
old <- options(stringsAsFactors = FALSE)
on.exit(options(old), add = TRUE)

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tm))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(wesanderson))
suppressPackageStartupMessages(library(simplifyEnrichment))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(igraph))

# for parallel processing
# suppressPackageStartupMessages(library(doSNOW))
# suppressPackageStartupMessages(library(doMC))
# suppressPackageStartupMessages(library(doParallel))
# suppressPackageStartupMessages(library(parallel))
# suppressPackageStartupMessages(library(foreach))
# suppressPackageStartupMessages(library(iterators))

# for profiling
# suppressPackageStartupMessages(library(pryr))
# suppressPackageStartupMessages(library(profvis))


#' reduceSimMatrix
#' Reduce a set of based on their similarity and scores.
#' Code derived from https://github.com/ssayols/rrvgo/blob/master/R/rrvgo.R
#'
#' @details
#' Currently, rrvgo uses the similarity between pairs of terms to compute a
#' distance matrix, defined as (1-simMatrix). The terms are then hierarchically
#' clustered using complete linkage, and the tree is cut at the desired threshold,
#' picking the term with the highest score as the representative of each group.
#'
#' Therefore, higher thresholds lead to fewer groups, and the threshold should be
#' read as the expected similarity of terms within a group (though this is not
#' entirely correct, and you'll see similarities below this threshold being put in
#' the same group).
#'
#' @param simMatrix a (square) similarity matrix
#' @param scores *named* vector with scores (weights) assigned to each term.
#'   Higher is better. Note: if you have
#'   p-values as scores, consider -1*log-transforming them (`-log(p)`)
#' @param cluster_vector vector where values = cluster and names = term_ids.
#'  If NULL then will estimate using threshold.
#' @param threshold similarity threshold (0-1). Some guidance:
#'   Large (allowed similarity=0.9), Medium (0.7), Small (0.5), Tiny (0.4)
#'   Defaults to Medium (0.7)
#'
#' @return a data.frame with all terms and it's "reducer"
#'   (NA if the term was not reduced)
#'
#' @importFrom stats cutree hclust
#' @export
reduceSimMatrix <- function(
    simMatrix,
    scores,
    cluster_vector=NULL,
    threshold=0.7
) {

  # check function arguments
  if(!is.null(scores) && !all(rownames(simMatrix) %in% names(scores))) {
    stop("Scores vector does not contain all terms in the similarity matrix")
  }

  scores <- scores[match(rownames(simMatrix), names(scores))]

  # reorder the similarity matrix as in the scores, just in case they don't
  # come in the same order
  orows <- match(rownames(simMatrix), names(scores))
  ocols <- match(colnames(simMatrix), names(scores))
  simMatrix <- simMatrix[orows, ocols]

  # sort matrix based on the score
  o <- rev(order(scores, na.last=FALSE))
  simMatrix <- simMatrix[o, o]

  # cluster terms and cut the tree at the desired threshold.
  # Then find the term with the highest score as the representative of each
  # cluster
  if (is.null(cluster_vector)) {
      cluster_vector <- cutree(hclust(as.dist(1-simMatrix)), h = threshold)
  }
  # Do not use binary cut as a clustering method:
  # https://doi.org/10.1101/2020.10.27.312116
  # We compared binary cut clustering on the similarity matrices from different
  # similarity measurements and  we found the semantic similarity worked well
  # with binary cut while the similarity matrices based on gene overlap
  # showed less consistent patterns and they were not recommended to work
  # with binary cut.
  # if (cluster_simplifyEnrichment) {
  #     cluster <- simplifyEnrichment::simplifyEnrichment(
  #         simMatrix,
  #         method = "binary_cut",
  #         plot = FALSE
  #     )
  # }
  clusterRep <- tapply(
      rownames(simMatrix), cluster_vector, function(x) x[which.max(scores[x])]
  )

  # return
  df_results <- data.frame(
    term_id=rownames(simMatrix),
    cluster=cluster_vector,
    parent_term_id=clusterRep[cluster_vector],
    parentSimScore=unlist(Map(
        seq_len(nrow(simMatrix)),
        clusterRep[cluster_vector],
        f=function(i, j) simMatrix[i, j]
    )),
    score=scores[match(rownames(simMatrix), names(scores))]
  )
  return(df_results)
}


#' Combines dataframes in a list
#'
#' @importFrom data.table rbindlist
# rbindlist_df <- function(...) {
#     return(data.table::rbindlist(list(...), use.names = TRUE, fill = TRUE))
# }


#' Combines dataframes in a list
#'
#' @importFrom data.table rbindlist
# rbindlist_dflist <- function(...) {
#     df_coloc <- data.table::rbindlist(
#         lapply(list(...), FUN = function(x) {return(x[["df_coloc"]])}),
#         use.names = TRUE,
#         fill = TRUE
#     )
#     df_harmonized_ss <- data.table::rbindlist(
#         lapply(
#             list(...),
#             FUN = function(x) {return(x[["df_harmonized_ss"]])}
#         ),
#         use.names = TRUE,
#         fill = TRUE
#     )
#
#     return(list(
#         "df_coloc" = df_coloc,
#         "df_harmonized_ss" = df_harmonized_ss
#     ))
# }


#' Command line interface wrapper
#'
#' @importFrom optparse make_option
#' @importFrom optparse OptionParser
#' @importFrom optparse parse_args
command_line_interface <- function() {
    optionList <- list(
        optparse::make_option(c("--gsea_results"),
            type = "character",
            help = paste0(
                "Gene set enrichment results file.",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--gsets_gene_matrix"),
            type = "character",
            help = "TSV file containing gene sets with genes."
        ),

        optparse::make_option(c("--term_distance_metric"),
            type = "character",
            default = "kappa",
            help = paste0(
                "Metric to caclulate distance between terms.",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--clustering_method"),
            type = "character",
            default = "louvain",
            help = paste0(
                "Method to cluster terms.",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--output_file_tag"),
            type = "character",
            default = "fgsea_reduced",
            help = paste0(
                "Tag for output files.",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--verbose"),
            type = "logical",
            action = "store_true",
            default = FALSE,
            help = paste0(
                "Verbose mode (write extra infro to std.err).",
                " [default %default]"
            )
        )
    )

    parser <- optparse::OptionParser(
        usage = "%prog",
        option_list = optionList,
        description = paste0(
            "Runs colocalization test for tabix indexed files."
        )
    )

    # a hack to fix a bug in optparse that won"t let you use positional args
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

    # read in the parameters
    param <- list()
    for (i in names(arguments$options)) {
        param[[i]] <- arguments$options[[i]]
    }

    return(run_analysis(param))
}


#' Run analysis function
#'
#' @importFrom data.table fread
#' @importFrom parallel detectCores
#' @importFrom parallel makeForkCluster
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @importFrom doMC registerDoMC
#' @importFrom iterators isplit
#' @importFrom parallel stopCluster
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom coloc finemap.abf
#' @importFrom plyr rename
run_analysis <- function(param) {

    # set the parallel backend: fork, snow, domc
    param[["parallel_type"]] <- "domc"

    # read in the gsea results file
    f <- param[["gsea_results"]]
    df_gsea <- data.frame(data.table::fread(
        cmd = paste("gunzip -c", f),
        sep = "\t",
        header = T,
        stringsAsFactors = F
    ))
    # Rename the anno_id column
    df_gsea$term_id <- df_gsea$annot_id
    df_gsea$annot_id <- NULL

    # Read in the gsea annotations
    f <- param[["gsets_gene_matrix"]]
    df_gsea_genes <- data.frame(
        data.table::fread(
            cmd = paste("gunzip -c", f),
            sep = "\t",
            header = T,
            stringsAsFactors = F
        ),
        row.names=1
    )
    # df_gsea_genes <- read.csv(
    #   f,
    #   sep='\t',
    #   header=T,
    #   row.names=1
    # )

    # for parallel processing with foreach
    # more info here: https://ljdursi.github.io/beyond-single-core-R
    time_init <- Sys.time()
    # if (param[["threads"]] > parallel::detectCores()) {
    #     warning(paste0(
    #         "User requested more threads than available. Setting number of",
    #         " threads to doParallel::detectCores() - 1, that is",
    #         " ", parallel::detectCores() - 1, "."
    #     ))
    #     param[["threads"]] <- parallel::detectCores() - 1
    # }
    # if (param[["parallel_type"]] == "fork") {
    #     sink("/dev/null") # use sink to grab annoying output when outfile = ""
    #         # (v1) Forking = copy the R session in its current state. Fast
    #         # because only copies objects if they are modified; however it takes
    #         # much more memory as files are copied
    #         cl <- parallel::makeForkCluster(param[["threads"]], outfile = "")
    #     sink()
    #     doParallel::registerDoParallel(cl)
    # } else if (param[["parallel_type"]] == "snow") {
    #     # since snow could be on any machine, set outfile to /dev/null rather
    #     # than the parent process (unlink fork)
    #     # (v2) Spawn new subprocess: benefits can be done on any machine,
    #     # new process will not have unneeded data
    #     # if (param[["verbose"]]) {
    #     #     message(paste0(
    #     #         "[", SCRIPT_NAME, "]:\t",
    #     #         " Running parallel with snow. Any std.out or std.err (apart",
    #     #         " from a failed exit) will be lost."
    #     #     ))
    #     # }
    #     # cl <- parallel::makeCluster(param[["threads"]], outfile = "/dev/null")
    #     # parallel::clusterExport(cl, c(
    #     #     "SCRIPT_NAME", "param", "df_gsea", "df_gsea_genes",
    #     #     "df_ss_trait_type", "df_ss_sdY", "rbindlist_dflist"
    #     # ))
    #     # doSNOW::registerDoSNOW(cl)
    #     stop(paste0("[", SCRIPT_NAME, "]:\tERROR invalid parallel_type."))
    # } else if (param[["parallel_type"]] == "domc") {
    #     # (v3) Shared memory process
    #     # doesn't use persistent workers like the snow package or ths
    #     # now-derived functions in the parallel package
    #     doMC::registerDoMC(param[["threads"]])
    # } else {
    #     stop(paste0("[", SCRIPT_NAME, "]:\tERROR invalid parallel_type."))
    # }

    # print out parameters
    if (param[["verbose"]]) {
        message(paste0("Parameters [", SCRIPT_NAME, "]:"))
        for (i in names(param)) {
            message("\t", i, " = ", param[[i]])
        }
    }

    # This is just a short hand way to print the model in pdfs
    MODEL_ID_NUMBER <- 1

    # Run the analysis in three ways:
    # (a) all annotated terms. EDIT: no this makes no sense
    # (b) all nominally associated terms (pvalue < 0.05)
    # (c) all statistically signficant terms (qvalue_bh < 0.05)
    #
    # For each iteration:
    # (a) generate clustering of enrichment terms, PCs, and
    #     summary term per cluster
    # (b) plot PCs
    # (c) plot heatmap with wordcloud of clusters
    term_summary_types <- c(
        #"all",
        "pvalue_less_than_0pt05",
        "qvalue_bh_less_than_0pt05"
    )
    # Make the table of iterations
    # iterate over each df_key in df_gsea_i
    # df_iter_table <- expand.grid(
    #     df_gsea_iter_var = unique(df_gsea$df_key),
    #     term_summary_type = term_summary_types,
    #     stringsAsFactors = FALSE
    # )

    # Iterate over each index item
    #
    # note that if using foreach and combining the results into a dataframe,
    # errors will not stop execution, but will not necessarily be returned to
    # the user. Will show up as:
    # Error in { : task 1 failed - "replacement has 1 row, data has 0"
    #
    # df_list_results <- foreach::foreach(
    #     iter_info = iterators::iter(term_summary_types),
    #     #iter_info = iterators::iter(df_iter_table, by = "row"), # iter row
    #     #.combine = rbindlist_dflist,
    #     .inorder = FALSE,
    #     .multicombine = FALSE,
    #     .errorhandling = "stop"
    # ) %do% {
    for (term_summary_type_i in term_summary_types) {
        if (param[["verbose"]]) {
            message(paste0("(iter_i) ", term_summary_type_i, "\n"))
        }

        # Subset the results data accoring to term_summary_type
        df_gsea_i <- df_gsea
        if (term_summary_type_i == "pvalue_less_than_0pt05") {
            df_gsea_i <- subset(df_gsea_i, pvalue < 0.05)
        } else if (term_summary_type_i == "qvalue_bh_less_than_0pt05") {
            df_gsea_i <- subset(df_gsea_i, qvalue_bh < 0.05)
        }
        df_gsea_i$term_summary_type <- term_summary_type_i
        if (nrow(df_gsea_i) <= 2) {
            message(paste0(
                "Skipping because too few terms ",
                term_summary_type_i,
                "\n"
            ))
            next
        }

        # Init the pdf to dump figures
        pdf(
            file = paste0(
                param[["output_file_tag"]],
                "-",
                term_summary_type_i,
                ".pdf",
                sep = ""
            ),
            height = 8,
            width = 8
        )

        # Init the final dataframe with reduced dims
        df_results_list <- list()

        # Iterate over each gsea result
        # foreach::foreach(
        #     df_gsea_ij = iterators::isplit(df_gsea_i, df_gsea_i$df_key),
        #     .inorder = FALSE,
        #     .multicombine = FALSE,
        #     .errorhandling = "stop"
        # ) %do% {
        for (j in unique(df_gsea_i$df_key)) {
            #message(paste0("(iter_j) ", j, "\n"))

            # Subset to the data we want
            df_gsea_ij <- subset(df_gsea_i, df_key == j)
            if (nrow(df_gsea_ij) <= 2) {
                message(paste0(
                    "Skipping because too few terms ",
                    j,
                    "\n"
                ))
                next
            }

            # Get the model label
            df_key_id <- paste0(
                "df_key_",
                sprintf("%03d", as.numeric(MODEL_ID_NUMBER))
            )
            MODEL_ID_NUMBER <- MODEL_ID_NUMBER + 1
            df_gsea_ij$df_key_id <- df_key_id

            # Subset df_gsea_genes to just the terms in df_gsea and vice versa
            common_terms <- intersect(
                colnames(df_gsea_genes),
                df_gsea_ij$term_id
            )
            df_gsea_genes_ij <- df_gsea_genes[common_terms]
            # if not all df_gsea in common terms, raise warning
            missing_terms <- df_gsea_ij$term_id[
                !df_gsea_ij$term_id %in% common_terms
            ]
            if (length(missing_terms) > 0) {
                warning("Missing terms:", missing_terms)
                df_gsea_ij <- subset(df_gsea_ij, term_id %in% common_terms)
            }

            # Format gene sets data to be a list with term_id and gene
            # lists as value
            annot_data <- df_gsea_genes_ij
            annot_data <- sapply(colnames(annot_data), function(x) {
              return(list(
                rownames(annot_data)[which(annot_data[[x]] == 1)]
              ))
            })

            # Calculate the similarity between terms based on overlap of gene lists
            # This is OK and only option for KEGG and REACTOME.
            # For GO terms, better to use semantic distance between terms which
            # takes into account the embeded structure
            #
            # From https://doi.org/10.1101/2020.10.27.312116:
            # The similarity between two gene sets is based on gene overlap and is
            # calculated as Jaccard coefficient, Dice coefficient or overlap
            # coefficient. Kappa coefficient is also suggested as more robust,
            # which is used in the widely used DAVID tool.
            # Unless GO terms, one should use kappa
            mtx_similarity <- simplifyEnrichment::term_similarity(
                annot_data,
                method = param[["term_distance_metric"]]
            )
            # plt <- simplifyEnrichment::select_cutoff(
            #     mtx_similarity
            # )
            # print(plt)
            # NOTE: Notes on clustering methods
            # * binary_cut really good for GO term semantic similarity
            # * louvain is ok, but tends to produce bigger clusters
            # * mclust is also ok
            term_clusters <- simplifyEnrichment::simplifyEnrichment(
                    mtx_similarity,
                    method = param[["clustering_method"]],
                    plot = FALSE
                ) %>%
                dplyr::arrange(cluster) %>%
                as.data.frame()
            # tmp_cluster_obj <- M3C::M3C(
            #     as.matrix(as.dist(1 - mtx_similarity)),
            #     method = 2
            # )
            # term_clusters <- data.frame(
            #     id = colnames(mtx_similarity),
            #     cluster = tmp_cluster_obj$assignments
            # )
            rownames(term_clusters) <- term_clusters$id
            term_clusters$cluster <- paste0(
                "cluster_",
                sprintf("%02d", as.numeric(term_clusters$cluster))
            )
            # Make a list of terms in each cluster for the word cloud
            term_clusters_wc <- term_clusters %>%
                dplyr::group_by(cluster) %>%
                dplyr::summarize(
                    term = stringr::str_squish(
                        tm::removeWords(
                            stringr::str_replace_all(
                                tolower(paste(id, sep = " ", collapse = " ")),
                                "_",
                                " "
                            ),
                            c(tm::stopwords("english"), "kegg", "reactome")
                        )
                    )
                ) %>%
                dplyr::ungroup() %>%
                dplyr::arrange(cluster) %>%
                as.data.frame()
            term_clusters_wc_list <- as.list(term_clusters_wc$term)
            names(term_clusters_wc_list) <- term_clusters_wc$cluster
            # Make a plot of each cluster as heatmap
            # NOTE: for simplifyEnrichment::anno_word_cloud function you may need
            # github version of simplifyEnrichment
            plt <- ComplexHeatmap::Heatmap(
                mtx_similarity[term_clusters$id, term_clusters$id],
                row_split = term_clusters$cluster,
                show_row_names = FALSE,
                show_column_names = FALSE,
                row_title = NULL,
                row_title_rot = 0,
                right_annotation = ComplexHeatmap::rowAnnotation(
                    wc = simplifyEnrichment::anno_word_cloud(
                        term_clusters$cluster,
                        term_clusters_wc_list
                    )
                ),
                col = wesanderson::wes_palette(
                    "Zissou1",
                    25,
                    type = c("continuous")
                ),
                heatmap_legend_param = list(
                    title = "Similarity"
                ),
                column_title = paste0(
                    unique(df_gsea_ij$cell_label),
                    " ",
                    df_key_id
                )
            )
            print(plt)

            # Make a plot of the PCs
            # Reduce a set of terms based on their similarity and encrichment pvalues.
            # Also select a term as representative of all other terms in a cluster.
            #
            # Prep scores to obtain term summary for each term
            enrichment_scores <- setNames(
                -log10(df_gsea_ij$pvalue),
                df_gsea_ij$term_id
            )
            # format clusters
            term_clusters_vector <- term_clusters$cluster
            names(term_clusters_vector) <- term_clusters$id
            reducedTerms <- reduceSimMatrix(
                mtx_similarity,
                enrichment_scores,
                term_clusters_vector,
                threshold = 0.75
            )
            # Code from https://github.com/ssayols/rrvgo/blob/master/R/plotlib.R
            x <- cmdscale(
                as.matrix(as.dist(1 - mtx_similarity)),
                eig = TRUE,
                k = min(c(nrow(mtx_similarity)-1, 15))
            )
            # Calculate a UMAP
            # Found this works better on the PCs rather than calling UMAP
            # directly
            umap_config <- umap::umap.defaults
            # umap_config$n_neighbors <- min(c(nrow(mtx_similarity), 15))
            umap_config$n_neighbors <- min(c(nrow(x$points), 15))
            umap_fit <- umap::umap(
                # scale(as.matrix(as.dist(1 - mtx_similarity))),
                x$points,
                config = umap_config
            )
            umap_df <- umap_fit$layout %>%
                as.data.frame() %>%
                dplyr::rename(
                    term_UMAP1 = "V1",
                    term_UMAP2 = "V2"
                ) %>%
                as.data.frame()
            df_plt <- cbind(
                data.frame(
                    "term_PC1" = x$points[,1],
                    "term_PC2" = x$points[,2]
                ),
                umap_df[
                    match(rownames(x$points), rownames(umap_df)),
                    colnames(umap_df)
                ],
                reducedTerms[
                    match(rownames(x$points), reducedTerms$term_id),
                    colnames(reducedTerms)
                ]
            )
            df_plt <- merge(
                df_plt,
                df_gsea_ij,
                all.x = TRUE,
                all.y = FALSE
            )

            # plot the data
            plt <- ggplot2::ggplot(df_plt, ggplot2::aes(
                color = as.factor(cluster),
                x = term_UMAP1,
                y = term_UMAP2
            ))
            if (term_summary_type_i != "qvalue_bh_less_than_0pt05") {
                plt <- plt + ggplot2::geom_point(ggplot2::aes(
                        alpha = qvalue_bh < 0.05,
                        size = score
                    )
                )
                plt <- plt + ggplot2::scale_alpha_discrete(
                    name = "FDR<5%",
                    range = c(0.3, 0.9)
                )
            } else {
                plt <- plt + ggplot2::geom_point(ggplot2::aes(
                        size = score
                    ),
                    alpha = 0.5
                )
            }
            plt <- plt + ggplot2::scale_color_discrete(guide = "none")
            # plt <- plt + ggplot2::scale_size_continuous(
            #     guide="none",
            #     range=c(0, 25)
            # )
            plt <- plt + ggplot2::labs(
                x = "UMAP1",
                y = "UMAP2",
                size = bquote(-log[10](italic(P))),
                title = paste0(
                    unique(df_plt$cell_label),
                    " ",
                    df_key_id
                )
            )
            plt <- plt + ggplot2::theme_bw()
            plt <- plt + ggplot2::theme(
                axis.text.x = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank()
            )
            # Add label
            labelSize <- 2
            plt <- plt + ggrepel::geom_label_repel(
                ggplot2::aes(label = parent_term_id),
                data = subset(df_plt, term_id == df_plt$parent_term_id),
                box.padding = grid::unit(1, "lines"),
                max.overlaps = Inf,
                size = labelSize
            )
            print(plt)

            df_results_list[[j]] <- df_plt

            # TODO: On the bases of GO term semantic similarity, calculate the
            # similarity between any two gene lists.
            #
            # NOTE: This could also be done based on disease ontology similarity
            # DOSE::geneSim(g1[1], g2[1], measure="Wang", combine="BMA")
            #
            # NOTE: This could also be done based on MESH terms linked to genes through
            #       textual analysis
            # library(AnnotationHub)
            # library(MeSHDbi)
            # ah <- AnnotationHub(localHub=FALSE)
            # hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
            # file_hsa <- hsa[[1]]
            # db <- MeSHDbi::MeSHDb(file_hsa)
            # hsamd <- meshes::meshdata(
            #     db, category='A', computeIC=T, database="gendoo")
            # meshes::geneSim("241", "251", semData=hsamd,
            #     measure="Wang", combine="BMA")
            #
            # similarity_via_go_semantic <- GOSemSim::clusterSim(
            #     annot_data[['REACTOME_ACTIVATION_OF_NF_KAPPAB_IN_B_CELLS']],
            #     annot_data[['REACTOME_REGULATION_OF_COMPLEMENT_CASCADE']],
            #     semData=semdata,
            #     measure="Rel"
            # )

            # # Calculate a similarity matrix between terms
            # semdata <- GOSemSim::godata(
            #     OrgDb="org.Hs.eg.db",
            #     keytype="SYMBOL",
            #     ont="BP"
            # )
            # simMatrix <- rrvgo::calculateSimMatrix(
            #     go_analysis$ID,
            #     semdata=semdata,
            #     method="Rel"
            # )

        } # End inner loop

        # Close the pdf output
        dev.off()

        # Save the data matrix of these results
        df_results <- data.table::rbindlist(df_results_list)
        fh <- paste0(
            param[["output_file_tag"]],
            "-",
            term_summary_type_i,
            ".tsv",
            sep = ""
        )
        gzfh <- gzfile(paste0(fh, ".gz"), "w", compression = 9)
        write.table(
            df_results,
            gzfh,
            row.names = FALSE,
            col.names = TRUE,
            quote = TRUE,
            sep = "\t",
            na = ""
        )
        close(gzfh)

        #return(1)
    } # end outer foreach loop


    # stop the cluster
    # if (param[["parallel_type"]] != "domc") {
    #     parallel::stopCluster(cl)
    # }
    time_end <- Sys.time()
    if (param[["verbose"]]) {
        message(paste0(
            "\nForeach loop execution time", " [", SCRIPT_NAME, "]:\t",
            difftime(time_end, time_init, units = c("hours")),
            " hours."
        ))
    }


    # df_list_results[["df_harmonized_ss"]] <- data.table::rbindlist(
    #     list(df_ss, df_list_results[["df_harmonized_ss"]]),
    #     use.names = TRUE,
    #     fill = TRUE
    # )

    # add param to return list of results
    #df_list_results[["param"]] <- param


    return(0)
}


main <- function() {
    # run analysis
    run_time <- system.time(df_results <- command_line_interface())
    message(paste0(
        "Analysis execution time", " [", SCRIPT_NAME, "]:\t",
        run_time[["elapsed"]]/3600, # proc.time sec to hours
        " hours."
    ))

    return(0)
}


dev <- function() {

    param <- list()
    param[["verbose"]] <- TRUE
    #param[["threads"]] <- 2
    param[["output_file_tag"]] <- "test"
    param[["term_distance_metric"]] <- "kappa"
    param[["clustering_method"]] <- "louvain"
    param[["gsea_results"]] <- c(
        "~/temp/data_temp/time_point_smoothed_glucose_concentration_merged-gsea_results.tsv.gz"
    )
    param[["gsea_results"]] <- c(
        "~/temp/dge/disease_status_merged-gsea_results.tsv.gz"
    )
    param[["gsets_gene_matrix"]] <- c(
        "~/repo/sc_nf_diffexpression/data/gene_set_gene_matrix-fixed.tsv.gz"
    )
    run_analysis(param)

    # for profiling
    # suppressPackageStartupMessages(library(profvis))
    #
    # base <- paste0("profile-", gsub(" ", "_", Sys.time()))
    # # run the analysis to understand memory usage
    # prof <- profvis::profvis(command_line_interface())
    # saveRDS(prof, paste0(base, ".Rds.gz"), compress = TRUE)
    # #prof <- readRDS(paste0(base, ".Rds.gz"))
    # #print(prof)
    #
    # # only run if snow
    # prof <- snow::snow.time(run_analysis())
    # pdf(file = paste0(base, ".pdf"), height = 5, width = 6)
    #     print(plot(prof))
    # dev.off()
}


main()
#dev()
