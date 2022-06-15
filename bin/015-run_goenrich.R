#!/usr/bin/env Rscript

SCRIPT_NAME <- "run_goenrich.R"

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
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(GOSemSim))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(clusterProfiler))

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
        optparse::make_option(c("--dge_results"),
            type = "character",
            help = paste0(
                "Differential gene expression results file.",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--go_ontology"),
            type = "character",
            default = "CC",
            help = paste0(
                "Ontology type. Either MF (Molecular Function)",
                "CC (Cellular Component)",
                "BP (Biological Process)",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--clustering_method"),
            type = "character",
            default = "binary_cut",
            help = paste0(
                "Method to cluster terms.",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--output_file_tag"),
            type = "character",
            default = "go_enrichment",
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
            "Runs GO enrichment analysis."
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

    # Set the parallel backend: fork, snow, domc
    #param[["parallel_type"]] <- "domc"

    # Set the output file
    output_file_tag <- paste0(
        param[["output_file_tag"]],
        "-go_",
        param[["go_ontology"]]
    )

    # Load the GO semnatic similarity matrix so we do not need to do this
    # multiple times
    go_similary <- GOSemSim::godata(
        "org.Hs.eg.db",
        keytype = "ENSEMBL",
        ont = param[["go_ontology"]]
    )

    # Read in the dge results file
    f <- param[["dge_results"]]
    df_dge <- data.frame(data.table::fread(
        cmd = paste("gunzip -c", f),
        sep = "\t",
        header = T,
        stringsAsFactors = F
    ))
    # Strip out cell type for the df_key so then we can iterate over all
    # models in the output file.
    df_dge$df_key_model <- mapply(
        function(x,y) { gsub(x, "", y) },
        df_dge$cell_label,
        df_dge$df_key
    )
    # Label each df model with a shorter identifier so we can write that
    # in plots and output
    model_id_map <- list()
    MODEL_ID_NUMBER <- 1
    for (key in unique(df_dge$df_key_model)) {
        model_id_map[[key]] <- paste0(
            "df_key_",
            sprintf("%03d", as.numeric(MODEL_ID_NUMBER))
        )
        MODEL_ID_NUMBER <- MODEL_ID_NUMBER + 1
    }
    df_dge$df_key_model_id <- unlist(model_id_map[df_dge$df_key_model])
    #df_dge <- subset(df_dge, df_key == head(df_dge$df_key, n = 1))
    df_dge <- subset(
        df_dge,
        cell_label %in% head(unique(df_dge$cell_label), n = 10)
    )

    # Loop over each df_key and:
    # (1) Run enrichment
    # (2) Make joint plots
    # (3) Run cell type specific plots
    for (model_key_i in unique(df_dge$df_key_model)) {
        df_dge_i <- subset(df_dge, df_key_model == model_key_i)

        # Prep the data for enrichment analysis
        df_dge_theme <- df_dge_i %>%
            dplyr::filter(!is.na(pvalue)) %>%
            dplyr::filter(qvalue_bh_percelltype < 0.05) %>%
            dplyr::mutate(
                group_direction = if_else(test_statistic > 0, "up", "down")
            )
        
        if (nrow(df_dge_theme) == 0) {                               
            print("No DGE results for this model. Exiting.")                                                                                      
            return(df_dge_i)
        }                                


        # https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html
        # Perform biological theme comparison.
        # This is only available for enrichment analysis not GSEA
        # NOTE: There are minimal size of genes annotated by Ontology term for
        # testing and maximal. That means that if there is an intersection
        # less or greater than the two values of the DGE hits and the genes
        # in the term, then there will be no test. By default, in
        # clusterProfiler, the cutoffs are 10 and 500
        results_enrich <- clusterProfiler::compareCluster(
            # gene ~ cell_label + group_direction,
            gene ~ cell_label,
            data = df_dge_theme,
            fun = "enrichGO",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            qvalueCutoff = 1,
            OrgDb = org.Hs.eg.db::org.Hs.eg.db,
            keyType = "ENSEMBL",
            ont = param[["go_ontology"]],
            # minGSSize = 10,
            # maxGSSize = 500,
            readable = FALSE # maps ENG to gene names NOTE: causes errors if T
        )

        # Save the full enrichment results
        # TODO: manipulate dataframe to look like other fgsea results
        df_enrich <- results_enrich@compareClusterResult
        df_enrich$df_key_model <- model_key_i
        df_enrich$df_key_model_id <- unique(df_dge_i$df_key_model_id)
        # Add back in the full model key with cell type listed
        df_enrich$df_key <- mapply(
            function(x,y) { gsub("cell_label=", paste0("cell_label=", x), y) },
            df_enrich$cell_label,
            df_enrich$df_key_model
        )
        fh <- paste0(
            output_file_tag,
            "-",
            unique(df_dge_i$df_key_model_id),
            ".tsv",
            sep = ""
        )
        gzfh <- gzfile(paste0(fh, ".gz"), "w", compression = 9)
        write.table(
            df_enrich,
            gzfh,
            row.names = FALSE,
            col.names = TRUE,
            quote = TRUE,
            sep = "\t",
            na = ""
        )
        close(gzfh)

        # clusterProfiler::cnetplot(results_enrich) # other plot we don't need
        # Add term similarities to the @termsim slot
        term_similarity <- enrichplot::pairwise_termsim(
             results_enrich,
             method = "Rel",
             semData = go_similary,
             showCategory = length(unique(
                 results_enrich@compareClusterResult$ID
             ))
        )

        # Make global network plots across all cell types at two cutoffs:
        # pvalue < 0.05 and qvalue < 0.05
        filters <- list(
            "pvalue_less_than_0pt05" = (
                term_similarity@compareClusterResult$pvalue < 0.05
            ),
            "qvalue_bh_less_than_0pt05" = (
                term_similarity@compareClusterResult$qvalue < 0.05
            )
        )
        for (filter_i_name in names(filters)) {
            if (param[["verbose"]]) {
                cat(fh, filter_i_name, "\n")
            }
            # Set the output file handle for all plots
            fh <- paste0(
                output_file_tag,
                "-",
                unique(df_dge_i$df_key_model_id),
                "-",
                filter_i_name,
                sep = ""
            )

            # Get the appropriate filter
            # NOTE: for qvalue sometimes there are NAs in the matrix
            filter_i <- filters[[filter_i_name]]
            filter_i[is.na(filter_i)] <- FALSE

            # Subset down to the data
            term_similarity_i <- term_similarity
            term_similarity_i@compareClusterResult <- (
                term_similarity_i@compareClusterResult[filter_i,]
            )
            term_similarity_i@termsim <- term_similarity_i@termsim[
                term_similarity_i@compareClusterResult$Description,
                term_similarity_i@compareClusterResult$Description
            ]

            # Save the results
            pdf(
                paste0(fh, "-emapplot.pdf", sep = ""),
                height = 25,
                width = 25
            )
                # NOTE: points and pies do not show when
                #       gene ~ cell_label + group_direction
                # NOTE: if showCategory too big, throws the following error:
                #       Error in seq.default(min(radius), max(radius),
                #       length.out = n) : 'from' must be a finite number
                plt <- clusterProfiler::emapplot(
                    term_similarity_i,
                    showCategory = min(20, length(unique(
                        term_similarity_i@compareClusterResult$ID
                    ))), # how many terms on the plot
                    node_label = "category",
                    # group_category = TRUE,
                    # group_legend = TRUE,
                    # with_edge = TRUE,
                    min_edge = 0.75,
                    cex_line = 0.01,
                    # cex_category = 3,
                    pie = "count" # make the size correspond to clusters
                )
                print(plt)

                plt <- clusterProfiler::emapplot(
                    term_similarity_i,
                    showCategory = 20, # how many terms on the plot
                    node_label = "category",
                    min_edge = 0.75,
                    cex_line = 0.01,
                    # cex_category = 1,
                    pie = "count" # make the size correspond to clusters
                )
                print(plt)

            dev.off()

            # Init the pdf to dump figures for each cell type
            pdf(
                file = paste0(fh, "-per_celltype_plots.pdf", sep = ""),
                height = 8,
                width = 8
            )
            # Init the final dataframe with reduced dims
            df_results_list <- list()

            # Make a shorthand of the dataframe that we are going to loop over
            #df_enrich_i <- subset(df_enrich, filter_i)
            df_enrich_i <- term_similarity_i@compareClusterResult
            df_enrich_i$df_key_model <- model_key_i
            df_enrich_i$df_key_model_id <- unique(df_dge_i$df_key_model_id)
            # Add back in the full model key with cell type listed
            df_enrich_i$df_key <- mapply(
                function(x,y) {
                    gsub("cell_label=", paste0("cell_label=", x), y)
                },
                df_enrich_i$cell_label,
                df_enrich_i$df_key_model
            )
            df_enrich_i$term_id <- df_enrich_i$ID # for code below

            # Now loop over each cell type
            for (cell_label_j in unique(df_enrich_i$cell_label)) {
                if (param[["verbose"]]) {
                    cat(fh, cell_label_j, "\n")
                }
                df_enrich_ij <- subset(df_enrich_i, cell_label == cell_label_j)
                if (nrow(df_enrich_ij) <= 2) {
                    message(paste0(
                        fh,
                        cell_label_j,
                        ":\tskipping because too few terms\n"
                    ))
                    next
                }

                # Get the semantic similarity between GO terms
                mtx_similarity <- simplifyEnrichment::GO_similarity(
                    df_enrich_ij$term_id,
                    ont = param[["go_ontology"]],
                    db = "org.Hs.eg.db",
                    measure = "Rel",
                    remove_orphan_terms = FALSE
                )
                #print(head(mtx_similarity))
                print("mtx_similarity")

                # Get a dataframe of GO term and cluster pairs
                term_clusters <- simplifyEnrichment::simplifyGO(
                        mtx_similarity,
                        method = param[["clustering_method"]],
                        plot = FALSE,
                        verbose = TRUE
                    ) %>%
                    dplyr::arrange(cluster) %>%
                    as.data.frame()
                rownames(term_clusters) <- term_clusters$id
                term_clusters$cluster <- paste0(
                    "cluster_",
                    sprintf("%02d", as.numeric(term_clusters$cluster))
                )
                #print(head(term_clusters))
                print("term_clusters")

                # Make a GO term word cloud
                # https://jokergoo.github.io/simplifyEnrichment/articles/word_cloud_anno.html
                term_clusters_wc_list <- list()
                for (cluster_i in unique(term_clusters$cluster)) {
                    term_clusters_wc_list[[cluster_i]] <- subset(
                        term_clusters, cluster == cluster_i
                    )$id
                }
                # Make a plot of each cluster as heatmap
                # NOTE: for simplifyEnrichment::anno_word_cloud function you
                # may need github version of simplifyEnrichment
                plt <- ComplexHeatmap::Heatmap(
                    mtx_similarity[term_clusters$id, term_clusters$id],
                    row_split = term_clusters$cluster,
                    show_row_names = FALSE,
                    show_column_names = FALSE,
                    row_title = NULL,
                    row_title_rot = 0,
                    right_annotation = ComplexHeatmap::rowAnnotation(
                        wc = simplifyEnrichment::anno_word_cloud_from_GO(
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
                    column_title = cell_label_j
                )
                print(plt)
                if (param[["verbose"]]) {
                    cat(fh, cell_label_j, ":\tsaved heatmap\n")
                }

                # Make a plot of the PCs
                # Reduce a set of terms based on their similarity and
                # encrichment pvalues.
                # Also select a term as representative of all other terms in a
                # cluster.
                #
                # Prep scores to obtain term summary for each term
                enrichment_scores <- setNames(
                    -log10(df_enrich_ij$pvalue),
                    df_enrich_ij$term_id
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
                umap_pcs <- FALSE
                umap_config <- umap::umap.defaults
                if (umap_pcs) {
                    umap_config$n_neighbors <- min(c(nrow(x$points), 15))
                    umap_fit <- umap::umap(
                        x$points,
                        config = umap_config
                    )
                } else {
                    umap_config$n_neighbors <- min(c(nrow(mtx_similarity), 15))
                    umap_fit <- umap::umap(
                        scale(as.matrix(as.dist(1 - mtx_similarity))),
                        config = umap_config
                    )

                }
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
                    df_enrich_ij,
                    all.x = TRUE,
                    all.y = FALSE
                )
                # Add the dataframe to the output list
                df_results_list[[cell_label_j]] <- df_plt
                if (param[["verbose"]]) {
                    cat(fh, ":\tsaved semantic similarity df\n")
                }

                term_summary_type_i = ""
                # plot the data
                plt <- ggplot2::ggplot(df_plt, ggplot2::aes(
                    color = as.factor(cluster),
                    x = term_PC1,
                    y = term_PC2
                ))
                if (term_summary_type_i != "qvalue_bh_less_than_0pt05") {
                    plt <- plt + ggplot2::geom_point(ggplot2::aes(
                            alpha = qvalue < 0.05,
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
                    title = cell_label_j
                )
                plt <- plt + ggplot2::theme_bw()
                plt <- plt + ggplot2::theme(
                    axis.text.x = ggplot2::element_blank(),
                    axis.text.y = ggplot2::element_blank()
                )
                # Add label
                labelSize <- 2
                plt <- plt + ggrepel::geom_label_repel(
                    ggplot2::aes(label = Description),
                    data = subset(df_plt, term_id == parent_term_id),
                    box.padding = grid::unit(1, "lines"),
                    max.overlaps = Inf,
                    size = labelSize
                )
                print(plt)
                if (param[["verbose"]]) {
                    cat(
                        fh, cell_label_j,
                        ":\tsaved semantic similarity pca/umap plot\n"
                    )
                }

                # Make a treemapPlot to better see the terms in each cluster
                rownames(df_plt) <- df_plt$ID
                reducedTerms$term <- df_plt[reducedTerms$term_id, ]$Description
                reducedTerms$parentTerm <- df_plt[
                    reducedTerms$parent_term_id,
                ]$Description
                rrvgo::treemapPlot(reducedTerms)
                if (param[["verbose"]]) {
                    cat(fh, cell_label_j, ":\tsaved treemap plot\n")
                }
            } # End per cell type loop

            # Clouse out the plots for all cell types
            dev.off()

            # Save the results
            df_results <- data.table::rbindlist(df_results_list)
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

        } # End filter pvalue_less_than_0pt05 or qvalue_bh_less_than_0pt05
    } # End model_key_i loop


    # stop the cluster
    # if (param[["parallel_type"]] != "domc") {
    #     parallel::stopCluster(cl)
    # }
    # time_end <- Sys.time()
    # if (param[["verbose"]]) {
    #     message(paste0(
    #         "\nForeach loop execution time", " [", SCRIPT_NAME, "]:\t",
    #         difftime(time_end, time_init, units = c("hours")),
    #         " hours."
    #     ))
    # }

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
    param[["clustering_method"]] <- "binary_cut"
    param[["go_ontology"]] <- "MF"
    param[["dge_results"]] <- fs::path_expand(
        "~/temp/dge/disease_status_merged-de_results.tsv.gz"
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
