#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
  optparse::make_option(c('-i', '--GSEA_data_xml'),
                        type = 'character',
                        default = 'de_file',
                        help = 'XML containing GSEA pathways. Downloaded from
                        http://www.gsea-msigdb.org/gsea/downloads.jsp.'
  ),

  optparse::make_option(c('-o', '--output_file_base'),
                        type = 'character',
                        default = '',
                        help = 'Base output name.'
  )
)

parser <- optparse::OptionParser(
  usage = '%prog',
  option_list = optionList,
  description = paste0(
    'Calculates differentially expressed genes using MAST.'
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
suppressPackageStartupMessages(library(xml2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
################################################################################

################################ Functions #####################################

get_category <- function(gsea_cat, gsea_subcat) { 
  if (gsea_cat %in% c('H')) {
    return('h')
  } else if (gsea_subcat %in% c('C2_CGP', 'CGP')) {
    return('c2.cgp')
  } else if (gsea_subcat %in% c('C2_CP:BIOCARTA', 'CP:BIOCARTA')) {
    return('c2.cp.biocarta')
  } else if (gsea_subcat %in% c('CP:KEGG')) {
    return('c2.cp.kegg')
  } else if (gsea_subcat %in% c('C2_CP:REACTOME', 'CP:REACTOME')) {
    return('c2.cp.reactome')
  } else if (gsea_subcat %in% c('CP', 'CP:PID')) {
    return('c2.cp')
  } else if (gsea_subcat %in% c('CGN')) {
    return('c4.cgn')
  } else if (gsea_subcat %in% c('CM')) {
    return('c4.cm')
  } else if (gsea_subcat %in% c('C5_GO:BP', 'GO:BP')) {
    return('c5.bp')
  } else if (gsea_subcat %in% c('C5_GO:CC', 'GO:CC')) {
    return('c5.cc')
  } else if (gsea_subcat %in% c('C5_GO:MF', 'GO:MF')) {
    return('c5.mf')
  } else if (gsea_cat %in% c('C6')) {
    return('c6.all')
  } else if (gsea_cat %in% c('C7')) {
    return('c7.all')
  } else {
    return(gsea_subcat)
  }
}
################################################################################


######################## Read Data & Manipulate ################################
output_file_base <- arguments$options$output_file_base


print('Reading in the data...')
xml <- xml2::read_xml(arguments$options$GSEA_data_xml)

# First, get gene list and create matrix
genes <- unique(unlist(lapply(
  xml2::xml_attr(xml2::xml_children(xml), 'MEMBERS_SYMBOLIZED'), 
  function(x) {
    return(strsplit(x, ','))
  }
)))
genes <- sort(genes, decreasing = F)

gene_mtx <- data.frame(row.names=genes)
gene_set_mtx <- NULL

# Iterate through the list and add elements
for (pathway in xml2::xml_children(xml)) {
  pathway_genes <- strsplit(xml2::xml_attr(pathway, 'MEMBERS_SYMBOLIZED'), ',')[[1]]
  if (length(pathway_genes) == 0) {
    next
  }
  
  # Get data
  pathway_name <- xml2::xml_attr(pathway, 'STANDARD_NAME')
  pathway_cat <- get_category(
    xml2::xml_attr(pathway, 'CATEGORY_CODE'),
    xml2::xml_attr(pathway, 'SUB_CATEGORY_CODE')
  )
  
  if (is.null(gene_set_mtx)) {
    gene_set_mtx <- data.frame(
      gene_set = pathway_name,
      category = pathway_cat
    )
  } else {
    gene_set_mtx <- rbind(gene_set_mtx, c(pathway_name, pathway_cat))
  }
  
  # Assign member genes to pathway
  gene_mtx[pathway_name] <- as.numeric(rownames(gene_mtx) %in% pathway_genes)
}


# Save file
gene_mtx <- tibble::rownames_to_column(gene_mtx, var='gene')
gz_file <- gzfile(sprintf('%s/gene_set_genes.tsv.gz', output_file_base),
                  'w',
                  compression = 9)
write.table(x=gene_mtx,
            file = gz_file,
            sep='\t',
            col.names=T,
            row.names=F,
            quote=F)
close(gz_file)

gz_file <- gzfile(sprintf('%s/gene_set_info.tsv.gz', output_file_base),
                  'w',
                  compression = 9)
write.table(x=gene_set_mtx,
            file = gz_file,
            sep='\t',
            col.names=T,
            row.names=F,
            quote=F)
close(gz_file)

################################################################################
