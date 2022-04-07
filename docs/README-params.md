
# Parameters

## Grouping parameters

- `experiment_key_column`: Column in cell metadata used to generate pseudobulk datasets and calculate cell type proportions (combined with term `proportion_covariate_column`)
- `anndata_cell_label`: Column in cell metadata representing cell types to perform differential expression and gene set enrichment

## Differential expression parameters

Parameters are applied within each cell type as denoted by `anndata_cell_label`

- `mean_cp10k_filter`: Remove genes with mean counts per 10,000 (CP10K) expression <= `mean_cp10k_filter`
- `models`: List of configurations for differential expression methods to run
  - `method`: String in the format `package::resolution::model`. Possible values for each:
    * `package`: "mast", "edger", "deseq"
    * `resolution`: "singlecell" or "pseudobulk"
    * `model`: Options:
      * For MAST, "bayesglm", "glmer" (needed for random effect models), or "glm"
      * For edgeR, "glmQLFit" or "glmLRT"
      * For DESeq2, "glmGamPoi" (recommended), "parametric", "local", or "mean"
  - `formula`: Formula to model the gene expression. Terms should be columns in cell metadata (Ex: "~ sex + age + disease_status"). The pipeline also supports the following operations:
    * R formula functions, such as "I(age^2)" and interaction effects, denoted with ":" (e.g., "time_point:disease_status"). See: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/formula
    * Random effects as denoted by "(1|participant_id)"
  - `variable_target`: Term in `formula` to test (e.g., disease_status)
  - `variable_continuous`: Terms in `formula` that should be cast as continuous covariates
  - `variable_discrete`: Terms in `formula` that should be cast as discrete (i.e., categorical) covariates
  - `variable_discrete_level`: Reference information for discrete covariates, formatted as "cov_1::ref,alt_1,alt_2,...,alt_n;;cov_2::ref,alt_1"
  - `pre_filter_genes`: Logical (e.g., true or false) to apply `mean_cp10k_filter` before or after performing differential expression
  - `proportion_covariate_column`: Column in cell metadata to calculate the proportion of cells from each experiment (defined by `experiment_key_column`) representing each value. For instance, if same as `anndata_cell_label`,  pipeline will calculate the proportion of each cell type for each experiment key.
  - `include_proportion_covariates`: Logical (e.g., true or false) to include proportions from `proportion_covariate_column` in `formula`
- `de_merge_config`: Configuration for merge settings
  - `ihw_correction`: Configuration for IHW correction
    - `covariates`: Comma-separated list of covariates to include in IHW correction (e.g., "cell_label,disease_status")
    - `alpha`: See [IHW](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4930141/) documentation
- `de_plot_config`: Parameters for plotting differential expression results
  - `mean_expression_filter`: List of mean expression thresholds to drop for plots for each group in `anndata_cell_label`. For example: if gene A expression is 0 counts in cluster 1 and 10 in cluster 2, it will be dropped from cluster 1 but not cluster 2.
- `fgsea_config`: Parameters for running gene set enrichment using fGSEA
  - `sample_size`: See [fGSEA](https://www.biorxiv.org/content/10.1101/060012v3) documentation
  - `score_type`: See [fGSEA](https://www.biorxiv.org/content/10.1101/060012v3) documentation
  - `value`: List of alternate configurations
    - `min_set_size`: See [fGSEA](https://www.biorxiv.org/content/10.1101/060012v3) documentation
    - `max_set_size`: See [fGSEA](https://www.biorxiv.org/content/10.1101/060012v3) documentation
    - `eps`: See [fGSEA](https://www.biorxiv.org/content/10.1101/060012v3) documentation
    - `database`: Comma-separated list of databases to test for enrichments. Detailed descriptions of databases can be found [here](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). Options:
      * `c2.cgp`: Chemical and genetic perturbations
      * `c2.cp.biocarta`: BioCarta
      * `c2.cp.kegg`: KEGG
      * `c2.cp.reactome`: Reactome
      * `c2.cp`: PID
      * `c5.bp`: GO biological process
      * `c5.cc`: GO cellular component
      * `c5.mf`: GO molecular function
      * `c6.all`: Oncogenic signatures
      * `c7.all`: Immunologic signatures
      * `all`: All gene sets (c2.cp.reactome, c2.cp.kegg, c5.bp, c5.cc, c5.mf)
