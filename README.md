
# Description

This pipeline performs differential gene expression and gene set enrichment in droplet single cell RNA-seq. As input, this pipeline takes an H5AD single cell experiment, modeled after the output of the complimentary pipeline [nf_scrna_qc](https://github.com/wtsi-hgi/nf_scrna_qc).

For differential gene expression, we support the following packages:
1. [MAST](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5)
2. [edgeR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/)
3. [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)

For gene set enrichment, we support the following packages:
1. [fGSEA](https://www.biorxiv.org/content/10.1101/060012v3)

For multiple correction, we support:
1. [Benjamini-Hochberg](https://www.jstor.org/stable/2346101?seq=1#metadata_info_tab_contents)
2. [IHW](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4930141/)

For a full description of possible parameters, see [README-params.md](docs/README-params.md).

For more information on running the pipeline, see:
* [Common problems](docs/README-common_problems.md)
* [Output directory structure](docs/README-outdir_structure.md)


# Quickstart

Quickstart for deploying this pipeline locally and on a high performance compute cluster.


## 1. Set up the environment

See [environment README](env/README.md) to set up environment. Once the environment is set up, activate the `conda` environment:

```bash
source activate sc_diff_expr
```

Alternatively, if using singularity or docker, one can pull the image from [letaylor/sc_nf_diffexpression:latest](https://hub.docker.com/layers/sc_nf_diffexpression/letaylor/sc_nf_diffexpression/latest/images/sha256-3e37b79b64aa5b899f0bbb683e206a777f9d4afc9216968a2fe5dac1785fbfe4?context=explore).


## 2. Prepare the input files

Generate and/or edit input files for the pipeline.

The pipeline takes as input:
1. **--file_anndata**:  H5AD file containing sequencing data. Required.
2. **-params-file**:  YAML file containing analysis parameters. Required.

Examples of these files can be found in `demo/`.


## 3. Run pipeline

**NOTE**: All input file paths should be full paths.

To run:
```bash
nextflow run \
  "/path/to/repo/dir/main.nf" \
  -profile "local" \ # could also be configuration for cluster
  --file_anndata "/path/to/data/data.h5ad" \
  -params-file "/path/to/config/params.yml"
```

Examples:
* [Local with docker](demo/run_differential_expression_demo-docker.sh)
* [Local with singularity](demo/run_differential_expression_demo-singularity.sh)
* [Cluster with conda](demo/run_differential_expression_demo-conda.sh)


Authors: Henry Taylor and Leland Taylor
