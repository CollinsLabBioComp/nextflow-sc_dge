#!/bin/sh

# check nextflow is in PATH
nextflow -version

# check docker is in PATH
docker --version

# create and use local temporary directory
export TMPDIR="$(pwd)/tmp"
mkdir -p ${TMPDIR}

# Nextflow settings - one may need to set these
# export JAVA_HOME="/path/to/java/jre1.8.0_251"
# export JAVA_CMD="/path/to/java/jre1.8.0_251/bin/java"
# export NXF_OPTS="-Xms25G -Xmx25G"

# Run the DGE model demo
nextflow run \
  ../main.nf \
  -profile "local_docker" \
  --file_anndata "$(pwd)/demo_data.h5ad" \
  --output_dir "$(pwd)/results" \
  -params-file "$(pwd)/params-dge_test_type_demo.yml" \
  -resume

# Run the GO and GSEA demo
nextflow run \
  ../main.nf \
  -profile "local_docker" \
  --file_anndata "$(pwd)/demo_data.h5ad" \
  --output_dir "$(pwd)/results" \
  -params-file "$(pwd)/params-gsea_go_demo.yml" \
  -resume
