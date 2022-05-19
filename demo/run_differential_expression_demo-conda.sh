#!/bin/sh

# Activate the conda environment
#source activate sc_diff_expr

# Remove old logs but not the most previous run
rm -r *html.*;
rm .nextflow.log.*;
rm flowchart.png.*;
rm trace.txt.*;

# Set repo details. Below we assume the repo is in ${HOME}/repo.
# TODO: User edit this section to fit workstation
export NF_PATH="/path/to/nextflow"
export REPO_MODULE="${HOME}/repo/sc_nf_diffexpression"
export STUDY_DIR="${HOME}/repo/sc_nf_diffexpression/demo"
export OUTPUT_DIR="$(pwd)/results"

# Nextflow settings
export NXF_HOME=${OUTPUT_DIR}
export NXF_WORK="${NXF_HOME}/.nextflow_work"
export NXF_TEMP="${NXF_HOME}/.nextflow_temp"
export NXF_CONDA_CACHEDIR="${NXF_HOME}/.nextflow_conda"
export NXF_SINGULARITY_CACHEDIR="${NXF_HOME}/.nextflow_singularity"

# Cluster specific settings
# TODO: User edit this section to fit workstation
export JAVA_HOME="/path/to/java/jre1.8.0_251"
export JAVA_CMD="/path/to/java/jre1.8.0_251/bin/java"
export NXF_OPTS="-Xms25G -Xmx25G"

# Uncomment this if get strange bus errors
# export NXF_OPTS="${NXF_OPTS} -Dleveldb.mmap=false" # No resume functionality

# To stop TclError
export QT_QPA_PLATFORM='offscreen'

# Run nextflow
# TODO: User edit to set the proper profile
${NF_PATH}/nextflow run \
    "${REPO_MODULE}/main.nf" \
     -profile "lsf" \
     --file_anndata "${STUDY_DIR}/demo_data.h5ad" \
     --output_dir "${OUTPUT_DIR}" \
     -params-file "${STUDY_DIR}/params-dge_test_type_demo.yml" \
     -with-report \
     -with-trace \
     -with-timeline \
     -with-dag flowchart.png \
     -resume
