#!/usr/bin/env nextflow


def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}


process get_cell_label_list {
    // Get all of the cell labels in an anndata file
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    input:
        path(anndata)
        val(anndata_cell_label)

    output:
        path("cell_labels.csv", emit: cell_labels)

    script:
        runid = random_hex(16)
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "get_cell_label_list: ${process_info}"
        010-get_cell_label_list.py \
            --h5_anndata ${anndata} \
            --cell_label ${anndata_cell_label}
        """
}


process run_differential_expression {
    // Run differential expression
    // ------------------------------------------------------------------------
    label 'long_job'
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(anndata)
        val(cell_label_column)
        val(experiment_id)
        val(mean_cp10k_filter)
        each cell_label
        each model

    output:
        val(outdir, emit: outdir)
        tuple(
            val(runid), //need random hex to control grouping
            val(variable_target),
            val(cell_label),
            val(formula_clean),
            val(model.method),
            path("${outfile}_unfiltered-de_results.tsv.gz"),
            path("${outfile}_filtered-de_results.tsv.gz"),
            val(outdir),
            optional: true,
            emit: results
        )
        path("*_ruvseq_factors.tsv.gz") optional true
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        // on first call, cell_label comes in array-format. Need to check if
        // list because Nextflow has terrible retry behavior where it passes
        // as string on second call
        cell_label = cell_label instanceof List ? cell_label[0] : cell_label
        prop_cov_col = model.proportion_covariate_column
        formula = "${model.formula}"
        formula_clean = "${model.formula}".replaceAll("_", "")
        formula_clean = "${formula_clean}".replaceAll(" ", "_")
        formula_clean = "${formula_clean}".replaceAll("~", "")
        formula_clean = "${formula_clean}".replaceAll("\\+", "_plus_")
        formula_clean = "${formula_clean}".replaceAll("_", "")
        formula_clean = "${formula_clean}".replaceAll("\\)", "")
        formula_clean = "${formula_clean}".replaceAll("\\(.*\\|", "ra_")
        formula_clean = "${formula_clean}".replaceAll("I\\(", "")
        formula_clean = "${formula_clean}".replaceAll("\\^", "power")
        formula_clean = "${formula_clean}".replaceAll("\\.", "pt")
        formula_clean = "${formula_clean}".replaceAll("\\/", "_div_")
        variable_target = "${model.variable_target}"
        variable_target_clean = "${model.variable_target}".replaceAll(
            "\\)", ""
        )
        variable_target_clean = "${variable_target_clean}".replaceAll(
            "I\\(", ""
        )
        variable_target_clean = "${variable_target_clean}".replaceAll(
            "\\^", "power"
        )
        variable_target_clean = "${variable_target_clean}".replaceAll(
            "\\.", "pt"
        )
        // Get optional flags and add to formula if needed
        cmd__options = ""
        if (model.pre_filter_genes) {
            cmd__options = "--pre_filter_genes"
        }
        if (model.include_proportion_covariates) {
            cmd__options = "${cmd__options} --include_proportion_covariates"
            formula_clean = "${formula_clean}__proportion_covs-${prop_cov_col}"
        }
        if (model.ruvseq) {
            cmd__options = "${cmd__options} --run_ruvseq"
            formula_clean = "${formula_clean}__ruvseq-ngenes=${model.ruvseq_n_empirical_genes}"
            formula_clean = "${formula_clean}_min_pvalue=${model.ruvseq_min_pvalue}"
            formula_clean = "${formula_clean}_kfactors=${model.ruvseq_k}"
        }
        outdir = "${outdir_prev}/differential_expression/${variable_target_clean}"
        outdir = "${outdir}/cell_label=${cell_label}"
        outdir = "${outdir}/method=${model.method}___formula=${formula_clean}"
        // Sort out any variables that need to be cast
        cmd__varcast = ""
        if (model.variable_discrete != "") {  // add disc cov call
            cmd__varcast = "${cmd__varcast} --discrete_variables \"${model.variable_discrete}\""
        }
        if (model.variable_continuous != "") {  // add contin cov call
            cmd__varcast = "${cmd__varcast} --continuous_variables \"${model.variable_continuous}\""
        }
        // Make discrete levels command
        cmd__levels = ""
        if (model.variable_discrete_level != "") {  // add disc cov call
            cmd__levels = "--discrete_levels \"${model.variable_discrete_level}\""
        }
        outfile = "cell_label__${cell_label}"
        // Finally get the correct script
        base_method = model.method.split("::")[0]
        if (base_method == "edger") {
            method_script = "011-run_edger.R"
        } else if (base_method == "deseq") {
            method_script = "011-run_deseq.R"
        } else if (base_method == "mast") {
            method_script = "011-run_mast.R"
        }
        // Details on process
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "run_differential_expression: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        generate_experiment_covariates.py \
            --h5_anndata ${anndata} \
            --experiment_id "${experiment_id}" \
            --proportion_covariate_column "${prop_cov_col}" \
            --out_file tmp_anndata.h5ad
        convert_h5ad_R.py \
            --h5ad_file tmp_anndata.h5ad \
            --variable_target "${variable_target}" \
            --cell_label_column "${cell_label_column}" \
            --cell_label "${cell_label}" \
            --output_dir de_input
        011-run_differential_expression.R \
            --input_dir de_input \
            --cell_label_column "${cell_label_column}" \
            --cell_label "${cell_label}" \
            --experiment_key "${experiment_id}" \
            --formula "${model.formula}" \
            --variable_target "${variable_target}" \
            --method "${model.method}" \
            --method_script $baseDir/bin/${method_script} \
            --mean_cp10k_filter ${mean_cp10k_filter} \
            --ruvseq_n_empirical_genes ${model.ruvseq_n_empirical_genes} \
            --ruvseq_min_pvalue ${model.ruvseq_min_pvalue} \
            --ruvseq_k_factors ${model.ruvseq_k} \
            --out_file "${outfile}" \
            --cores_available ${task.cpus} \
            ${cmd__varcast} \
            ${cmd__levels} \
            ${cmd__options}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        rm tmp_anndata.h5ad
        """
}


process plot_dge_results {
    // Plot DGE results
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        tuple(
            val(runid),
            val(variable_target),
            val(cell_label),
            val(formula_clean),
            val(method),
            file(de_results_unfiltered),
            file(de_results_filtered),
            val(outdir_prev)
        )

    output:
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = outdir_prev
        outfile = "cell_label__${cell_label}"
        // Details on process
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_dge_results: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        012-plot_dge_results.R \
            --input_file ${de_results_unfiltered} \
            --target_var 'coef_value' \
            --sig_var 'qvalue_bh_percelltype' \
            --sig_threshold 0.05 \
            --sig_label 'FDR <= 0.05' \
            --out_file '${outfile}_unfiltered' \
            --verbose
        012-plot_dge_results.R \
            --input_file ${de_results_filtered} \
            --target_var 'coef_value' \
            --sig_var 'qvalue_bh_percelltype' \
            --sig_threshold 0.05 \
            --sig_label 'FDR <= 0.05' \
            --out_file '${outfile}_filtered' \
            --verbose
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

process serialize_de_files {
    // Serializes known markers for analysis
    // ------------------------------------------------------------------------
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    input:
        tuple(
            val(runid),
            val(variable_target),
            val(cell_label),
            val(formula),
            val(method),
            file(de_results_unfiltered),
            file(de_results_filtered),
            val(outdir_prev)
        )

    output:
        tuple(
            val(runid),
            val(variable_target),
            val(cell_label),
            val(formula),
            val(method),
            path("${runid}-${de_results_unfiltered}"),
            path("${runid}-${de_results_filtered}"),
            val(outdir),
            emit: results
        )

    script:
        //runid = random_hex(16)
        outdir = outdir_prev
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "serialize_de_files: ${process_info}"
        cp ${de_results_unfiltered} ${runid}-${de_results_unfiltered}
        cp ${de_results_filtered} ${runid}-${de_results_filtered}
        """
}

process merge_de_dataframes {
    // Merge resulting dataframes
    // NOTE: if this function is called more than once (e.g., first time fails)
    //       then nextflow bug will result in
    //       val(result_keys) and file(result_paths) being empty and the
    //       pipeline will fail.
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        val(run_correction)
        val(correction_config)
        tuple(
            val(condition),
            val(result_keys),
            file(result_paths)
        )

    output:
        val(outdir, emit: outdir)
        tuple(
            val(condition),
            path("${outfile}-de_results.tsv.gz"),
            emit: merged_results
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/differential_expression/${condition}/"
        outfile = "${condition}_merged"
        result_keys = result_keys.join(",")
        result_paths = result_paths.join(",")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "merge_de_dataframes: ${process_info}"
        echo "publish_directory: ${outdir}"
        sleep 5m
        merge_dataframes.py \
            --dataframe_keys '${result_keys}' \
            --dataframe_paths '${result_paths}' \
            --output_file  '${outfile}-de_results.tsv.gz'
        012-correct_pvals.R \
            --de_results '${outfile}-de_results.tsv.gz' \
            --grouping_cols 'de_method,formula_passed,coef_value,include_cell_proportions' \
            --output_file '${outfile}-de_results.tsv.gz'
        if [ "${run_correction}" = true ]; then
            012-ihw_correction.R \
                --de_results '${outfile}-de_results.tsv.gz' \
                --covariates '${correction_config.covariates}' \
                --alpha ${correction_config.alpha} \
                --grouping_cols 'de_method,formula_passed,coef_value,include_cell_proportions' \
                --output_file '${outfile}-de_results.tsv.gz'
        fi
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process plot_merged_dge {
    // Generate plots from the merged data frames to evaluate
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        tuple(
            val(condition),
            path(merged_df)
        )
        each mean_expression_filter

    output:
        val(outdir, emit: outdir)
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/differential_expression/${condition}/"
        outfile = "${condition}_merged_de"
        // script automatically adds expression filter
        // outfile = "${condition}-mean_expr_filt__${mean_expression_filter}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_merged_dge: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        013-compare_de_results.py \
            --dataframe ${merged_df} \
            --columns_to_compare de_method,formula_passed,include_cell_proportions \
            --mean_expression_filter ${mean_expression_filter} \
            --output_file '${outfile}'
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process run_fgsea {
    // Run fGSEA for each DE result
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script
    label 'long_job'

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        tuple(
            val(runid),
            val(condition_column),
            val(cell_label),
            val(covariate_columns),
            val(method),
            file(de_results_unfiltered),
            file(de_results_filtered),
            val(outdir_prev)
        )
        each model
        each signed
        file(gene_matrix)
        file(gene_info)

    output:
        tuple(
            val(runid),
            val(condition_column),
            val(cell_label),
            val(covariate_columns),
            val(method),
            path("${outfile}-gsea_results.tsv.gz"),
            val(fgsea_key),
            val(outdir),
            emit: results
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true
        path("plots/*.jpg") optional true

    script:
        runid = random_hex(16)
        // Generate key to represent run to group results
        fgsea_key = "sample_size=${model.sample_size}"
        fgsea_key = "${fgsea_key}::score_type=${model.score_type}"
        fgsea_key = "${fgsea_key}::min_set_size=${model.min_set_size}"
        fgsea_key = "${fgsea_key}::max_set_size=${model.max_set_size}"
        fgsea_key = "${fgsea_key}::eps=${model.eps}"
        fgsea_key = "${fgsea_key}::db=${model.database}"
        fgsea_key = "${fgsea_key}::${signed}"
        fgsea_key = fgsea_key.replaceAll(",", "_")
        outdir = outdir_prev
        outdir = "${outdir}/fgsea"
        outdir = "${outdir}-${fgsea_key}"
        outfile = "cell_label__${cell_label}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        cmd_null = "" // Don't pass null args in
        if (signed == "unsigned") {
            cmd_null = "${cmd_null} --unsigned_ranking"
        }
        """
        echo "run_fGSEA: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        015-run_fgsea.R \
            --de_results ${de_results_filtered} \
            --group_var 'coef_value' \
            --ranking_var 'test_statistic' \
            --sample_size ${model.sample_size} \
            --score_type ${model.score_type} \
            --min_set_size ${model.min_set_size} \
            --max_set_size ${model.max_set_size} \
            --eps '${model.eps}' \
            --gsets_gene_matrix ${gene_matrix} \
            --gsets_info_file ${gene_info} \
            --database '${model.database}' \
            --n_cores ${task.cpus} \
            --output_file '${outfile}' \
            --verbose \
            ${cmd_null}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        mv *jpg plots/ 2>/dev/null || true
        """
}


process serialize_gsea_files {
    // Serializes known markers for analysis
    // ------------------------------------------------------------------------
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    input:
        tuple(
            val(runid),
            val(condition_column),
            val(cell_label),
            val(covariate_columns),
            val(method),
            path(results_file),
            val(run_key),
            val(outdir_prev)
        )

    output:
        tuple(
            val(runid),
            val(condition_column),
            val(cell_label),
            val(covariate_columns),
            val(method),
            path("${runid}-${results_file}"),
            val(run_key),
            val(outdir),
            emit: results
        )

    script:
        //runid = random_hex(16)
        outdir = outdir_prev
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "serialize_gsea_results: ${process_info}"
        echo "publish_directory: ${outdir}"
        cp ${results_file} ${runid}-${results_file}
        """
}


process merge_gsea_dataframes {
    // Merge resulting dataframes from GSEA
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        tuple(
            val(condition),
            val(result_keys),
            file(result_paths)
        )

    output:
        val(outdir, emit: outdir)
        tuple(
            val(condition),
            path("${outfile}-gsea_results.tsv.gz"),
            emit: merged_results
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/differential_expression/${condition}/"
        outfile = "${condition}_merged"
        result_keys = result_keys.join(",")
        result_paths = result_paths.join(",")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "merge_gsea_dataframes: ${process_info}"
        echo "publish_directory: ${outdir}"
        merge_dataframes.py \
            --dataframe_keys '${result_keys}' \
            --dataframe_paths '${result_paths}' \
            --output_file '${outfile}-gsea_results.tsv.gz'
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process plot_gsea_results {
    // Generate plots from the merged data frames to evaluate
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        tuple(
            val(condition),
            path(merged_df)
        )

    output:
        val(outdir, emit: outdir)
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/differential_expression/${condition}/"
        outfile = "${condition}_merged_gsea"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_gsea_results: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        016-compare_gsea_results.py \
            --dataframe ${merged_df} \
            --columns_to_compare de_method,gsea_method,formula_passed,signed_ranking,include_cell_proportions \
            --output_file '${outfile}'
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process summarize_gsea_results {
    // Summarize merged GSEA data frame
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        tuple(
            val(condition),
            path(merged_df)
        )
        each params

    output:
        val(outdir, emit: outdir)
        path("${outfile}-*.tsv.gz")
        path("*.pdf")

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/differential_expression/${condition}/"
        outfile = "${condition}_merged_gsea_summary"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "summarize_gsea_results: ${process_info}"
        echo "publish_directory: ${outdir}"
        017-summarize_fgsea.R \
            --gsea_results ${merged_df} \
            --gsets_gene_matrix '$baseDir/data/gene_set_genes.tsv.gz' \
            --term_distance_metric ${params.distance_metric} \
            --clustering_method ${params.clustering_method} \
            --output_file_tag '${outfile}'
        """
}


process run_goenrich {
    // Run GO Enrich on DGE results
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        tuple(
            val(condition),
            path(merged_df)
        )
        each params

    output:
        val(outdir, emit: outdir)
        path("${outfile}*.tsv.gz")
        path("*.pdf")

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/differential_expression/${condition}/"
        outfile = "${condition}_goenrich_results"
        pathways = params.go_terms.replace(",", " ")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "run_goenrich: ${process_info}"
        echo "publish_directory: ${outdir}"
        for pathway in ${pathways}; do
          015-run_goenrich.R \
              --dge_results ${merged_df} \
              --go_ontology \$pathway \
              --clustering_method ${params.clustering_method} \
              --output_file_tag "${outfile}__\${pathway}_terms" \
              --verbose
        done
        """
}


workflow wf__differential_expression {
    take:
        outdir
        anndata
        anndata_cell_label
        experiment_key
        model
        de_merge_config
        de_plot_config
        goenrich_config
        gsea_config
    main:
        // Get a list of all of the cell types
        get_cell_label_list(
            anndata,
            anndata_cell_label
        )
        // For each cell type compute differential expression for that cell
        // type
        cell_labels = get_cell_label_list.out.cell_labels
            .splitCsv(header: false, sep: ',')
            // .map{row -> tuple(
            //     row[0]
            // )}

        // Run DGE
        run_differential_expression(
            outdir,
            anndata,
            anndata_cell_label,
            experiment_key,
            model.mean_cp10k_filter,
            // '1',  // just run on first cluster for development
            cell_labels,  // run for all clusters for run time
            model.value
        )

        de_results = run_differential_expression.out.results

        // Plot results from each DGE run
        plot_dge_results(
            de_results
        )

        // Serialize input files to prep for merge
        serialize_de_files(
            de_results
        )

        // Group results by formula for merge
        de_results_merged = serialize_de_files.out.results
            .reduce([:]) { map, tuple ->
                def dataframe_key = "cell_label=" + tuple[2]
                dataframe_key += "::formula=" + tuple[3].replaceAll(
                    ",",
                    "-"
                )
                dataframe_key += "::method=" + tuple[4]

                def map_key = tuple[1] // structure map by condition
                def key_list = map[map_key]
                // Right now, we only want to carry down the filtered data
                if (!key_list) {
                    key_list = [[dataframe_key, tuple[6]]]
                } else {
                    key_list.add([dataframe_key, tuple[6]])
                }
                map[map_key] = key_list
                return(map)
            }
            .flatMap()
            .map {  entry ->
                combined_data = [entry.key, [], []]
                entry.value.each {
                    combined_data[1].add(it[0])
                    combined_data[2].add(it[1])
                }
                return(combined_data)
            }
        merge_de_dataframes(
            outdir,
            de_merge_config.ihw_correction.run_process,
            de_merge_config.ihw_correction.value,
            de_results_merged
        )
        // Basic plots of the differential expression results across all models
        plot_merged_dge(
            outdir,
            merge_de_dataframes.out.merged_results,
            de_plot_config.mean_expression_filter.value
        )

        // First run GO Enrich on merged DGE results
        if (goenrich_config.run_process) {
            run_goenrich(
                outdir,
                merge_de_dataframes.out.merged_results,
                goenrich_config.value
            )
        }

        // Run fGSEA on DE results
        gsea_results = null
        if (gsea_config.fgsea_parameters.run_process) {
            run_fgsea(
                de_results,
                gsea_config.fgsea_parameters.value,
                ["unsigned", "signed"],
                file("$baseDir/data/gene_set_genes.tsv.gz"),
                file("$baseDir/data/gene_set_info.tsv.gz")
            )

            gsea_results = run_fgsea.out.results
        }

        if (gsea_results != null) {
            // Serialize input files to prep for merge
            serialize_gsea_files(
                gsea_results
            )

            // Combine and compare all of the enrichment analysis results
            gsea_results_merged = serialize_gsea_files.out.results
                .reduce([:]) { map, tuple ->
                    def dataframe_key = "cell_label=" + tuple[2]
                    dataframe_key += "::covariates=" + tuple[3].replaceAll(
                        ",",
                        "-"
                    )
                    dataframe_key += "::method=" + tuple[4]
                    dataframe_key += "__gsea_params=" + tuple[6]

                    def map_key = tuple[1] // structure map by condition
                    def key_list = map[map_key]
                    if (!key_list) {
                        key_list = [[dataframe_key, tuple[5]]]
                    } else {
                        key_list.add([dataframe_key, tuple[5]])
                    }
                    map[map_key] = key_list
                    return(map)
                }
                .flatMap()
                .map {  entry ->
                    combined_data = [entry.key, [], []]
                    entry.value.each {
                        combined_data[1].add(it[0])
                        combined_data[2].add(it[1])
                    }
                    return(combined_data)
                }

            merge_gsea_dataframes(
                outdir,
                gsea_results_merged
            )

            // Basic plots of the GSEA results across all models
            plot_gsea_results(
                outdir,
                merge_gsea_dataframes.out.merged_results
            )

            // Summarize GSEA results across all models
            summarize_gsea_results(
                outdir,
                merge_gsea_dataframes.out.merged_results,
                gsea_config.gsea_summarize_parameters
            )
        }

    emit:
        cell_labels = get_cell_label_list.out.cell_labels
}
