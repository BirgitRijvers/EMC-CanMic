process QIIME2_COREDIVERSITY {
    label 'process_single', 'error_ignore'
    conda "${moduleDir}/environment.yml"
    
    input:
    path freq_table
    val sampling_depth
    path metadata
    val name

    output:
    path "*rarefied_table.qza"       , emit: rarefied_table
    path "*observed_features.qza"    , emit: observed_features
    path "*shannon.qza"              , emit: shannon
    path "*evenness.qza"             , emit: evenness
    path "*jaccard_distance.qza"     , emit: jaccard_distance
    path "*bray_curtis_distance.qza" , emit: bray_curtis_distance
    path "*jaccard_pcoa.qza"         , emit: jaccard_pcoa
    path "*bray_curtis_pcoa.qza"     , emit: bray_curtis_pcoa
    path "*jaccard_plot.qzv"         , emit: jaccard_plot
    path "*bray_curtis_plot.qzv"     , emit: bray_curtis_plot
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    qiime diversity core-metrics \\
        $args \\
        --i-table ${freq_table} \\
        --p-sampling-depth ${sampling_depth} \\
        --m-metadata-file ${metadata} \\
        --p-n-jobs $task.cpus \\
        --o-rarefied-table ${name}_rarefied_table.qza \\
        --o-observed-features-vector ${name}_observed_features.qza \\
        --o-shannon-vector ${name}_shannon.qza \\
        --o-evenness-vector ${name}_evenness.qza \\
        --o-jaccard-distance-matrix ${name}_jaccard_distance.qza \\
        --o-bray-curtis-distance-matrix ${name}_bray_curtis_distance.qza \\
        --o-jaccard-pcoa-results ${name}_jaccard_pcoa.qza \\
        --o-bray-curtis-pcoa-results ${name}_bray_curtis_pcoa.qza \\
        --o-jaccard-emperor ${name}_jaccard_plot.qzv \\
        --o-bray-curtis-emperor ${name}_bray_curtis_plot.qzv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$(qiime --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """

    // stub:
    // def args = task.ext.args ?: ''
    
    // // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    // //               Have a look at the following examples:
    // //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    // //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    // """
    // touch ${prefix}.bam

    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     qiime2: \$(samtools --version |& sed '1!d ; s/samtools //')
    // END_VERSIONS
    // """
}
