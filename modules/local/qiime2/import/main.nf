process QIIME2_IMPORT {
    label 'process_single','error_ignore'
    conda "${moduleDir}/environment.yml"

    input:
    path biom
    
    output:
    path "frequency.qza"                  , emit: frequency
    path "taxonomy.qza"                  , emit: taxonomy
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    

    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    """
    qiime tools import \\
        $args \\
        --input-path ${biom} \\
        --type FeatureTable[Frequency] \\
        --input-format BIOMV210Format \\
        --output-path frequency.qza

    qiime tools import \\
        $args \\
        --input-path ${biom} \\
        --type FeatureData[Taxonomy] \\
        --input-format BIOMV210Format \\
        --output-path taxonomy.qza

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
