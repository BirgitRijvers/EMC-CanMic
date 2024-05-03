process KRAKENBIOM_INDIVIDUAL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kraken-biom:1.2.0--pyh5e36f6f_0':
        'biocontainers/kraken-biom:1.2.0--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(kreport)
    val(tool) // To differentiate between outputs when using module twice in same workflow

    output:
    path "*_classified.biom" , emit: biom
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    """
    kraken-biom \\
        $kreport \\
        -o ${tool}_${prefix}_classified.biom \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakenbiom: \$(kraken-biom --version | sed 's/^.*version //; s/,.*\$//')
    END_VERSIONS
    """

    // stub:
    // def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    // // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    // //               Have a look at the following examples:
    // //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    // //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    // """
    // touch ${prefix}.bam

    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     krakenbiom: \$(samtools --version |& sed '1!d ; s/samtools //')
    // END_VERSIONS
    // """
}
