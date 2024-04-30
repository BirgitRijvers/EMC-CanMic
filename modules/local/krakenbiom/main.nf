process KRAKENBIOM {
    tag "$meta.id"
    // label 'process_single'
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kraken-biom:1.2.0--pyh5e36f6f_0':
        'biocontainers/kraken-biom:1.2.0--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(kreport)
    val(tool) // Usefull for when using module twice
    // path  kreports, stageAs: "?/*"
    // Input looks like: 
    // [['id':'1643_sub', 'single_end':false], '/home/birgit/work/25/7b0c04a7c6c86d6a0c41986e3b220f/1643_sub.kraken2.report.txt'],
    // [['id':'1312_sub', 'single_end':false], '/home/birgit/work/a9/a11d9f84b4669d5fda5ea7aed11ec1/1312_sub.kraken2.report.txt'],
    // [['id':'1592_sub', 'single_end':false], '/home/birgit/work/8e/a86503f289b00e60862e83643b86d9/1592_sub.kraken2.report.txt']

    output:
    path "*_classified.biom" , emit: biom
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def kreport_paths = kreport.collect { it[1] }.join(' ')
    // Needs to convert input to space separated list of only paths to reports (second part of input tuple)

    """
    kraken-biom \\
        $kreport \\
        $args \\
        -o ${prefix}_${tool}_classified.biom \\
        $args

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
