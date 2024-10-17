process KRAKENBIOM_COMBINED {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kraken-biom:1.2.0--pyh5e36f6f_0':
        'biocontainers/kraken-biom:1.2.0--pyh5e36f6f_0' }"

    input:
    path kreports
    val(tool) // To differentiate between outputs when using module twice in same workflow

    output:
    path "*.biom"           , emit: biom
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    kraken-biom \\
        $kreports \\
        $args \\
        -o ${tool}_combined_classified.biom \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakenbiom: \$(kraken-biom --version | sed 's/^.*version //; s/,.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    touch ${tool}_combined_classified.biom

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakenbiom: \$(kraken-biom --version | sed 's/^.*version //; s/,.*\$//')
    END_VERSIONS
    """
}
