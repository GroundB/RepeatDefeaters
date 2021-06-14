params.options = [:]

process BLAST_BLASTX {

    conda (params.enable_conda ? 'bioconda::blast=2.10.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/blast:2.10.1--pl526he19e7b1_3'
    } else {
        container 'quay.io/biocontainers/blast:2.10.1--pl526he19e7b1_3'
    }

    input:
    val query
    path db

    output:
    path "*.blastx.txt" , emit: txt
    path "*.version"    , emit: version

    script:
    """
    blastx -num_threads ${task.cpus} \\
        -query $query \\
        -db $db \\
        $options.args \\
        -out ${prefix}.blastx.txt
    blastx -version | sed -e '/^blastx:/!d; s/^.*blastx: //' > blastx.version
    """

}
