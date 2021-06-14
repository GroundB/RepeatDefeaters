params.options = [:]

process BLAST_MAKEBLASTDB {

    conda (params.enable_conda ? 'bioconda::blast=2.10.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/blast:2.10.1--pl526he19e7b1_3'
    } else {
        container 'quay.io/biocontainers/blast:2.10.1--pl526he19e7b1_3'
    }

    input:
    path fasta

    output:
    path "blastdb"    , emit: db
    path "*.version"  , emit: version

    script:
    """
    makeblastdb \\
        -in $fasta \\
        $options.args
    mkdir blast_db
    mv ${fasta}* blast_db
    makeblastdb -version | sed -e '/^makeblastdb:/!d; s/^.*makeblastdb: //' > ${software}.version.txt
    """

}
