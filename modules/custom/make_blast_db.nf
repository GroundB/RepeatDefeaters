process MAKE_BLAST_DB {

    conda (params.enable_conda ? 'bioconda::blast=2.10.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.10.1--pl526he19e7b1_3' :
        'quay.io/biocontainers/blast:2.10.1--pl526he19e7b1_3' }"

    input:
    path fasta

    output:
    path "*_blast_db"    , emit: db
    path "versions.yml"  , emit: versions

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: fasta.baseName
    """
    makeblastdb \\
        $args \\
        -in $fasta
    mkdir ${prefix}_blast_db
    mv ${fasta}* ${prefix}_blast_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        makeblastdb: \$(makeblastdb -version | sed -e '/^makeblastdb:/!d; s/^.*makeblastdb: //' )
    END_VERSIONS
    """

}
