process BLASTN {

    conda (params.enable_conda ? 'bioconda::blast=2.10.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.10.1--pl526he19e7b1_3' :
        'quay.io/biocontainers/blast:2.10.1--pl526he19e7b1_3' }"

    input:
    path query
    path db

    output:
    path "*.blastn.tsv"       , emit: tsv
    path "versions.yml"       , emit: versions

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: query.baseName
    """
    BLASTDB=\$( find -L ./ -name "*.ndb" | sed 's/.ndb//' )
    blastn -num_threads ${task.cpus} \\
        -query $query \\
        -db \$BLASTDB \\
        $args \\
        -out ${prefix}.blastn.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$( blastn -version | sed -e '/^blastn:/!d; s/^.*blastn: //' )
    END_VERSIONS
    """

}
