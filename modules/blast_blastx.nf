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
    // path "*.blastx_fmt14.xml" , emit: xml
    path "*.blastx.tsv"       , emit: tsv
    path "*.predicted.fasta"  , emit: fasta
    path "*.version"          , emit: version

    script:
    def prefix = query.baseName
    """
    blastx -num_threads ${task.cpus} \\
        -query $query \\
        -db $db \\
        $options.args \\
        -outfmt "6 qseqid qseq" \\
        -out ${prefix}.blastx.tsv
    awk '{
        seq = ""
        i = 1
        if (seq == $1) {
            i++
        } else {
            i = 1
            seq = $1
        }
        gsub(/[-X*]/,"",$2)
        print ">"$1"_"i"\n"$2
    }' ${prefix}.blastx.tsv > ${prefix}.predicted.fasta
    blastx -version | sed -e '/^blastx:/!d; s/^.*blastx: //' > blastx.version
    """

}
