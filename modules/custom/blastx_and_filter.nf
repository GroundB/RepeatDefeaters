process BLASTX_AND_FILTER {

    conda (params.enable_conda ? 'bioconda::blast=2.10.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.10.1--pl526he19e7b1_3' :
        'quay.io/biocontainers/blast:2.10.1--pl526he19e7b1_3' }"

    input:
    path query
    path db
    each strand

    output:
    path "*.blastx.tsv"       , emit: tsv
    path "*.predicted.fasta"  , emit: fasta
    path "versions.yml"       , emit: versions

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: query.baseName
    """
    BLASTDB=\$( find -L ./ -name "*.pdb" | sed 's/.pdb//' )
    blastx -num_threads ${task.cpus} \\
        -query $query \\
        -db \$BLASTDB \\
        $args \\
        -outfmt "6 qseqid qseq" \\
        -strand $strand \\
        -out ${prefix}.${strand}.blastx.tsv
    # Blast output is a two column table
    # awk: Format two column output to fasta
    #      Label sequence header with strand and
    #      accumulator to form unique sequence headers
    awk 'BEGIN {
        seq = ""
        i = 1
    }
    {
        if (seq == \$1) {
            i++
        } else {
            i = 1
            seq = \$1
        }
        gsub(/[-X*]/,"",\$2)
        print ">"\$1"_${strand}_qseq_"i"\\n"\$2
    }' ${prefix}.${strand}.blastx.tsv > ${prefix}.${strand}.predicted.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastx: \$( blastx -version | sed -e '/^blastx:/!d; s/^.*blastx: //' )
        awk   : \$( awk |& sed '1!d; s/.*\\(v[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/' )
    END_VERSIONS
    """

}
