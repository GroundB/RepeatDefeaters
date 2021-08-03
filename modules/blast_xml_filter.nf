params.options = [:]

process FILTER_BLAST_XML {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img'
    } else {
        container 'biocontainers/biocontainers:v1.2.0_cv1'
    }

    input:
    // path blast_xml // -outfmt 14 = Multiple-file BLAST XML2
    path blast_tsv // -outfmt "6 qseqid qseq"
    val prefix // prefix - cclaro ?

    output:
    path "predicted.fasta", emit: fasta

    script:
    """
    #plus strand cleanup
    echo "perform cleanup on plus strand predictions"
    for xml in $blast_xml; do
        # hash : number of query sequences
        hash=$( grep -c qseq \$xml )
        # bname : capture query-title
        bname=$( grep query-title \$xml | sed "s/.*${prefix}/${prefix}/" | sed 's/<.*//g' )
        if [ \$hash != 0 ]; then
            echo \$bname \$hash
            for ((i=1;i<=\$hash;i++)); do
                # TODO:: Change PLUS in filenames
                echo '>'\${bname}_plus_qseq_\$i'' >> \${bname}_\${hash}.plus.pre
                # Extract protein sequence and remove certain chacters {-,X,*}
                awk '/qseq/{i++}i=='\$i'' \$xml | \\
                   sed -n '/hseq/q;p' | \\
                   sed s'/qseq>//g' | \\
                   sed s'/<//g' | \\
                   sed s'/\///g' | \\
                   sed s'/-//'g | \\
                   sed s'/\*//'g | \\
                   sed s'/X//'g | \\
                   tr -d " \t" >> \${bname}_\${hash}.plus.pre
            done
        fi
    done
    cat *.plus.pre > ${prefix}.predicted.fasta
    """
}
