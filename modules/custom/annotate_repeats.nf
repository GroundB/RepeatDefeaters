process ANNOTATE_REPEATS {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path repeat_library            // Renamed repeat consensus library
    path pfam_table                // Pfam output (*.pfamtbl)
    path pfam_keyword_accession    // Pfam accessions from keyword search
    val prefix                     // prefix

    output:
    path "${prefix}.renamed.fasta"              , emit: fasta
    path "${prefix}.Unclassified_consensus_TEs" , emit: unclassified_with_te_domains
    path "${prefix}.consensus.both.strand"      , emit: unclassified_with_non_te_domains_both_strands
    path "versions.yml"                         , emit: versions

    script:
    """
    # Find unclassified consensus with TE domains
    for TBL in $pfam_table; do
        # grep #1    : Find unclassified consensus
        # grep #2    : which have TE domains
        # cut + uniq : and extract their id's
        grep -i "#unknown" "\$TBL" | \\
            tee "\${TBL}.unclassified" | \\
            grep -i -w -f $pfam_keyword_accession | \\
            tee -a ${prefix}.Unclassified_consensus_TEs | \\
            cut -f1 -d"#" | uniq > "\${TBL/.pfamtbl/.unclassified_ids}"
    done

    # Concatenate ids of consensus with TE domains from both strands
    cat *.unclassified_ids | uniq > ${prefix}.Unclassified_consensus_TEs.ids

    # Find unclassified consensus without TE domains
    for UNCLASSIFIED in *.unclassified; do
        # grep       : Remove consensus which have TE domains
        # cut + uniq : and extract their id's
        grep -v -f ${prefix}.Unclassified_consensus_TEs.ids "\$UNCLASSIFIED" | \\
            tee "\${UNCLASSIFIED}.TEpurged" | \\
            cut -f1 -d'#' | uniq > "\${UNCLASSIFIED}.TEpurged.ids"
    done

    # Use shell expansion to expand plus and minus strand files for unsorted inner join
    grep -f *.TEpurged.ids > ${prefix}.consensus.both.strand
    # ${prefix}.consensus.both.strand : Unclassified consensus sequences that have
    # non-TE domains detected in both strands.
    # These are tricky to annotate.

    # In consensus without TE domains, remove consensus with non-TE domains on both strands
    # (leaving consensus with non-TE domains on a single-strand)
    for TEPURGED in *.TEpurged; do
        # grep       : Remove consensus with non-TE domains on both strands
        # awk        : then remove consensus shorter than 100 amino acids
        # cut + uniq : and extract their id's
        grep -v -f ${prefix}.consensus.both.strand "\$TEPURGED" | \\
            awk '\$11 >= 100' | tee "\$TEPURGED.mono" | \\
            cut -f1 -d'#' | uniq > "\$TEPURGED.mono.ids"
    done

    # Make a copy of repeat library to be modified.
    cp $repeat_library ${prefix}.renamed.fasta

    # Rename repeat model based on strand evidence.
    for CONSENSUS in *.mono.ids; do
        # while       : for each consensus id
        # echo        : record id as renamed
        # NAMEHASH    : create a name suffix from the pfam domain table
        # OLDNAME     : find old name from repeat consensus library
        # sed         : replace Unknown with NAMEHASH
        while read -r SEQID; do
            echo "\$SEQID" >> ${prefix}.renamed
            NAMEHASH=\$( grep "\${SEQID}#" "\${CONSENSUS/.ids/}" | \\
                tr -s " " "\t" | cut -f7 | \\
                sort | uniq | \\
                paste -s -d '-' )
            OLDNAME=\$( grep "\${SEQID}#" $repeat_library | cut -c2- )
            sed -i "s|\$OLDNAME|\${OLDNAME%Unknown}\$NAMEHASH|g" ${prefix}.renamed.fasta
        done < "\$CONSENSUS"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk  : \$( awk  -W version |& head -n1 )
        cat  : \$( cat   --version |& head -n1 )
        cut  : \$( cut   --version |& head -n1 )
        grep : \$( grep  --version |& head -n1 )
        paste: \$( paste --version |& head -n1 )
        sed  : \$( sed   --version |& head -n1 )
        sort : \$( sort  --version |& head -n1 )
        tee  : \$( tee   --version |& head -n1 )
        uniq : \$( uniq  --version |& head -n1 )
    END_VERSIONS
    """
}
