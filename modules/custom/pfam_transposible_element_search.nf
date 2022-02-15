process PFAM_TRANSPOSIBLE_ELEMENT_SEARCH {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path uniprot_db             // Compressed db
    path keywords
    path blacklist

    output:
    path "Pfam.Proteins_wTE_Domains.seqid", emit: te_domain_proteins
    path "versions.yml"                   , emit: versions

    script:
    """
    # Search for keywords and IDs
    zgrep -i -e "^#=GF ID" -f $keywords $uniprot_db > pattern_matches.txt
    # Remove blacklisted keywords
    grep -i -v -f $blacklist pattern_matches.txt > pattern_matches_revised.txt
    # Print closest ID above keyword match
    awk '{
        if (\$0 ~ /#=GF ID/) {
            id_line = \$0
        } else {
            print id_line
        }
    }' pattern_matches_revised.txt | \\
        uniq | cut -c11- > Pfam.Proteins_wTE_Domains.seqid

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk  : \$( awk  -W version |& head -n1 )
        cut  : \$( cut   --version |& head -n1 )
        uniq : \$( uniq  --version |& head -n1 )
        zgrep: \$( zgrep --version |& head -n1 )
    END_VERSIONS
    """
}
