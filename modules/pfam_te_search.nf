params.options = [:]

process PFAM_TRANSPOSIBLE_ELEMENT_SEARCH {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img'
    } else {
        container 'biocontainers/biocontainers:v1.2.0_cv1'
    }

    input:
    path uniprot_db             // Compressed db
    path keywords

    output:
    path "Pfam.TE.accessions", emit: pfam_accessions

    script:
    """
    zgrep -i -f $keywords $uniprot_db > pattern_matches.txt
    awk '{ if (\$0 ~ /#=GF ID/) { id_line = \$0 } else { print id_line } }' pattern_matches.txt | uniq | cut -c11- > Pfam.TE.accessions
    """
}
