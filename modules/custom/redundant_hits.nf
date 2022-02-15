process REDUNDANT_HITS {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path blast_tsv

    output:
    path "self_comparison.tsv", emit: tsv
    path "versions.yml"       , emit: versions

    script:
    """
    # Show redundancy in the annotated library
    ## For queries and subjects that are both Unknown
    ## Get total matched bases/query length, which
    ## indicates the coverage of an overlap
    ## Filter coverage using 0.6 as cut-off then dedup
    grep Unknown $blast_tsv | \\
        awk 'NR >= 1 {
            \$8=(\$6)/(\$2)
        }
        \$8 > 0.6 && \$1 != \$3 && !a[\$5\$6]++ {
            print \$0
        }' > self_comparison.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep : \$( grep  --version |& head -n1 )
        awk  : \$( awk  -W version |& head -n1 )
    END_VERSIONS
    """
}