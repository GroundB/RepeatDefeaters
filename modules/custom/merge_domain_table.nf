process MERGE_DOMAIN_TABLE {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple path ( hmmscan_tbl ), path ( hmmscan_pfamtbl )

    output:
    path '*.domtbl',     emit: domain_table
    path "versions.yml", emit: versions

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: hmmscan_tbl.baseName
    """
    # Combine tbl and pfamtbl data to
    # unclassified_domain_table format
    # col field               file-index   join-output-col
    #   1 <seq id>            1.4          1
    #   2 <alignment start>   2.9          2
    #   3 <alignment end>     2.10         3
    #   4 <envelope start>    2.7          4
    #   5 <envelope end>      2.8          5
    #   6 <hmm acc>           "-"
    #   7 <hmm name>          1.2          6
    #   8 <type>              "-"
    #   9 <hmm start>         2.11         7
    #  10 <hmm end>           2.12         8
    #  11 <hmm length>        `calculate`
    #  12 <bit score>         1.10         9
    #  13 <E-value>           1.9          10
    #  14 <significance>      "1"
    #  15 <clan>              "-"

    join -1 1 -2 1 -o1.4,2.9,2.10,2.7,2.8,1.2,2.11,2.12,1.10,1.9 \\
        <( grep -v -e '^[[:space:]]*\$' -e '^#' "$hmmscan_tbl" | \\
            awk '{ print \$1"-"\$8"-"\$9" "\$0 } ' | sort -k1,1 ) \\
        <( grep -v -e '^[[:space:]]*\$' -e '^#' "$hmmscan_pfamtbl" | \\
            awk 'NF > 10 { print \$1"-"\$3"-"\$2" "\$0 } ' | sort -k1,1 ) | \\
        awk '{ \$5=\$5 "\\t-\\t"; \$6=\$6 "\\t-\\t"; \$8=\$8 "\\t-\\t"; \$10=\$10 "\\t1\\t-"; print \$0 }' \\
        > ${prefix}.domtbl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        join: \$( join --version | sed '1 !d;s/[^0-9]*\\(\\([0-9]\\.\\)\\{0,4\\}[0-9][^.]\\).*/\\1/' )
    END_VERSIONS
    """

}
