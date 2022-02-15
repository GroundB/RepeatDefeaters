process HMMSCAN {

    conda (params.enable_conda ? 'bioconda::pfam_scan==1.6' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pfam_scan:1.6--hdfd78af_4' :
        'quay.io/biocontainers/pfam_scan:1.6--hdfd78af_4' }"

    input:
    path fasta
    path hmm_db

    output:
    path '*.hmmscan.out',                   emit: hmmscan_out
    tuple path('*.tbl'), path('*.pfamtbl'), emit: hmmscan_tables
    path "versions.yml",                    emit: versions

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: fasta.baseName
    """
    # Stage local copy of HMM DB
    mkdir -p HMM_DB
    gzip -cdf $hmm_db > HMM_DB/\$( basename $hmm_db .gz )
    find HMM_DB -name "*.hmm" -exec hmmpress {} \\;

    hmmscan $args \\
        --cpu $task.cpus \\
        --pfamtblout ${prefix}.pfamtbl \\
        --tblout ${prefix}.tbl \\
        -o ${prefix}.hmmscan.out \\
        HMM_DB/*.hmm \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmscan: \$( hmmscan -h | sed '2 !d;s/[^0-9]*\\(\\([0-9]\\.\\)\\{0,4\\}[0-9][^.]\\).*/\\1/' )
    END_VERSIONS
    """

}
