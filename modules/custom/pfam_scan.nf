process PFAM_SCAN {

    conda (params.enable_conda ? 'bioconda::pfam_scan==1.6' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pfam_scan:1.6--hdfd78af_4' :
        'quay.io/biocontainers/pfam_scan:1.6--hdfd78af_4' }"

    input:
    path fasta
    path hmm_db

    output:
    path '*.pfamtbl'    , emit: pfam_table
    path "versions.yml" , emit: versions

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: fasta.baseName
    """
    # Stage local copy of HMM DB
    mkdir -p HMM_DB
    printf "%s\\n" $hmm_db | \\
        xargs -P $task.cpus -I {} \\
        sh -c 'gzip -cdf {} > HMM_DB/\$( basename {} .gz )'
    find HMM_DB -name "*.hmm" -exec hmmpress {} \\;

    pfam_scan.pl \\
        -fasta $fasta \\
        -dir HMM_DB \\
        -cpu ${task.cpus} \\
        -outfile ${prefix}.pfamtbl \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pfam_scan: 1.6
    END_VERSIONS
    """

}
