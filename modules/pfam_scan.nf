params.options = [:]

process PFAM_SCAN {

    conda (params.enable_conda ? 'bioconda::pfam_scan==1.6' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/pfam_scan:1.6--hdfd78af_4'
    } else {
        container 'quay.io/biocontainers/pfam_scan:1.6--hdfd78af_4'
    }

    input:
    path fasta
    path db

    output:
    path '', emit: pfam_table
    path '*.version', emit: version

    script:
    prefix = ''
    """
    pfam_scan.pl \\
        -fasta $fasta \\
        -dir $PFAM \\
        -outfile ${prefix}.plus.predicted.pfamtbl \\
        $options.args
    """

}
