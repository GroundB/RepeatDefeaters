process RENAME_SEQUENCES {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path fasta      // Repeat Modeler fasta
    val sci_name    // Short name species identifier

    output:
    path '*.fasta'      , emit: fasta
    path "versions.yml" , emit: versions

    script:
    """
    renameRMDLconsensi.pl $fasta $sci_name ${sci_name}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl : \$( perl -v |& sed '/This is perl/!d; s/.*\\(v[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/' )
    END_VERSIONS
    """
}
