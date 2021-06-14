params.options = [:]

process DIVSUM {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img'
    } else {
        container 'biocontainers/biocontainers:v1.2.0_cv1'
    }

    input:
    path div

    output:
    path '', emit:

    script:
    """
    DIV=cclaro_1_nanopore_idrenamed.align.devil.divsum
    cat $DIV|sed '1,/Coverage\ for/d' > $DIV.simple
    sed '/Coverage\ for/q' $DIV|sed '/^$/d'|sed '$d'|sed '1,4d'|sed '2d' > $DIV.clean
    head -1 $DIV.clean >> $DIV.clean.nonTE
    grep -wf RhiirA1_1 RhiirA1_1_AssemblyScaffolds.fasta.align.divsum.clean >> RhiirA1_1_AssemblyScaffolds.fasta.align.divsum.clean.nonTE
    grep -vwEf RhiirA1_1 RhiirA1_1_AssemblyScaffolds.fasta.align.divsum.clean > destinationfile
    sed 's/Unknown/Non-TE-Repeat/g' RhiirA1_1_AssemblyScaffolds.fasta.align.divsum.clean.nonTE |sed '1d' >> destinationfile
    """
}
