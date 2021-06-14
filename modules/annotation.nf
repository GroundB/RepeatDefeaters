params.options = [:]

process ANNOTATION {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img'
    } else {
        container 'biocontainers/biocontainers:v1.2.0_cv1'
    }

    input:
    path blast_xml // list of files?
    val prefix // prefix

    output:
    path "", emit: xml

    script:
    """
    printf "Annotation initiated.\nBacking up the original library file.\n`cp $RRLQUERY $RRLQUERY.bk`"
    cat $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged.mono.final|awk '{print $1}'|sed 's/\#.*//g'|uniq > 005.tem
    cat $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged.mono.final|awk '{print $1}'|sed 's/\#.*//g'|uniq > 006.tem
    printf "Following repeat models were renamed based on plus strand evidence:\n"
    for i in `cat 005.tem`;do echo $i >> $RENAME.renamed; grep $i $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged.mono.final|awk '{print $7}' >> $i.plus.annotate.tem; done
    for i in *.plus.annotate.tem; do namehash=$(cat $i|sort|uniq|paste -s -d '_'); oldname=$(grep ${i%.plus.annotate.tem}\# $RRLQUERY |sed 's/>//g'); sed -i "s|$oldname|${oldname%Unknown}$namehash|g" $RRLQUERY; done
    printf "Following repeat models were renamed based on minus strand evidence:\n"
    for i in `cat 006.tem`;do echo $i >> $RENAME.renamed; grep $i $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged.mono.final|awk '{print $7}' >> $i.minus.annotate.tem; done
    for i in *.minus.annotate.tem; do namehash=$(cat $i|sort|uniq|paste -s -d '_'); oldname=$(grep ${i%.minus.annotate.tem}\# $RRLQUERY |sed 's/>//g'); sed -i "s|$oldname|${oldname%Unknown}$namehash|g" $RRLQUERY; done
    """
}
