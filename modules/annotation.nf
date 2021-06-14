params.options = [:]

process ANNOTATION {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img'
    } else {
        container 'biocontainers/biocontainers:v1.2.0_cv1'
    }

    input:
    path pfam_table // list of pfam_tables
    val prefix // prefix

    output:
    path "", emit: xml

    script:
    """
    printf "subsetting the pfamtbl output to unclassified consensus only. \n"
    grep -i "#unknown" $RENAME.plus.predicted.pfamtbl > $RENAME.plus.predicted.pfamtbl.unclassified
    grep -i "#unknown" $RENAME.minus.predicted.pfamtbl > $RENAME.minus.predicted.pfamtbl.unclassified

    printf "Detect unclassified consensi that have TE domains. Plus strand:\n`grep -i -w -f $PFAMTE $RENAME.plus.predicted.pfamtbl.unclassified` \n"
    printf "Detect unclassified consensi that have TE domains. Minus strand:\n`grep -i -w -f $PFAMTE $RENAME.minus.predicted.pfamtbl.unclassified` \n"

    grep -i -w -f $PFAMTE $RENAME.plus.predicted.pfamtbl.unclassified  | awk '{print $1}' | sed 's/\#.*//g' | uniq > 001.tem
    grep -i -w -f $PFAMTE $RENAME.minus.predicted.pfamtbl.unclassified | awk '{print $1}' | sed 's/\#.*//g' | uniq > 002.tem

    cat 001.tem 002.tem | uniq > $RENAME.Unclassified_consensus_TEs.ids #can be used to annotate unknown consensus to TEs

    grep -i -w -f $PFAMTE $RENAME.plus.predicted.pfamtbl.unclassified > $RENAME.Unclassified_consensus_TEs
    grep -i -w -f $PFAMTE $RENAME.minus.predicted.pfamtbl.unclassified >> $RENAME.Unclassified_consensus_TEs
    rm *.tem

    grep -vf  $RENAME.Unclassified_consensus_TEs.ids $RENAME.plus.predicted.pfamtbl.unclassified > $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged
    grep -vf  $RENAME.Unclassified_consensus_TEs.ids $RENAME.minus.predicted.pfamtbl.unclassified > $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged
    printf "Gathering unclassified consensus sequences that have domains detected in both strands.\n"
    printf "These consensus sequences are tricky to annotate.\nPut aside for anyone that are interested. Check file consensus.proteins.both.strand"

    grep -i "#unknown" $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged | \\
        awk '{print $1}' | \\
        sed 's/\#.*//g' | \\
        uniq > 003.tem
    grep -i "#unknown" $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged | \\
        awk '{print $1}' | \\
        sed 's/\#.*//g' | \\
        uniq > 004.tem

    # TODO:: This appears to be a inner join.
    printf "These consensus sequences are:\n`grep -f 003.tem 004.tem` \n"
    grep -f 003.tem 004.tem | uniq > $RENAME.consensus.proteins.both.strand
    rm  *.tem

    # TODO:: This appears to be set difference.
    printf "Now short matches that are shorter than 100AA are removed"
    grep -vf $RENAME.consensus.proteins.both.strand $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged > $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged.mono
    awk '{ if ($11 >= 100) { print } }' $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged.mono > c$RENAME.plus.predicted.pfamtbl.unclassified.TEpurged.mono.final
    grep -vf $RENAME.consensus.proteins.both.strand $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged > $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged.mono
    awk '{ if ($11 >= 100) { print } }' $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged.mono > $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged.mono.final

    # TODO:: Why is there a copy command in the print statement?
    printf "Annotation initiated.\nBacking up the original library file.\n`cp $RRLQUERY $RRLQUERY.bk`"

    # TODO:: unneeded cat
    cat $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged.mono.final | \\
        awk '{print $1}' | \\
        sed 's/\#.*//g' | \\
        uniq > 005.tem
    cat $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged.mono.final | \\
        awk '{print $1}' | \\
        sed 's/\#.*//g' | \\
        uniq > 006.tem
    printf "Following repeat models were renamed based on plus strand evidence:\n"
    for i in `cat 005.tem`; do
        echo \$i >> $RENAME.renamed;
        grep "\$i" $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged.mono.final | \\
            awk '{print \$7}' >> $i.plus.annotate.tem;
    done
    for i in *.plus.annotate.tem; do
        namehash=$( cat \$i | sort | uniq | paste -s -d '_');
        oldname=$( grep "\${i%.plus.annotate.tem}#" $RRLQUERY | sed 's/>//g' );
        sed -i "s|$oldname|${oldname%Unknown}$namehash|g" $RRLQUERY;
    done
    printf "Following repeat models were renamed based on minus strand evidence:\n"
    for i in `cat 006.tem`;do
        echo \$i >> $RENAME.renamed;
        grep "\$i" $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged.mono.final | \\
            awk '{print $7}' >> $i.minus.annotate.tem;
    done
    for i in *.minus.annotate.tem; do
        namehash=$( cat $i | sort | uniq | paste -s -d '_' );
        oldname=$( grep "\${i%.minus.annotate.tem}#" $RRLQUERY | sed 's/>//g' );
        sed -i "s|$oldname|${oldname%Unknown}$namehash|g" $RRLQUERY;
    done
    """
}
