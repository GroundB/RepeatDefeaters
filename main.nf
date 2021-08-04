#! /usr/bin/env nextflow

// Enable DSL2 syntax for Nextflow
nextflow.enable.dsl = 2

// Print parameters to screen before running workflow.
log.info("""
NBIS support 5861

 Annotation of unclassified TEs
 ===================================
""")

// Check a project allocation is given for running on Uppmax clusters.
if(workflow.profile == "uppmax" && !params.project){
    exit 1, "Please provide a SNIC project number ( --project )!\n"
}

// The main workflow
workflow {

    main:
        RENAME_REPEAT_MODELER_SEQUENCES(
            file(params.repeat_modeler_fasta, checkIfExists:true)
            params.species_short_name)
        PFAM_TRANSPOSIBLE_ELEMENT_SEARCH(
            file(params.pfam_a_db, checkIfExists:true)
            file(params.transposon_keywords, checkIfExists:true))
        MAKEBLASTDB(
            file(params.protein_reference, checkIfExists:true))
        BLASTX(
            RENAME_REPEAT_MODELER_SEQUENCES.out.fasta,
            BLAST_MAKEBLASTDB.out.db,
            ['plus','minus'])
        PFAM_SCAN(BLASTX.out.fasta,params.pfam_hmm_db)
        ANNOTATION(PFAM_SCAN.out.pfam_table.collect(),
            PFAM_TRANSPOSIBLE_ELEMENT_SEARCH.out.pfam_accessions,
            params.species_short_name)

}

process RENAME_REPEAT_MODELER_SEQUENCES {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img'
    } else {
        container 'biocontainers/biocontainers:v1.2.0_cv1'
    }

    input:
    path fasta      // Repeat Modeler fasta
    val sci_name    // Short name species identifier

    output:
    path '*.fasta', emit: fasta

    script:
    """
    renameRMDLconsensi.pl $fasta $sciname ${sci_name}.fasta
    """
}

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
    # Search for keywords and IDs
    zgrep -i -e "^#=GF ID" -f $keywords $uniprot_db > pattern_matches.txt
    # Print closest ID above keyword match
    awk '{ if (\$0 ~ /#=GF ID/) { id_line = \$0 } else { print id_line } }' pattern_matches.txt | uniq | cut -c11- > Pfam.TE.accessions
    """
}

process MAKEBLASTDB {

    conda (params.enable_conda ? 'bioconda::blast=2.10.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/blast:2.10.1--pl526he19e7b1_3'
    } else {
        container 'quay.io/biocontainers/blast:2.10.1--pl526he19e7b1_3'
    }

    input:
    path fasta

    output:
    path "blastdb"    , emit: db
    path "*.version"  , emit: version

    script:
    """
    makeblastdb \\
        ${params.modules['makeblastdb'].args} \\
        -in $fasta \\
        -dbtype prot
    mkdir blast_db
    mv ${fasta}* blast_db
    makeblastdb -version | sed -e '/^makeblastdb:/!d; s/^.*makeblastdb: //' > ${software}.version.txt
    """

}

process BLASTX {

    conda (params.enable_conda ? 'bioconda::blast=2.10.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/blast:2.10.1--pl526he19e7b1_3'
    } else {
        container 'quay.io/biocontainers/blast:2.10.1--pl526he19e7b1_3'
    }

    input:
    val query
    path db
    each val strand

    output:
    path "*.blastx.tsv"       , emit: tsv
    path "*.predicted.fasta"  , emit: fasta
    path "*.version"          , emit: version

    script:
    def prefix = query.baseName
    """
    blastx -num_threads ${task.cpus} \\
        -query $query \\
        -db $db \\
        ${params.modules['blastx'].args} \\
        -outfmt "6 qseqid qseq" \\
        -strand $strand \\
        -out ${prefix}.${strand}.blastx.tsv
    # Blast output is a two column table
    # awk: Format two column output to fasta
    #      Label sequence header with strand and
    #      accumulator to form unique sequence headers
    awk '{
        seq = ""
        i = 1
        if (seq == \$1) {
            i++
        } else {
            i = 1
            seq = \$1
        }
        gsub(/[-X*]/,"",\$2)
        print ">"\$1"_${strand}_qseq_"i"\n"\$2
    }' ${prefix}.${strand}.blastx.tsv > ${prefix}.${strand}.predicted.fasta
    blastx -version | sed -e '/^blastx:/!d; s/^.*blastx: //' > blastx.version
    """

}

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
    path '*.pfamtbl', emit: pfam_table
    path '*.version', emit: version

    script:
    def prefix = fasta.baseName
    """
    pfam_scan.pl \\
        -fasta $fasta \\
        -dir $PFAM \\
        -outfile ${prefix}.pfamtbl \\
        ${params.modules['pfam'].args}
    """

}

process ANNOTATION {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img'
    } else {
        container 'biocontainers/biocontainers:v1.2.0_cv1'
    }

    input:
    path pfam_table                // Pfam output (*.pfamtbl)
    path pfam_keyword_accession    // Pfam accessions from keyword search
    val prefix                     // prefix

    output:
    path "", emit: xml

    script:
    """
    for TBL in $pfam_table; do
        grep -i "#unknown" \$TBL | \\
        grep -i -w -f $pfam_keyword_accession | \\
        tee -a ${prefix}.Unclassified_consensus_TEs | \\
        cut -f1 -d"#" | uniq > \${TBL/.pfamtbl/.unclassified_ids}
    done

    ## printf "subsetting the pfamtbl output to unclassified consensus only. \n"
    ## grep -i "#unknown" $RENAME.plus.predicted.pfamtbl > $RENAME.plus.predicted.pfamtbl.unclassified
    ## grep -i "#unknown" $RENAME.minus.predicted.pfamtbl > $RENAME.minus.predicted.pfamtbl.unclassified

    ## printf "Detect unclassified consensi that have TE domains. Plus strand:\n`grep -i -w -f $PFAMTE $RENAME.plus.predicted.pfamtbl.unclassified` \n"
    ## printf "Detect unclassified consensi that have TE domains. Minus strand:\n`grep -i -w -f $PFAMTE $RENAME.minus.predicted.pfamtbl.unclassified` \n"

    ## grep -i -w -f $PFAMTE $RENAME.plus.predicted.pfamtbl.unclassified  | awk '{print $1}' | sed 's/\#.*//g' | uniq > 001.tem
    ## grep -i -w -f $PFAMTE $RENAME.minus.predicted.pfamtbl.unclassified | awk '{print $1}' | sed 's/\#.*//g' | uniq > 002.tem

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

