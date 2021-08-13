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
            file(params.repeat_modeler_fasta, checkIfExists:true),
            params.species_short_name)
        PFAM_TRANSPOSIBLE_ELEMENT_SEARCH(
            file(params.pfam_a_db, checkIfExists:true),
            file(params.transposon_keywords, checkIfExists:true))
        MAKEBLASTDB(
            file(params.protein_reference, checkIfExists:true))
        BLASTX(
            RENAME_REPEAT_MODELER_SEQUENCES.out.fasta,
            BLAST_MAKEBLASTDB.out.db,
            ['plus','minus'])
        PFAM_SCAN(BLASTX.out.fasta,
            file(params.pfam_hmm_db, checkIfExists:true))
        ANNOTATION(RENAME_REPEAT_MODELER_SEQUENCES.out.fasta,
            PFAM_SCAN.out.pfam_table.collect(),
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
    path repeat_library            // Renamed repeat consensus library
    path pfam_table                // Pfam output (*.pfamtbl)
    path pfam_keyword_accession    // Pfam accessions from keyword search
    val prefix                     // prefix

    output:
    path "${prefix}.renamed.fasta"              , emit: fasta
    path "${prefix}.Unclassified_consensus_TEs" , emit: unclassified_with_te_domains
    path "${prefix}.consensus.both.strand"      , emit: unclassified_with_non_te_domains_both_strands

    script:
    """
    # Find unclassified consensus with TE domains
    for TBL in $pfam_table; do
        # grep #1    : Find unclassified consensus
        # grep #2    : which have TE domains
        # cut + uniq : and extract their id's
        grep -i "#unknown" "\$TBL" | \\
            tee "\${TBL}.unclassified" | \\
            grep -i -w -f $pfam_keyword_accession | \\
            tee -a ${prefix}.Unclassified_consensus_TEs | \\
            cut -f1 -d"#" | uniq > "\${TBL/.pfamtbl/.unclassified_ids}"
    done

    # Concatenate ids of consensus with TE domains from both strands
    cat *.unclassified_ids | uniq > ${prefix}.Unclassified_consensus_TEs.ids

    # Find unclassified consensus without TE domains
    for UNCLASSIFIED in *.unclassified; do
        # grep       : Remove consensus which have TE domains
        # cut + uniq : and extract their id's
        grep -v -f ${prefix}.Unclassified_consensus_TEs.ids "\$UNCLASSIFIED" | \\
            tee "\${UNCLASSIFIED}.TEpurged" | \\
            cut -f1 -d'#' | uniq > "\${UNCLASSIFIED}.TEpurged.ids"
    done

    # Use shell expansion to expand plus and minus strand files for unsorted inner join
    grep -f *.TEpurged.ids > ${prefix}.consensus.both.strand
    # ${prefix}.consensus.both.strand : Unclassified consensus sequences that have
    # non-TE domains detected in both strands.
    # These are tricky to annotate.

    # In consensus without TE domains, remove consensus with non-TE domains on both strands
    # (leaving consensus with non-TE domains on a single-strand)
    for TEPURGED in *.TEpurged; do
        # grep       : Remove consensus with non-TE domains on both strands
        # awk        : then remove consensus shorter than 100 amino acids
        # cut + uniq : and extract their id's
        grep -v -f ${prefix}.consensus.both.strand "\$TEPURGED" | \\
            awk '\$11 >= 100' | tee "\$TEPURGED.mono" | \\
            cut -f1 -d'#' | uniq > "\$TEPURGED.mono.ids"
    done

    # Make a copy of repeat library to be modified.
    cp $repeat_library ${prefix}.renamed.fasta

    # Rename repeat model based on strand evidence.
    for CONSENSUS in *.mono.ids; do
        # while       : for each consensus id
        # echo        : record id as renamed
        # NAMEHASH    : create a name suffix from the pfam domain table
        # OLDNAME     : find old name from repeat consensus library
        # sed         : replace Unknown with NAMEHASH
        while read -r SEQID; do
            echo "\$SEQID" >> ${prefix}.renamed
            NAMEHASH=\$( grep "\$SEQID" "\${CONSENSUS/.ids/}" | cut -f7 | sort | uniq | paste -s -d '_' )
            OLDNAME=\$( grep "\$SEQID" $repeat_library | cut -c2- )
            sed -i "s|\$OLDNAME|\${OLDNAME%Unknown}\$NAMEHASH|g" ${prefix}.renamed.fasta
        done < "\$CONSENSUS"
    done
    """
}

