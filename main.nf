#! /usr/bin/env nextflow

// Enable DSL2 syntax for Nextflow
nextflow.enable.dsl = 2

log.info("""
NBIS support 5861

 Annotation of unclassified TEs
 ===================================
""")

// Check a project allocation is given for running on Uppmax clusters.
if(workflow.profile == "uppmax" && !params.project){
    exit 1, "Please provide a SNIC project number ( --project )!\n"
}

// The primary analysis workflow
workflow {

    main:
        // Step 1: Change Fasta headers
        RENAME_REPEAT_MODELER_SEQUENCES(
            file(
                params.repeat_modeler_fasta,
                checkIfExists: true
            ),
            params.species_short_name
        )
        ch_versions = RENAME_REPEAT_MODELER_SEQUENCES.out.versions

        // Step 2: Get ID's of PFAM Proteins with TE domains
        protein_te_domain_list = Channel.empty()
        if ( params.pfam_proteins_with_te_domain_list ){
            protein_te_domain_list = file(
                params.pfam_proteins_with_te_domain_list,
                checkIfExists: true
            )
        } else {
            PFAM_TRANSPOSIBLE_ELEMENT_SEARCH(
                file(
                    params.pfam_a_db,
                    checkIfExists: true
                ),
                file(
                    params.transposon_keywords,
                    checkIfExists: true
                ),
                file (
                    params.transposon_whitelist,
                    checkIfExists: true
                )
            )
            protein_te_domain_list = PFAM_TRANSPOSIBLE_ELEMENT_SEARCH.out.te_domain_proteins
            ch_versions.mix(
                PFAM_TRANSPOSIBLE_ELEMENT_SEARCH.out.versions
            )
        }

        // Step 3: Strand specific Blast search of Repeats against
        // a protein reference database
        MAKEBLASTDB(
            file(
                params.protein_reference,
                checkIfExists: true
            )
        )
        BLASTX(
            RENAME_REPEAT_MODELER_SEQUENCES.out.fasta,
            MAKEBLASTDB.out.db,
            ['plus','minus']
        )

        // Step 4: Scan Repeats for PFAM domains
        PFAM_SCAN(
            BLASTX.out.fasta,
            Channel.fromPath(
                [params.pfam_hmm_db, params.pfam_hmm_dat],
                checkIfExists: true
            ).collect()
        )

        // Step 5: Reannoate Fasta headers to emphasize
        // repeats with single-stranded non-TE domains
        ANNOTATION(
            RENAME_REPEAT_MODELER_SEQUENCES.out.fasta,
            PFAM_SCAN.out.pfam_table.collect(),
            protein_te_domain_list,
            params.species_short_name
        )

        // Report: Record software versions
        ch_versions.mix(
            MAKEBLASTDB.out.versions,
            BLASTX.out.versions.first(),
            PFAM_SCAN.out.versions.first(),
            ANNOTATION.out.versions
        ).collectFile(
            name: "software_versions.yml",
            cache: false,
            sort: true,
            storeDir: "${params.results}/pipeline_info"
        )

}

process RENAME_REPEAT_MODELER_SEQUENCES {

    publishDir "${params.results}/01_Renamed_Repeat_modeler_sequences", mode: params.publish_mode

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
    path '*.fasta'      , emit: fasta
    path "versions.yml" , emit: versions

    script:
    """
    renameRMDLconsensi.pl $fasta $sci_name ${sci_name}.fasta

    cat <<-END_VERSIONS > versions.yml
    RENAME_REPEAT_MODELER_SEQUENCES:
        perl : \$( perl -v |& sed '/This is perl/!d; s/.*\\(v[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/' )
    END_VERSIONS
    """
}

process PFAM_TRANSPOSIBLE_ELEMENT_SEARCH {

    publishDir "${params.results}/02_Pfam_TE_IDs", mode: params.publish_mode

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img'
    } else {
        container 'biocontainers/biocontainers:v1.2.0_cv1'
    }

    input:
    path uniprot_db             // Compressed db
    path keywords
    path whitelist

    output:
    path "Pfam.Proteins_wTE_Domains.seqid", emit: te_domain_proteins
    path "versions.yml"                   , emit: versions

    script:
    """
    # Search for keywords and IDs
    zgrep -i -e "^#=GF ID" -f $keywords $uniprot_db > pattern_matches.txt
    # Remove whitelisted keywords
    grep -i -f $whitelist pattern_matches.txt > pattern_matches_revised.txt
    # Print closest ID above keyword match
    awk '{
        if (\$0 ~ /#=GF ID/) {
            id_line = \$0
        } else {
            print id_line
        }
    }' pattern_matches_revised.txt | \\
        uniq | cut -c11- > Pfam.Proteins_wTE_Domains.seqid

    cat <<-END_VERSIONS > versions.yml
    PFAM_TRANSPOSIBLE_ELEMENT_SEARCH:
        awk  : \$( awk  -W version |& head -n1 )
        cut  : \$( cut   --version |& head -n1 )
        uniq : \$( uniq  --version |& head -n1 )
        zgrep: \$( zgrep --version |& head -n1 )
    END_VERSIONS
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
    path "blast_db"      , emit: db
    path "versions.yml"  , emit: versions

    script:
    """
    makeblastdb \\
        ${params.modules['makeblastdb'].args} \\
        -in $fasta \\
        -dbtype prot
    mkdir blast_db
    mv ${fasta}* blast_db

    cat <<-END_VERSIONS > versions.yml
    MAKEBLASTDB:
        makeblastdb: \$(makeblastdb -version | sed -e '/^makeblastdb:/!d; s/^.*makeblastdb: //' )
    END_VERSIONS
    """

}

process BLASTX {

    publishDir "${params.results}/03_Blastx", mode: params.publish_mode

    conda (params.enable_conda ? 'bioconda::blast=2.10.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/blast:2.10.1--pl526he19e7b1_3'
    } else {
        container 'quay.io/biocontainers/blast:2.10.1--pl526he19e7b1_3'
    }

    input:
    val query
    path db
    each strand

    output:
    path "*.blastx.tsv"       , emit: tsv
    path "*.predicted.fasta"  , emit: fasta
    path "versions.yml"       , emit: versions

    script:
    def prefix = query.baseName
    """
    BLASTDB=\$( find -L ./ -name "*.pdb" | sed 's/.pdb//' )
    blastx -num_threads ${task.cpus} \\
        -query $query \\
        -db \$BLASTDB \\
        ${params.modules['blastx'].args} \\
        -outfmt "6 qseqid qseq" \\
        -strand $strand \\
        -out ${prefix}.${strand}.blastx.tsv
    # Blast output is a two column table
    # awk: Format two column output to fasta
    #      Label sequence header with strand and
    #      accumulator to form unique sequence headers
    awk 'BEGIN {
        seq = ""
        i = 1
    }
    {
        if (seq == \$1) {
            i++
        } else {
            i = 1
            seq = \$1
        }
        gsub(/[-X*]/,"",\$2)
        print ">"\$1"_${strand}_qseq_"i"\\n"\$2
    }' ${prefix}.${strand}.blastx.tsv > ${prefix}.${strand}.predicted.fasta

    cat <<-END_VERSIONS > versions.yml
    BLASTX:
        blastx: \$( blastx -version | sed -e '/^blastx:/!d; s/^.*blastx: //' )
        awk   : \$( awk |& head -n1 | sed 's/(.*//' )
    END_VERSIONS
    """

}

process PFAM_SCAN {

    publishDir "${params.results}/04_Pfam_scan", mode: params.publish_mode

    conda (params.enable_conda ? 'bioconda::pfam_scan==1.6' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/pfam_scan:1.6--hdfd78af_4'
    } else {
        container 'quay.io/biocontainers/pfam_scan:1.6--hdfd78af_4'
    }

    input:
    path fasta
    path hmm_db

    output:
    path '*.pfamtbl'    , emit: pfam_table
    path "versions.yml" , emit: versions

    script:
    def prefix = fasta.baseName
    """
    # Stage local copy of Pfam HMM DB
    mkdir -p HMM_DB
    for FILE in $hmm_db; do
        PREFIX=\$( basename "\$FILE" .gz )
        zcat "\$FILE" > "HMM_DB/\${PREFIX}"
    done
    for HMMDB in HMM_DB/*.hmm; do
        hmmpress \$HMMDB
    done

    pfam_scan.pl \\
        -fasta $fasta \\
        -dir HMM_DB \\
        -cpu ${task.cpus} \\
        -outfile ${prefix}.pfamtbl \\
        ${params.modules['pfam'].args}

    cat <<-END_VERSIONS > versions.yml
    PFAM_SCAN:
        pfam_scan: 1.6
    END_VERSIONS
    """

}

process ANNOTATION {

    publishDir "${params.results}/05_Reannotated_Repeat_modeler_sequences", mode: params.publish_mode

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
    path "versions.yml"                         , emit: versions

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
            NAMEHASH=\$( grep "\${SEQID}#" "\${CONSENSUS/.ids/}" | \\
                tr -s " " "\t" | cut -f7 | \\
                sort | uniq | \\
                paste -s -d '-' )
            OLDNAME=\$( grep "\${SEQID}#" $repeat_library | cut -c2- )
            sed -i "s|\$OLDNAME|\${OLDNAME%Unknown}\$NAMEHASH|g" ${prefix}.renamed.fasta
        done < "\$CONSENSUS"
    done

    cat <<-END_VERSIONS > versions.yml
    ANNOTATION:
        awk  : \$( awk  -W version |& head -n1 )
        cat  : \$( cat   --version |& head -n1 )
        cut  : \$( cut   --version |& head -n1 )
        grep : \$( grep  --version |& head -n1 )
        paste: \$( paste --version |& head -n1 )
        sed  : \$( sed   --version |& head -n1 )
        sort : \$( sort  --version |& head -n1 )
        tee  : \$( tee   --version |& head -n1 )
        uniq : \$( uniq  --version |& head -n1 )
    END_VERSIONS
    """
}

