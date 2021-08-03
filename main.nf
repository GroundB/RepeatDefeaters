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

include { BLAST_MAKEBLASTDB                     } from './modules/blast_makeblastdb', addParams(options:parmas.modules['blast_makeblastdb'])
include { BLAST_BLASTX                          } from './modules/blast_blastx'     , addParams(options:params.modules['blast_positive_strand'])
include { FILTER_BLAST_XML                      } from './modules/blast_xml_filter' , addParams(options:[:])
include { PFAM_SCAN                             } from './modules/pfam_scan'        , addParams(options:params.modules['pfam_scan'])
include { FILTER_PFAM                           } from './modules/pfam_filter'      , addParams(options:[:])
include { ANNOTATION                            } from './modules/annotation'       , addParams(options:[:])
include { PFAM_TRANSPOSIBLE_ELEMENT_SEARCH      } from './modules/pfam_te_search'   , addParams(options:[:])

// The main workflow
workflow {

    main:
        RENAME_REPEAT_MODELER_SEQUENCES(
            file(params.repeat_modeler_fasta, checkIfExists:true)
            params.species_short_name)
        PFAM_TRANSPOSIBLE_ELEMENT_SEARCH(
            file(params.pfam_a_db, checkIfExists:true)
            file(params.transposon_keywords, checkIfExists:true))
        BLAST_MAKEBLASTDB(
            file(params.protein_reference, checkIfExists:true))
        BLASTX(
            RENAME_REPEAT_MODELER_SEQUENCES.out.fasta,
            BLAST_MAKEBLASTDB.out.db,
            ['plus','minus'])
        PFAM_SCAN(BLAST_BLASTX.out.fasta)
        ANNOTATION(PFAM_SCAN.out.pfam_table.collect())

}

process RENAME_REPEAT_MODELER_SEQUENCES {

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
        $options.args \\
        -outfmt "6 qseqid qseq" \\
        -strand $strand \\
        -out ${prefix}.blastx.tsv
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
    }' ${prefix}.blastx.tsv > ${prefix}.predicted.fasta
    blastx -version | sed -e '/^blastx:/!d; s/^.*blastx: //' > blastx.version
    """

}
