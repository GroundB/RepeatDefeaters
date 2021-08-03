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
        BLAST_BLASTX(
            RENAME_REPEAT_MODELER_SEQUENCES.out.fasta,
            BLAST_MAKEBLASTDB.out.db)
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
