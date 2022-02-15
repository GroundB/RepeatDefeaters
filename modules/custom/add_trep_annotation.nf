process ADD_TREP_ANNOTATION {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path repeat_library            // Reannotated repeat consensus library
    path trep_blast_hits           // TREP blast results ( *.blastn.tsv )

    output:
    path "*.trep.fasta", emit: fasta
    path "versions.yml", emit: versions

    script:
    def prefix = repeat_library.baseName
    """
    # Expected blast outfmt 6
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    # Find queries ending in #Unknown and filter for smallest evalue

    # Find hits with smallest e-value
    awk 'BEGIN {
            prev_seq = ""    # sequence name
            prev_sfam = ""   # transposon superfamily
            prev_evalue = 0
        }
        \$1 ~ /#Unknown\$/ {
            # If current sequence is the same as previous sequence name
            # and evalue is smaller, update evalue
            if ( \$1 == prev_seq && \$11 < prev_evalue ) {
                prev_evalue = \$11
            }
            # If current sequence has a different name to previous, then
            # print replacement info
            if ( \$1 != prev_seq && prev_seq != "") {
                print prev_seq "\\t" prev_sfam
            }
            prev_seq = \$1
            prev_sfam = substr(\$2,1,3)
        }
        # Process last sequence.
        END {
            if ( prev_seq != "" ) {
                print prev_seq "\\t" prev_sfam
            }
        }' $trep_blast_hits > best_hits.tsv

    if [ -s best_hits.tsv ]; then
        # Read in annotations to memory, and use dictionary
        # to update fasta headers with transposon superfamily
        awk '
            # Load best_hits.tsv into dictionary
            FNR == NR {
                annotation[\$1] = \$2
                next
            }
            # Update annotation
            # If a line begins with >
            # and an annotation to add exists.
            /^>/ {
                seq_head = substr(\$1,2)
                if ( seq_head in annotation ) {
                    \$1 = \$1 "/" annotation[seq_head]
                }
            }
            # print each line
            1
            ' best_hits.tsv $repeat_library > ${prefix}.trep.fasta
    else
        cp $repeat_library ${prefix}.trep.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk  : \$( awk  -W version |& head -n1 )
    END_VERSIONS
    """

}