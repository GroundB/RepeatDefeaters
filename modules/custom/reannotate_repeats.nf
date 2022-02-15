process REANNOTATE_REPEATS {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path repeat_library             // Renamed repeat consensus library
    path unclassified_domain_table  // Pfam A unclassifed output
    path hmmscan_domain_table       // hmmscan output (*.tbl + *.pfamtbl = *.domtbl )

    output:
    path "*_reannotated.fasta", emit: fasta
    path "domain_table.tsv"   , emit: domain_table
    path "versions.yml"       , emit: versions

    script:
    def prefix = repeat_library.baseName
    """
    #! /usr/bin/env python

    import os
    import re
    import platform

    # Checklist determines if an annotation uses multiple hits or superfamily name
    superfamily = [
        "Academ",
        "CACTA",
        "DIRS",
        "Ginger",
        "Harbinger",
        "Helitron",
        "Kolobok",
        "MULE",
        "Mariner",
        "Merlin",
        "Novosib",
        "P",
        "Piggybac",
        "Sola1",
        "Sola2",
        "Sola3",
        "Transib",
        "Zator"
    ]

    # domain_table format
    # col field
    #   1 <seq id>
    #   2 <alignment start>
    #   3 <alignment end>
    #   4 <envelope start>
    #   5 <envelope end>
    #   6 <hmm acc>
    #   7 <hmm name>
    #   8 <type>
    #   9 <hmm start>
    #  10 <hmm end>
    #  11 <hmm length>
    #  12 <bit score>
    #  13 <E-value>
    #  14 <significance>
    #  15 <clan>

    header_lines = {}

    # Concatenate the domain table and cut out relevant columns
    with os.popen(
        'cat $hmmscan_domain_table $unclassified_domain_table | tee domain_table.tsv | awk \\'{ print \$1 "\\t" \$7 "\\t" \$12 }\\''
    ) as domain_table:

        # Parse domain table to form a lookup table.
        for line in domain_table:
            cols = line.split()
            header = re.split("_(plus|minus)", cols[0])[
                0
            ]  # e.g. cclaro4-779#Unknown_minus_qseq_1
            hit = cols[1].partition(".")[0]  # PARA-PROT.uniq.fasta -> PARA-PROT
            score = cols[2]  # bitscore
            print("Lookup: " + header + " : " + hit + " : " + score)

            # make entries for both concatenated names and best hit for each header
            if header in header_lines:
                header_lines[header]["concat"].add(hit)
                if score < header_lines[header]["best_hit_score"]:
                    header_lines[header]["best_hit_score"] = score
                    header_lines[header]["best_hit"] = hit
            else:
                header_lines.setdefault(header, {})["concat"] = {hit}
                header_lines[header]["best_hit_score"] = score
                header_lines[header]["best_hit"] = hit

    # Iterate over repeat library and update header information
    with open("$repeat_library", "r") as repeat_lib:
        with open("${prefix}_reannotated.fasta", "w") as renamed_repeats:
            for line in repeat_lib:

                # Get lines matching header e.g., cclaro4-1265#Unknown(/XXX)
                line_match = re.match("^>(.+#Unknown)(/[A-Z]{3})?\$", line)
                if line_match:
                    print("Old header: " + line_match.group(0))

                    # Update header and print to file
                    # Scan keys for matching substring
                    key_found = None
                    for key in header_lines:
                        header_key = re.search(line_match.group(1), key)
                        if header_key:
                            key_found = key  # matches updated to be unique key

                            # if best hit is superfamily use that, otherwise use concat
                            if header_lines[key]["best_hit"] in superfamily:
                                # Write superfamily name
                                header = (
                                    ">"
                                    + key.partition("#Unknown")[0]
                                    + "#"
                                    + header_lines[key]["best_hit"]
                                    + str(line_match.group(2) or "")
                                )
                                renamed_repeats.write(header + "\\n")
                                print("New header: " + header)
                            else:
                                # Write concatenation of names
                                header = (
                                    ">"
                                    + key.partition("#Unknown")[0]
                                    + "#"
                                    + "-".join(
                                        [str(hmm) for hmm in header_lines[key]["concat"]]
                                    )
                                    + str(line_match.group(2) or "")
                                )
                                renamed_repeats.write(header + "\\n")
                                print("New header: " + header)

                    # If header is unmatched, write out header again.
                    if not key_found:
                        renamed_repeats.write(line)
                        print("No match")

                else:
                    # print to file
                    renamed_repeats.write(line)

    with open("versions.yml","w") as versions:
        versions.write('"${task.process}":\\n    python: ' + platform.python_version() + '\\n')
    """
}
