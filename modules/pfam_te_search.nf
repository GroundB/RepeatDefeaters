params.options = [:]

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
    zgrep "#" $uniprot_db > pfam-a.uniprot.desc
    grep -n -i -f $keywords pfam-a.uniprot.desc | cut -f1 -d: | sort -n > all.lines.id

    nl pfam-a.uniprot.desc > pfam-a.uniprot.desc.nl

    for i in `cat all.line.id`; do
        a=\$( expr \$i - 200 )
        s=\$( cat Pfam-A.full.uniprot.des.nl | awk -v vari="\$i" -v vara="\$a" 'NR >= vara && NR <= vari' | grep STOCKHOLM | tail -1 | sed 's/#.*//g' )
        echo \$s \$i >> ranges
    done

    for i in \$( awk '{print \$1}' ranges | uniq |sort -n ); do
        grep "\$i" ranges | tail -n 1 >> ranges.uniq
    done
    sort -k1,2n ranges.uniq > ranges.uniq.sorted

    length=\$( cat ranges.uniq.sorted | wc -l )
    for i in \$( seq 1 \$length ); do
        d=\$( awk -v row="\$i" 'FNR == row {print \$1}' ranges.uniq.sorted )
        e=\$( awk -v row="\$i" 'FNR == row {print \$2}' ranges.uniq.sorted )
        cat Pfam-A.full.uniprot.des | awk -v vard="\$d" -v vare="\$e" 'NR >= vard && NR <= vare' >> final
    done

    grep "GF ID" final | uniq | sed 's/.*\ //g' > ../Resources/Pfam.TE.accessions_Release32
    # cat ../Resources/Pfam.TE.accessions_Release32
    # rm ranges.uniq.sorted ranges all.line.id Pfam-A.full.uniprot final Pfam-A.full.uniprot.des ranges.uniq *.lines Pfam-A.full.uniprot.des.nl
    """
}
