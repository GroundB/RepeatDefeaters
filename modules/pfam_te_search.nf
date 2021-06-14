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

    grep -n -i "transpos" pfam-a.uniprot.desc | \\
        grep "#=GF CC" | \\
        sed 's/:.*//g' | \\
        sort -n > tranpos.lines
    grep -n -i "Aspartyl protease" Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > Aspartyl.lines
    grep -n -i Asp_protease Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > Asp_protease.lines
    grep -n -i -w gag Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > gag.lines
    grep -n -i RNase_H Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > RNase_H.lines
    grep -n -i virus Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > virus.lines
    grep -n -i viral Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > viral.lines
    grep -n -i integrase Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > integrase.lines
    grep -n -i transcriptase Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > RT.lines
    grep -n -i Helicase Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > RT.lines
    grep -n -i envelope Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > envelope.lines
    grep -n -i capsid Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > capsid.lines
    grep -n -i endonuclease Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > endonuclease.lines
    grep -n -i -w CENP Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > CENP.lines
    grep -n -i Recombin Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > Recombin.lines

    grep -n -i transpos Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> tranpos.lines
    grep -n -i Aspartyl\ protease Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> Aspartyl.lines
    grep -n -i Asp_protease Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> Asp_protease.lines
    grep -n -i -w gag Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> gag.lines
    grep -n -i RNase_H Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> RNase_H.lines
    grep -n -i virus Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> virus.lines
    grep -n -i viral Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> viral.lines
    grep -n -i integrase Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> integrase.lines
    grep -n -i transcriptase Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> RT.lines
    grep -n -i Helicase Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> RT.lines
    grep -n -i envelope Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> envelope.lines
    grep -n -i capsid Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> capsid.lines
    grep -n -i endonuclease Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> endonuclease.lines
    grep -n -i -w CENP Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> CENP.lines
    grep -n -i Recombin Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> Recombin.lines

    cat *.lines | uniq | sort -n > all.line.id

    nl uniprot.des > uniprot.des.nl

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