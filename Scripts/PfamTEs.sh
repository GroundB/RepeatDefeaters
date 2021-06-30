#!/bin/bash
#title           :PfamTEs
#description     :Find PFAM accessions that may associate with TE activities
#author		 :wenbo.xu@botinst.uzh.ch
#date            :20201012
#version         :0.1
#usage		 :./PfamTEs.sh

#some families (e.g. MULE) the transposon info is only present in
#the DES section not in CC section.
#I was planning on using ONLY #=GF CC to filter...
#but decided not to.. ".lines" files will be VERY redundant and slow to run
#If you want to use unique lines
#do grep -n -i transpos Pfam-A.full.uniprot.des|grep \#\=GF\ CC

#get pfam files
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.full.uniprot.gz
gunzip Pfam-A.full.uniprot.gz
grep \# Pfam-A.full.uniprot > Pfam-A.full.uniprot.des
rm Pfam-A.full.uniprot
#key word greping
# This is fetching the line numbers of the match.
# grep -n -i "^#GF CC.*transpos.*$" <des> | = <line_num>:<grep_match> = | sed = cut -f1 -d:
grep -n -i transpos Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > tranpos.lines
grep -n -i Aspartyl\ protease Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > Aspartyl.lines
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
#grep -n -i Replication Pfam-A.full.uniprot.des|grep \#\=GF\ CC|sed 's/:.*//g'|sort -n > replication.lines
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
#grep -n -i Replication Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> replication.lines
grep -n -i Recombin Pfam-A.full.uniprot.des|grep \#\=GF\ DE|sed 's/:.*//g'|sort -n >> Recombin.lines

# This doesn't correctly make a uniquely sorted file.
cat *.lines|uniq|sort -n > all.line.id

# Adds line numbers to description
nl Pfam-A.full.uniprot.des > Pfam-A.full.uniprot.des.nl

# variable s = Find line number of nearest occurrance of string STOCKHOLM
for i in `cat all.line.id`; do
    a=$(expr $i - 200)
    s=$(cat Pfam-A.full.uniprot.des.nl | awk -v vari="$i" -v vara="$a" 'NR >= vara && NR <= vari' | grep STOCKHOLM | tail -1 | sed 's/#.*//g')
    echo $s $i >> ranges
done

# For each stockholm, grep last occurence of number?
# Remove duplicate entries
for i in $(awk '{print $1}' ranges | uniq | sort -n );do
    grep $i ranges | tail -n 1 >> ranges.uniq;
done
sort -k1,2n ranges.uniq > ranges.uniq.sorted

length=$(cat ranges.uniq.sorted|wc -l)
for i in $(seq 1 $length)
do d=$(awk -v row="$i" 'FNR == row {print $1}' ranges.uniq.sorted)
e=$(awk -v row="$i" 'FNR == row {print $2}' ranges.uniq.sorted)
cat Pfam-A.full.uniprot.des|awk -v vard="$d" -v vare="$e" 'NR >= vard && NR <= vare' >> final
done

grep GF\ ID final|uniq|sed 's/.*\ //g' > ../Resources/Pfam.TE.accessions_Release32
cat ../Resources/Pfam.TE.accessions_Release32
rm ranges.uniq.sorted ranges all.line.id Pfam-A.full.uniprot final Pfam-A.full.uniprot.des ranges.uniq *.lines Pfam-A.full.uniprot.des.nl
