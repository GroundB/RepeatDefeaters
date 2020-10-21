#!/bin/bash
#title           :Detection and annotation of non-TE-repeats
#description     :filter out host proteins from unclassified repeat consensus
#author		 :wenbo.xu@botinst.uzh.ch
#date            :20201011
#version         :0.3
#usage		 :./hostfinder.sh
#==============================================================================
PROREF=~/conc_sprot.fa
RRLQUERY=~/cclaro.classified.short
PFAM=/mnt/d/pfam32/
RENAME=cclaro
PFAMTE=Pfam.TE.accessions_Release32.amended
#========================START OF PROTEIN PREDITION============================
blastx -version
echo "running blastx to generate Amino acid predictions from plus strand of consensus sequences"
blastx -max_target_seqs 1 -num_threads 8 -db $PROREF -query $RRLQUERY -outfmt 14 -out $RRLQUERY.format14.plus -evalue 1e-10 -strand plus
echo "running blastx to generate Amino acid predictions from minus strand of consensus sequences"
blastx -max_target_seqs 1 -num_threads 8 -db $PROREF -query $RRLQUERY -outfmt 14 -out $RRLQUERY.format14.minus -evalue 1e-10 -strand minus
#plus strand cleanup
echo "perform cleanup on plus strand predictions"
for xml in $RRLQUERY.format14.plus*
do hash=$(grep qseq $xml|wc -l)
bname=$(grep query-title $xml|sed "s/.*$RENAME/$RENAME/"|sed 's/<.*//g')
if [ $hash != 0 ]; then
echo $bname $hash
for ((i=1;i<=$hash;i++));do echo '>'${bname}_plus_qseq_$i'' >> ${bname}_${hash}.plus.pre
awk '/qseq/{i++}i=='$i'' $xml|sed -n '/hseq/q;p'|sed s'/qseq>//g'|sed s'/<//g'|sed s'/\///g'|sed s'/-//'g|sed s'/\*//'g|sed s'/X//'g|tr -d " \t" >> ${bname}_${hash}.plus.pre
done
fi
done
cat *.plus.pre > ~/$RENAME.plus.predicted
#minus strand cleanup
echo "perform cleanup on minus strand predictions"
for xml in $RRLQUERY.format14.minus*
do hash=$(grep qseq $xml|wc -l)
bname=$(grep query-title $xml|sed "s/.*$RENAME/$RENAME/"|sed 's/<.*//g')
if [ $hash != 0 ]; then
echo $bname $hash
for ((i=1;i<=$hash;i++));do echo '>'${bname}_minus_qseq_$i'' >> ${bname}_${hash}.minus.pre
awk '/qseq/{i++}i=='$i'' $xml|sed -n '/hseq/q;p'|sed s'/qseq>//g'|sed s'/<//g'|sed s'/\///g'|sed s'/-//'g|sed s'/\*//'g|sed s'/X//'g|tr -d " \t" >> ${bname}_${hash}.minus.pre
done
fi
done
cat *.minus.pre > ~/$RENAME.minus.predicted
echo "cleaning up"
rm *.pre
rm *.xml
rm $RRLQUERY.format14.plus $RRLQUERY.format14.minus
#========================END OF PROTEIN PREDICTION============================

#========================START OF PFAM_SCAN===================================
pfam_scan.pl -fasta $RENAME.plus.predicted -dir $PFAM -outfile $RENAME.plus.predicted.pfamtbl -e_seq 0.001
pfam_scan.pl -fasta $RENAME.minus.predicted -dir $PFAM -outfile $RENAME.minus.predicted.pfamtbl -e_seq 0.001
#========================END OF PFAM_SCAN=====================================

#========================START OF ANNOTATION==================================
#========FILTERING==================================
printf "subsetting the pfamtbl output to unclassified consensus only. \n"
grep -i \#unknown $RENAME.plus.predicted.pfamtbl > $RENAME.plus.predicted.pfamtbl.unclassified
grep -i \#unknown $RENAME.minus.predicted.pfamtbl > $RENAME.minus.predicted.pfamtbl.unclassified
printf "Detect unclassified consensi that have TE domains. Plus strand:\n`grep -i -w -f $PFAMTE $RENAME.plus.predicted.pfamtbl.unclassified` \n"
printf "Detect unclassified consensi that have TE domains. Minus strand:\n`grep -i -w -f $PFAMTE $RENAME.minus.predicted.pfamtbl.unclassified` \n"
grep -i -w -f $PFAMTE $RENAME.plus.predicted.pfamtbl.unclassified |awk '{print $1}'|sed 's/\#.*//g'|uniq > 001.tem
grep -i -w -f $PFAMTE $RENAME.minus.predicted.pfamtbl.unclassified|awk '{print $1}'|sed 's/\#.*//g'|uniq > 002.tem
cat 001.tem 002.tem|uniq > $RENAME.Unclassified_consensus_TEs.ids #can be used to annotate unknown consensus to TEs
grep -i -w -f $PFAMTE $RENAME.plus.predicted.pfamtbl.unclassified > $RENAME.Unclassified_consensus_TEs
grep -i -w -f $PFAMTE $RENAME.minus.predicted.pfamtbl.unclassified >> $RENAME.Unclassified_consensus_TEs
rm *.tem
grep -vf  $RENAME.Unclassified_consensus_TEs.ids $RENAME.plus.predicted.pfamtbl.unclassified > $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged
grep -vf  $RENAME.Unclassified_consensus_TEs.ids $RENAME.minus.predicted.pfamtbl.unclassified > $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged
printf "Gathering unclassified consensus sequences that have domains detected in both strands.\n"
printf "These consensus sequences are tricky to annotate.\nPut aside for anyone that are interested. Check file consensus.proteins.both.strand"
grep -i \#unknown $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged|awk '{print $1}'|sed 's/\#.*//g'|uniq > 003.tem
grep -i \#unknown $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged|awk '{print $1}'|sed 's/\#.*//g'|uniq > 004.tem
printf "These consensus sequences are:\n`grep -f 003.tem 004.tem` \n"
grep -f 003.tem 004.tem|uniq > $RENAME.consensus.proteins.both.strand
rm  *.tem
printf "Now short matches that are shorter than 100AA are removed"
grep -vf $RENAME.consensus.proteins.both.strand $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged > $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged.mono
awk '{ if ($11 >= 100) { print } }' $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged.mono > c$RENAME.plus.predicted.pfamtbl.unclassified.TEpurged.mono.final
grep -vf $RENAME.consensus.proteins.both.strand $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged > $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged.mono
awk '{ if ($11 >= 100) { print } }' $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged.mono > $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged.mono.final
#========FILTERING_End==================================

#========ANNOTATION==================================
printf "Annotation initiated.\nBacking up the original library file.\n`cp $RRLQUERY $RRLQUERY.bk`"
cat $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged.mono.final|awk '{print $1}'|sed 's/\#.*//g'|uniq > 005.tem
cat $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged.mono.final|awk '{print $1}'|sed 's/\#.*//g'|uniq > 006.tem
printf "Following repeat models were renamed based on plus strand evidence:\n"
for i in `cat 005.tem`;do echo $i >> $RENAME.renamed; grep $i $RENAME.plus.predicted.pfamtbl.unclassified.TEpurged.mono.final|awk '{print $7}' >> $i.plus.annotate.tem; done
for i in *.plus.annotate.tem; do namehash=$(cat $i|sort|uniq|paste -s -d '_'); oldname=$(grep ${i%.plus.annotate.tem}\# $RRLQUERY |sed 's/>//g'); sed -i "s|$oldname|${oldname%Unknown}$namehash|g" $RRLQUERY; done
printf "Following repeat models were renamed based on minus strand evidence:\n"
for i in `cat 006.tem`;do echo $i >> $RENAME.renamed; grep $i $RENAME.minus.predicted.pfamtbl.unclassified.TEpurged.mono.final|awk '{print $7}' >> $i.minus.annotate.tem; done
for i in *.minus.annotate.tem; do namehash=$(cat $i|sort|uniq|paste -s -d '_'); oldname=$(grep ${i%.minus.annotate.tem}\# $RRLQUERY |sed 's/>//g'); sed -i "s|$oldname|${oldname%Unknown}$namehash|g" $RRLQUERY; done
#========ANNOTATION_END==================================
rm *.tem
#=====================================END====================================
