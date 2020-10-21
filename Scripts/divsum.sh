#!/bin/bash
#title           :trea
#description     :filter out host proteins from unclassified repeat consensus
#author		 :wenbo.xu@botinst.uzh.ch
#date            :20201011
#version         :0.3
#usage		 :./hostfinder.sh
 #==============================================================================
DIV=cclaro_1_nanopore_idrenamed.align.devil.divsum
cat $DIV|sed '1,/Coverage\ for/d' > $DIV.simple
sed '/Coverage\ for/q' $DIV|sed '/^$/d'|sed '$d'|sed '1,4d'|sed '2d' > $DIV.clean
head -1 $DIV.clean >> $DIV.clean.nonTE
grep -wf RhiirA1_1 RhiirA1_1_AssemblyScaffolds.fasta.align.divsum.clean >> RhiirA1_1_AssemblyScaffolds.fasta.align.divsum.clean.nonTE
grep -vwEf RhiirA1_1 RhiirA1_1_AssemblyScaffolds.fasta.align.divsum.clean > destinationfile
sed 's/Unknown/Non-TE-Repeat/g' RhiirA1_1_AssemblyScaffolds.fasta.align.divsum.clean.nonTE |sed '1d' >> destinationfile
