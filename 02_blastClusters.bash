#!/bin/bash 


#build an all-vs-all output with blastp and create an abc-file for mcl
#identify clusters with mcl

#usage: bash 02_blastClusters.bash proteinSequences.fasta


SeqFile=$1

#start BLAST
makeblastdb -in $SeqFile -dbtype prot -out $SeqFile.blastDB

blastp -db $SeqFile.blastDB -query $SeqFile -out $SeqFile.blastp.out -word_size 5 -outfmt '6 qseqid sseqid evalue bitscore' -num_threads 6 -threshold 50 #is threshold better than evalue? does it show bitscore?

cut -f1,2,3 $SeqFile.blastp.out > $SeqFile.blastp.abc


#start mcl analysis
mcxload -abc $SeqFile.blastp.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o $SeqFile.bp-protein_hits.mci -write-tab $SeqFile.bp-protein_hits.tab

#do the clustering with different I-values, which affect granularity of the clustering

mcl $SeqFile.bp-protein_hits.mci -I 1.4 -use-tab $SeqFile.bp-protein_hits.tab -o $SeqFile.bp-clustI14.out 
