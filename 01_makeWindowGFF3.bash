#!/bin/bash

#this script will extract all information from the gff file, which is annotated in the window around a specific gene of interest

#usage: bash 01_makeWindowGFF3.bash Sequencefile.fasta 50000

SEQFILE=$1 #sequence file or ordered list of titles...sequence file must be fasta format and the sequence header must start with four letter species identifier seperated with an underscore
WINDOW=$2 #size of genomic region to be included in the analysis

GFF_FOLDER='../../GFFs/keep' #location where all the gff-files are located. file name must begin with four letter species identifier and files are compressed SPEC_filename.gff3.gz
PROT_FOLDER='../../Proteins/modifiedTitles' #location where all protein files are located. file name must begin with four letter species identifier.  sequence header must also start with species id and the proteinsequence must be concatenated in one line after the header.

#output
PAIRS_output="$SEQFILE.pairs"
printf '' > $PAIRS_output
GFF_output="$SEQFILE.synt.gff"
printf '' >$GFF_output
PROT_output="$SEQFILE.synt.longest_prot.fa"
printf '' > $PROT_output

############
headerPresent=$(grep '>' $SEQFILE)

if [[ $headerPresent ]]; then
    TITLES=$( grep '>' $SEQFILE | sed 's/>//gi'  )  
    NUMBEROFSEQUENCES=$(grep -c '>' $SEQFILE)
    echo File contains $NUMBEROFSEQUENCES sequences...
    sleep 2
else
    TITLES=$(cat $SEQFILE)
    numberOfTITLE=$(cat $SEQFILE | wc -l)
    echo file contains $numberOfTITLE titles ...
fi

for SEQUENCE in $TITLES #get all the sequence titles
do 
echo $SEQUENCE
SPECIES=$( echo $SEQUENCE | sed 's/\(\w\w\w\w\)_\(.*\)/\1 \2/' | cut -d' ' -f1 ) #get the species identifier
SEQID=$( echo $SEQUENCE | sed 's/\(\w\w\w\w\)_\(.*\)/\1 \2/' | cut -d' ' -f2 ) #get the sequence identifier from the original dataset


THEGREP=$(zgrep "protein_id=$SEQID" $GFF_FOLDER/$SPECIES*) #go for protein_id= this should be present in all gff files.... but not in GBR....

if [[ $THEGREP ]] ; then #if the grep was successfull...

#echo we have a hit
mRNA_ID=$(echo $THEGREP | tr ";" "\n" | grep 'Parent=' | sort | uniq | sed 's/Parent=//' ) #identify the parent ID (the RNA-ID) of the CDS that the protein_id
numParents=$( echo $mRNA_ID | wc -l )
    
    if [[ "$numParents" > 1 ]]; then
    echo \#############################
    echo $SPECIES has in sequence $SEQID more than one parent this might lead to problems ...
    echo parents $mRNA_ID
    echo have a closer look at the manual annotation sometimes CDSs from different scaffolds are joined... #this happend in nasonia... partial genes....
    echo \############################
    
    else
        
    GENE_ID=$(zgrep "ID=$mRNA_ID;" $GFF_FOLDER/$SPECIES* | tr ";" "\n" | grep 'Parent=' | sort | uniq | sed 's/Parent=//') #do another grep and identify the Gene ID
    GENENAME=$(zgrep "ID=$GENE_ID;" $GFF_FOLDER/$SPECIES* | grep -P '\tgene\t' | sed -r 's/;/\n/gi' | grep 'Name=' | sed 's/Name=//gi' | sed 's/\\/\\\\/gi')
    SCAFFOLD=$(echo $THEGREP | cut -f1 -d' ' | head -1 )
    orientation=$(zgrep "ID=$GENE_ID;" $GFF_FOLDER/$SPECIES* | grep -P '\tgene\t' | cut -f7) 
    
    printf "$SEQUENCE\t${SPECIES}_${GENE_ID}\t${SPECIES}_${GENE_ID}__${SCAFFOLD}\t$orientation\n" >> $PAIRS_output
    
    Window_up=$(zgrep "ID=$GENE_ID;" $GFF_FOLDER/$SPECIES* | grep -P '\tgene\t' | cut  -f5)
    Window_down=$(zgrep "ID=$GENE_ID;" $GFF_FOLDER/$SPECIES* | grep -P '\tgene\t' | cut  -f4)
    Window_up_limit=$(($Window_up+$WINDOW))
    Window_down_limit=$(($Window_down-$WINDOW))
    if (( "$Window_down_limit" < 1 )); then
    Window_down_limit=0
    fi
    
    #go through all other genes on the scaffold....
    while read -r line ; 
    do
    up=$(echo $line | cut -f5 -d' ' )
    down=$(echo $line | cut -f4 -d' ')
    
        if (( $up < $Window_up_limit && $down > $Window_down_limit )); then 
        
            myGENE_ID=$(echo $line | cut -f9 -d' ' |  tr ";" "\n" | grep 'ID=' | sed 's/ID=//')
            mymRNA_ID=$(zgrep -P "^$SCAFFOLD\t" $GFF_FOLDER/$SPECIES* | grep -P '\tmRNA\t' | grep "Parent=$myGENE_ID" |  tr ";" "\n" | grep 'ID=' | sed 's/.*ID=//' ) #these can be multiple....
            
            if [[ -z $mymRNA_ID ]]; then echo there is no mRNA for gene $myGENE_ID need to continue... ;continue; fi #check if mRNA is present... if not continue in while-loop....
            
            numTranscripts=$( echo $mymRNA_ID | wc -w )
            
            echo $up $down gene $myGENE_ID mrna $numTranscripts $mymRNA_ID #
            
            if [[ "$numTranscripts" > 1 ]]; then #if we have multiple transcripts we have to find the one with the longest AA-sequence
                
                PROTEINS=''
                
               for Transcript in $( echo $mymRNA_ID); do #go through list of Transcripts and get the protein length
                
                myProt_ID=$(zgrep -P "^$SCAFFOLD\t" $GFF_FOLDER/$SPECIES* | grep -P '\tCDS\t' | grep "Parent=$Transcript" |  tr ";" "\n" | grep 'protein_id=' | sed 's/.*protein_id=//' | sort | uniq )
                PROTEINLENGTH=$( grep  -A 1 "$myProt_ID" $PROT_FOLDER/$SPECIES* | grep -v '>' | wc -m )
                
                PROTEINS=$(printf "$PROTEINS\n$myGENE_ID\t$Transcript\t$myProt_ID\t$PROTEINLENGTH\n")
                #echo gene $myGENE_ID mrna multi $mymRNA_ID protein $myProt_ID proteinlength $PROTEINLENGTH  
                
               done
                printf "$PROTEINS" | sort -nr -k4   #sort all proteinIDs after sequnces length
                echo longest:
                LONGEST=$(printf "$PROTEINS" | sort -nr -k4 | head -n1 | cut -f3)  #sort all proteinIDs after sequnces length
                grep  -A 1 "$LONGEST" $PROT_FOLDER/$SPECIES* | sed -r "s/(>.*)/\1_geneID_$SPECIES\_${myGENE_ID}_GOI_$SPECIES\_$SEQID/" >> $PROT_output
                
            else
            
                myProt_ID=$(zgrep -P "^$SCAFFOLD\t" $GFF_FOLDER/$SPECIES* | grep -P '\tCDS\t' | grep "Parent=$mymRNA_ID" |  tr ";" "\n" | grep 'protein_id=' | sed 's/.*protein_id=//' | sort | uniq )
                PROTEINLENGTH=$( grep  -A 1 "$myProt_ID" $PROT_FOLDER/$SPECIES* | grep -v '>' | wc -m )
                grep  -A 1 "$myProt_ID" $PROT_FOLDER/$SPECIES* | sed -r "s/(>.*)/\1_geneID_$SPECIES\_${myGENE_ID}_GOI_$SPECIES\_$SEQID/" >> $PROT_output
                
                echo gene $myGENE_ID mrna $mymRNA_ID protein $myProt_ID proteinlength $PROTEINLENGTH        
            fi
    
        #GFF OUTPUT of gene...
            echo $line | sed  "s/^$SCAFFOLD/$SPECIES\_$GENE_ID\__$SCAFFOLD/" - | sed "s/ID=$myGENE_ID/ID=$SPECIES\_$myGENE_ID/" - >> $GFF_output    #seperater species gene id and scaffoldID are two _ 
        fi
    
    done < <(zgrep -P "^$SCAFFOLD\t" $GFF_FOLDER/$SPECIES* | grep -P '\tgene\t' ) #all genes on the same scaffold...
    
    echo $GENE_ID $mRNA_ID Window_down: $Window_down limit $Window_down_limit Window_up: $Window_up limit $Window_up_limit
    fi

fi
done

###some times if two genes of same family appear on same scaffold in neighborhood this region will appear double this must be removed...

#also eliminate doubles from the fasta