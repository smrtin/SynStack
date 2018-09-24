#!/bin/bash

#part1: flip 
#our gene of interest might be located on the + or the - strand. 
#in order to produce a well-arranged graphical output all genes should be oriented in the same direction 

#part2: unify
#to simplify the graphical output all genes become the same size


#usage: bash 03_Flip.bash proteinfile.fasta

SEQFILE=$1

#files from previous steps:
GFF_output="$SEQFILE.synt.gff"
PAIRS_output="$SEQFILE.pairs"

#part1:
cat $SEQFILE.synt.gff | sort | uniq > $SEQFILE.synt2.GFF

OUTPUT="$SEQFILE.synt3.GFF"

printf '' >$OUTPUT

while read -r line ; 
    do
        echo $line
        sleep 1
    Scaffold1=$(echo $line | cut -f3 -d' ' )
    orient=$(echo $line | cut -f4 -d' ')
    
    if [[ "$orient" == '+' ]]; then
    grep "$Scaffold1" $SEQFILE.synt2.GFF  >>$OUTPUT
    
    else
    echo $Scaffold1 is in wrong direction... 
    HIGHEST=$(grep "$Scaffold1" $SEQFILE.synt2.GFF | sort -nr -k5 | head -1 | cut -f5 -d' ')
    
    while read -r secline ;
        do
        array=(${secline// / })	#split line into array  
        echo
        echo $secline
        echo -e "${array[0]} ${array[1]} ${array[2]} ${array[3]} ${array[4]} ${array[5]} ${array[6]} ${array[7]} ${array[8]}"
        newStart=$(expr $HIGHEST - ${array[4]} + 1 )
        newSTOP=$(expr $HIGHEST - ${array[3]} + 1 )
    
            if [ ${array[6]} == '+' ] ; then
            newStrand='-'
            else
            newStrand='+'
            fi

        echo -e "${array[0]} ${array[1]} ${array[2]} $newStart $newSTOP ${array[5]} $newStrand ${array[7]} ${array[8]}" >> $OUTPUT
    

        done < <(grep "$Scaffold1" $SEQFILE.synt2.GFF)
    
    fi
    
    
    done < <(cat $PAIRS_output)

    
#part2: to simplify the graphical output and give all genes a unified size 

OUT="$SEQ.synt3.GFF.unify"

printf '' > $OUT
#for each scaffold, start with the lowest position... adjust in whileloop 
for SCAFFOLD in $(cut -f1 -d' ' $SEQ.synt3.GFF | sort | uniq ); do
echo this is scaffold $SCAFFOLD
numGenes=$(grep "$SCAFFOLD" $SEQ.synt3.GFF | wc -l )

 while read -r line ; 
    do
    array=(${line// / })
    
    echo -e "${array[0]} ${array[1]} ${array[2]} ${array[3]} ${array[4]} ${array[5]} ${array[6]} ${array[7]} ${array[8]}"
    echo first
    newStart=$(($numGenes * 500)) ####100
    echo second
    newStop=$(expr $newStart + 400) #70
    
    echo -e "${array[0]} ${array[1]} ${array[2]} $newStart $newStop ${array[5]} ${array[6]} ${array[7]} ${array[8]}" >> $OUT
    
    numGenes=$(expr $numGenes - 1 )
    echo number of genes.... $numGenes
    
    done< <(grep "$SCAFFOLD" $SEQ.synt3.GFF  | sort -nr -k5 ) #wir fangen von hinten an....
done