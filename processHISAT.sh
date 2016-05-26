#!/bin/bash

# $1: input SAM file
# $2: output SAM prefix
$BOILER_HOME/removeUp.py $1 temp1.sam
$BOILER_HOME/enumeratePairs.py --input temp1.sam --output temp2.sam
$BOILER_HOME/inferXStags.py temp2.sam $2.sam
samtools view -bS $2.sam | samtools sort - $2
samtools view -h -o $2.sam $2.bam
rm temp1.sam
rm temp2.sam

