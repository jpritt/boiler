#!/bin/bash

# $1: input SAM file
# $2: output SAM prefix
$BOILER_HOME/inferXStags.py $1 $2.sam
samtools view -bS $2.sam | samtools sort - $2
samtools view -h -o $2.sam $2.bam

