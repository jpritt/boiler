compress-alignments
===================

'compress.py' is the main script that runs compression and decompression.
To run, first make a directory 'compressed' for compressed files.

> mkdir compressed
> ./compress --alignments path/to/alignments.sam

Compressed files will be generated in compressed/
The decompressed SAM file will appear in the current directory, named expanded.sam.

To sort and convert to BAM, run:
> samtools view -bS expanded.sam | samtools sort - expanded


To compare 2 cufflinks files, run:
> ./compareGTFs transcripts1.gtf transcripts2.gtf
 
