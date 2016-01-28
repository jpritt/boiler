compress-alignments
===================

'boiler.py' is the main script that runs compression and decompression.
To compress, run the following:.

> ./boiler.py --compress --alignments path/to/alignments.sam --output path/to/compressed.bin [--frag-len-z-cutoff 0.125] [--split-discordant] [--split-diff-strands]

--frag-len-z-cutoff sets the z-score for paired-end read lengths at which to set the cutoff for placing mates in different bundles. 0.125 seems to be a good z-score. Alternatively, you can use --frag-len-cutoff to set the cutoff directly.
If --split-discordant is present, discordant reads will be treated as unpaired reads.
If --split-diff-strands is present, reads with contradicting XS values will be treated as unpaired reads.

To decompress, run the following:

> ./boiler.py --decompress --output path/to/compressed.bin --expand-to path/to/expanded.sam

The decompressed SAM file will appear in the given directory, named expanded.sam.
Compression and decompression can also be combined into a single command:
> ./boiler --compress --alignments path/to/alignments.sam --output path/to/compressed.bin --expand-to path/to/expanded.sam [--frag-len-z-cutoff 0.125] [--split-discordant] [--split-diff-strands]


To sort and convert to BAM, run:
> samtools view -bS expanded.sam | samtools sort - expanded

To compare 2 cufflinks files, run:
> ./compareGTFs.py transcripts1.gtf transcripts2.gtf

To query a compressed file for bundles, coverage, or reads:

> ./query.py --compressed path/to/compressed.bin [--bundles | --coverage | --reads] --chrom c [--start s] [--end e]

If --start or --end is absent, Boiler will use the beginning or end of the chromsome as bounds on the query.
