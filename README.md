compress-alignments
===================

'boiler.py' is the main script that runs compression and decompression. Python 3 is required to run Boiler. The input SAM file must be sorted by read start position.
To compress, run the following:.

> ./boiler.py compress [--frag-len-z-cutoff 0.125] [--split-discordant] [--split-diff-strands] [--preprocess tophat | stringtie] path/to/alignments.sam path/to/compressed.bl

--frag-len-z-cutoff sets the z-score for paired-end read lengths at which to set the cutoff for placing mates in different bundles. 0.125 seems to be a good z-score. Alternatively, you can use --frag-len-cutoff to set the cutoff directly.
If --split-discordant is present, discordant reads will be treated as unpaired reads.
If --split-diff-strands is present, reads with contradicting XS values will be treated as unpaired reads.


To decompress, run the following:

> ./boiler.py decompress [--force-xs] path/to/compressed.bl path/to/expanded.sam

--force-xs will assign XS tags to all spliced reads, as required by Cufflinks. If spliced reads are found with XS tags, they will be assigned at random.
The decompressed SAM file will appear in the given directory, named expanded.sam.


To sort and convert to BAM, run:
> samtools view -bS expanded.sam | samtools sort - expanded

To compare 2 cufflinks files, run:
> ./compareGTFs.py transcripts1.gtf transcripts2.gtf


To query a compressed file for bundles, coverage, or reads:

> ./boiler.py query [--bundles | --coverage | --reads] --chrom c [--start s] [--end e] path/to/compressed.bl path/to/output

If no output argument is provided, standard output will be used
If --start or --end is absent, Boiler will use the beginning or end of the chromsome as bounds on the query.
