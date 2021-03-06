Tutorial
========

To begin, download the latest version of Boiler from https://github.com/jpritt/boiler. Add the main directory to your path and make sure you have Python version 3 or higher. You will also need SAMtools, which you can download from samtools.sourceforge.net.

Download the SAM dataset `here <https://drive.google.com/open?id=0B-3BXIpgKqSXNkdnVjl4TnhzdkE>`_ and move it to your working directory.

Run ::

    mkdir compressed
    python3 boiler.py compress --frag-len-z-cutoff 0.125 accepted_hits.sam compressed/compressed.bl

If all goes well, you should see something like this (exact output may change with future versions)::

    Set fragment length cutoff to z=0.125000 (33165) based on length distribution
    0.84 % of pairs are longer than the cutoff
    Using fragment length cutoff of 33165
    Not splitting mates on different strands
    Not splitting discordant
    0 cross-bundle reads unmatched
    Minimum bundle length: 12
    Maximum bundle length: 206957
    Average bundle length: 2514
    1097 cross-bundle buckets
    Compressed size: 29682
    Approximately 3979761 / 6972093 = 57.081295% of compressed file is coverage
    Finished compressing

You should now have a file ``compressed/compressed.bl`` roughly 4.3 MB in size.

Now let's query all of the bundles that Boiler found in chromosome 2L::

    python3 boiler.py query --bundles --chrom 2L compressed/compressed.bl bundles.txt

You should now have a file ``bundles.txt`` containing all of the bundles used by Boiler. Type ::

    head bundles.txt

to see the first few lines of this file::

    7478	9485
    9841	21430
    21825	23108
    23180	24034
    24856	25219
    25404	26251
    26333	33987
    34045	35094
    36182	37317
    37538	37931

To query the coverage in the first bundle, run ::

    python3 boiler.py query --coverage --chrom 2L --start 7478 --end 9485 compressed/compressed.bl coverage.txt

``coverage.txt`` should now contain a comma-separated vector containing the coverage at every base in the interval [7478, 9485). Finally, to query the reads in the first bundle, run ::

    python3 boiler.py query --reads --chrom 2L --start 7478 --end 9485 compressed/compressed.bl reads.sam

``reads.sam`` is a SAM file with no header, containing all the aligned reads in the interval [7478, 9485). Type ::

    head reads.sam

to see the first few reads in this bundle, which should look like this::

    2L:0	0	2L	7772	50	76M	*	0	0	*	*	NH:i:1
    2L:1	0	2L	7795	50	76M	*	0	0	*	*	NH:i:1
    2L:2	0	2L	7808	50	76M	*	0	0	*	*	NH:i:1
    2L:3	0	2L	7863	50	76M	*	0	0	*	*	NH:i:1
    2L:4	0	2L	8073	50	44M112N32M	*	0	0	*	*	XS:A:+	NH:i:1
    2L:5	0	2L	8595	50	76M	*	0	0	*	*	NH:i:1
    2L:6	0	2L	8781	50	76M	*	0	0	*	*	NH:i:1
    2L:7	0	2L	8852	50	76M	*	0	0	*	*	NH:i:1
    2L:8	0	2L	8963	50	76M	*	0	0	*	*	NH:i:1
    2L:9	0	2L	8969	50	76M	*	0	0	*	*	NH:i:1

Finally, let's decompress the compressed file by running ::

    python3 boiler.py decompress compressed/compressed.bl expanded.sam

The resulting SAM file is unsorted -- to sort and convert it to BAM, run ::

    samtools view -bS expanded.sam | samtools sort - expanded
