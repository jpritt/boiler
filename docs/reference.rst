Reference
=========

Boiler has three modes, each described in more detail below:

#. :ref:`compress` -- compress a SAM file. 
#. :ref:`query` -- query a compressed file.
#. :ref:`decompress` -- expand a compressed file, outputting a SAM file.

Boiler is invoked by entering ::

    python3 boiler.py <mode> <[args]>

To run with PyPy instead of Python, run ::

    pypy boiler.py <mode> <[args]>

To get help for a given mode, enter ::

    python3 boiler.py <mode> -h

.. _compress:

========
compress
========

Boiler requires a SAM file to compress. To convert a BAM file to SAM with SAMtools, run::

    samtools view -h -o path/to/alignments.bam path/to/alignments.sam


To compress a SAM file with Boiler, run the following command::

    python3 boiler.py compress <[args]> path/to/alignments.sam path/to/compressed.bl 

The following optional arguments are available:

``-c/--frag-len-cutoff <threshold>``

    As a first step in compressing, Boiler groups overlapping reads into 'bundles' using a similar method to Cufflinks (see the :ref:`bundles` query below for more details). In paired-end datasets, some mates are mapped millions of bases apart or even on different chromosomes. Boiler stores such pairs as bundle-spanning reads, rather than creating massively long bundles to suit them. If ``frag-len-cutoff`` is set, pairs longer than this cutoff will not contribute to determining bundles, so they will often be stored as bundle-spanning pairs.

    Changing the threshold for ``--frag-len-cutoff`` will not affect accuracy. Generally, decreasing the threshold will lead to faster compression, decompression, and query times, but also to larger file size, and vice versa. See ``--frag-len-z-cutoff`` below for an alternative to setting the threshold directly.

``-z/--frag-len-z-cutoff <z-score>``

    As an alternative to setting ``--frag-len-cutoff`` directly, if ``--frag-len-z-cutoff`` is set Boiler will perform a first pass over the reads to establish the average and standard distribution of all paired read lengths. Any reads with a `z-score <https://en.wikipedia.org/wiki/Standard_score>`_ greater than ``--frag-len-z-cutoff`` will not contribute to determining bundles, so they will often be stored as bundle-spanning pairs. If neither ``--frag-len-cutoff`` nor ``--frag-len-z-cutoff`` is set, Boiler will set ``--frag-len-z-cutoff`` to ``0.125``, which we have found to work well in practice.

``-s/--split-diff-strands``

    Sometimes a SAM file contains paired mates that lie on different chromosomes. Boiler preserves these pairs by default; use ``--split-diff-strands`` to convert them to unpaired reads.

``-d/--split-discordant``

    SAM files often contain discordant pairs, i.e. paires where one mate intersects an intron of the other mate. Boiler preserves these pairs by default; use ``--split-discordant`` to convert them to unpaired reads. 

``-p/--preprocess``

    This argument should be added for alignments produced by HISAT, which require an additional processing step before compression. For multimapped reads, HISAT outputs a list of ``n`` left mates and ``m`` right mates, where any left mate may be be paired with any right mate. In contrast, Cufflinks outputs a pair for each unique combination of possible left and right mates. Boiler requires pairs to be enumerated, as in Cufflinks output.

    Alternatively, you can preprocess HISAT alignments yourself by running::

        python enumeratePairs.py --input alignments.sam --output alignments.processed.sam
        samtools sort -bS alignments.processed.sam | samtools sort - alignments.processed
        samtools view -h -o alignments.processed.sam alignments.processed.bam

    Following which you can run Boiler as normal::

        python3 boiler.py compress <[args]> alignments.processed.sam path/to/compressed.bl

``-g/--gtf <path/to/transcripts.gtf>``

    Boiler offers the option of using a reference gtf file to guide compression. Boiler adds additional splice sites at every transcript splice site and endpoint in the gtf. This improves the accuracy of read recovery at the cost of a significant size increase.

``-v/--verbose``

    Print additional debug information.

.. _query:

=====
query
=====

Boiler currently supports the following queries:

#. :ref:`bundles`
#. :ref:`coverage`
#. :ref:`reads`
#. :ref:`counts`

The first 3 queries require a chromsome and optional start and end position. Boiler will return the results over the query over the given interval. If ``--start`` or ``--end`` are not specified, the endpoints of the given chromsome will be used. Results will be written to the given ``out_file``.  To run on of these queries, enter::

    python3 boiler.py [--bundles | --coverage | --reads] --chrom <c> --start <s> --end <e> path/to/compressed.bl out_file

The ``counts`` query takes a gtf file as input, but no chromosome or range. Results will be returned for the entire genome. To run this query, enter::

    python3 boiler.py --counts --gtf <path/to/transcripts.gtf> path/to/compressed.bl out_file

.. _bundles:

bundles
^^^^^^^

'Bundles' divide the genome into manageable chunks, roughly corresponding to potential gene boundaries. Boiler calculates bundles in a similar way to Cufflinks; as read are processed in order of position, the end of the current bundle is extended to the end of the current paired-end read. If the next read begins more than 50 bases after the end of the current bundle, a new bundle is creating beginning at the current read. 

.. _coverage:

coverage
^^^^^^^^

This query prints a vector containing the total coverage at each base in the given range.

.. _reads:

reads
^^^^^

This query outputs a SAM file containing all reads that overlap the given range.

.. _counts:

counts
^^^^^^

Boiler parses the gtf file and extracts a list of exons and junctions

.. _decompress:

==========
decompress
==========

To decompress a file with Boiler, run ::

    python3 boiler.py decompress <[args]> path/to/compressed.bl expanded.sam

The output SAM file is not sorted; to convert to a sorted BAM file, enter ::

    samtools view -bS expanded.sam | samtools sort - expanded

Then you can then sort the expanded SAM file by running ::

    samtools view -h -o expanded.sam expanded.bam

The following arguments are available for decompression:

``--force-xs``

    Boiler will assign a random XS value to any spliced reads that do not have an XS tag. This is meant for some tools such as Cufflinks which require that all spliced reads have an XS tag.

``-v/--verbose``

    Print additional debug information.