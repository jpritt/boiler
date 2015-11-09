================
 For Submission
================
* Compare query times to samtools (?)
    * Probably need to improve query times -- compress in smaller chunks
* Regenerate results
    * Run Boiler on all samples
    * Run all accuracy measures
* Revamp introduction

Questions
* Stringtie or Cufflinks primarily? (I'm leaning Stringtie)
* What about paired reads with different XS values? Are they discordant?
* What to do with discordant reads? (Currently we are unpairing them)
* Current datasets: 5 simulated Drosophila, 1 simulated human, 5 Geuvadis
    * How to present results across these datasets?

Most important figures:
* Compression timing comparision
* Memory Usage comparison
* Fragment length distribution
* Reference-based accuracy
* Non-reference-based accuracy
* Query timing comparision and memory footprint

Other relevant figures:
* SAM alignment precision/recall
* WKR across samples
* WKR for varying k
* Tripartite Accuracy

============
 Other TODO (not top priority)
============
* Die gracefully when input file lacks SAM header
* Handle HISAT output, where unaligned mates can appear in accepted_hits.sam
* Some scripts not compatible with Python 2 (e.g. use of multi-param `bytes`
  function), others are not compatible with Python 3 (e.g. use of `xrange`)
    * struct.pack_into(), unpack_from()
    * http://stackoverflow.com/questions/30402743/python-2-7-equivalent-of-built-in-method-int-from-bytes

* Check whether HISAT output contains discordant reads
