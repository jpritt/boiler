================
 For Submission
================
* Revamp introduction
* Compare query times to samtools (?)
    * Probably need to improve query times -- compress in smaller chunks
* Regenerate results
    * Run Boiler on all samples
    * Compute all accuracy measures
* Check whether HISAT output contains discordant reads

Questions
* Stringtie or Cufflinks primarily? (I'm leaning Stringtie)
* What about paired reads with different XS values? Are they discordant?
* What to do with discordant reads? (Currently we are unpairing them)


Most important figures:
* Compression timing comparision
* Fragment length distribution
* Reference-based accuracy
* Non-reference-based accuracy
* Query timing comparision

Other relevant figures:
* SAM alignment precision/recall
* WKR across samples
* WKR for varying k

============
 Other TODO
============
* Die gracefully when input file lacks SAM header
* Handle HISAT output, where unaligned mates can appear in accepted_hits.sam
* Some scripts not compatible with Python 2 (e.g. use of multi-param `bytes`
  function), others are not compatible with Python 3 (e.g. use of `xrange`)
