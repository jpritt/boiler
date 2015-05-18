* Die gracefully when input file lacks SAM header
* Handle HISAT output, where unaligned mates can appear in accepted_hits.sam
* Some scripts not compatible with Python 2 (e.g. use of multi-param `bytes`
  function), others are not compatible with Python 3 (e.g. use of `xrange`)
