class Read:
    def __init__(self, chrom, exons, strand=None, NH=1):
        ''' exon is a list of (start,end) tuples marking each exonic region of this read.
            xs indicates the strand on which this gene lies, as described in the SAM file format.
        '''
        self.chrom = chrom
        self.exons = exons
        self.strand = strand

        self.lenLeft = 0
        self.lenRight = 0

        self.NH = NH