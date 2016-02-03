#! /usr/bin/env python3

import resource
import subprocess
import threading

max_mem = 0
def finish():
    print('Maximum memory: %0.1f MB' % (max_mem / 1000.0))
    exit()

def run():
    proc = subprocess.Popen(['java', '-Xmx16g', '-jar', '/scratch0/langmead-fs1/shared/cramtools/cramtools-3.0.jar', 'cram', '-I', 'tophat_out/accepted_hits_fixed.bam', '-R', '/scratch0/langmead-fs1/user/jacob/genomes/Drosophila_melanogaster/Ensembl/BDGP5/Sequence/WholeGenomeFasta/genome.fa', '-O', 'compressed/compressed.cram'])
    proc.wait()
    finish()

    '''
    max_mem = 0
    while proc.poll() is None:
        mem = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        print(mem)
        if mem > max_mem:
            max_mem = mem

    print(max_mem)
    '''

thread = threading.Thread(target=run)
thread.start()

while True:
    mem = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if mem > max_mem:
        max_mem = mem

