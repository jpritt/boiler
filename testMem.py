#! /usr/bin/env python3

import resource
import subprocess
import threading
import sys

max_mem = 0
def finish():
    print('Maximum memory: %0.1f MB' % (max_mem / 1000.0))
    #exit()

sam_prefix = sys.argv[2]

def run():
    if sys.argv[1] == 'test':
        proc = subprocess.Popen(['./test.py'])
    elif sys.argv[1] == 'cram':
        print('Running CRAM')
        #proc = subprocess.Popen(['java', '-Xmx16g', '-jar', '/scratch0/langmead-fs1/shared/cramtools/cramtools-3.0.jar', 'cram', '-I', sam_prefix+'.bam', '-R', '/scratch0/langmead-fs1/user/jacob/genomes/Drosophila_melanogaster/Ensembl/BDGP5/Sequence/WholeGenomeFasta/genome.fa', '-O', 'compressed/compressed.cram'])
        proc = subprocess.Popen(['java', '-Xmx16g', '-jar', '/scratch0/langmead-fs1/shared/cramtools/cramtools-3.0.jar', 'cram', '-I', sam_prefix+'.bam', '-R', '/scratch0/langmead-fs1/shared/references/hg19/fasta/hg19.fa', '-O', 'compressed/compressed.cram'])
    elif sys.argv[1] == 'goby':
        print('Running Goby')
        proc = subprocess.Popen(['goby', '16g', 'sam-to-compact', '-i', sam_prefix+'.bam', '-o', 'compressed/compressed.goby'])
    else:
        print('Running Boiler')
        proc = subprocess.Popen(['/scratch0/langmead-fs1/shared/pypy3-2.4-linux_x86_64-portable/bin/pypy3', '/scratch0/langmead-fs1/user/jacob/compress-alignments/boiler.py', 'compress', '--frag-len-z-cutoff', '0.125', '--split-discordant', '--split-diff-strands', sam_prefix+'.sam', 'compressed/compressed.bin'])

    proc.wait()
    finish()

    return
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

while thread.isAlive():
    mem = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if mem > max_mem:
        max_mem = mem

