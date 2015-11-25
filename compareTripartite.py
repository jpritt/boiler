#! /usr/bin/env python
import sys
import logging
from transcript import Transcript

'''
Compare a Cufflinks GTF file to the .pro output by flux
'''


def read_xscripts(gtf_fn):
    transcriptsComp = dict()
    with open(gtf_fn, 'r') as tsv:
        for line in tsv:
            row = line.strip().split('\t')
            if len(row) < 5:
                continue

            covIndex = row[8].find('cov')
            covStart = row[8].find('"', covIndex) + 1
            covEnd = row[8].find('"', covStart)
            cov = float(row[8][covStart:covEnd])

            transcriptIdIndex = row[8].find('transcript_id')
            transcriptIdStart = row[8].find('"', transcriptIdIndex) + 1
            transcriptIdEnd = row[8].find('"', transcriptIdStart)
            transcriptId = row[8][transcriptIdStart:transcriptIdEnd]

            if row[2] == 'transcript':
                transcriptsComp[transcriptId] = Transcript(row[0], int(row[3]), int(row[4]), cov, transcriptId)
            elif row[2] == 'exon':
                transcriptsComp[transcriptId].exons.append((int(row[3]), int(row[4])))
    return transcriptsComp


def compareTripartite(proFile, truthGTF, compGTF1, compGTF2, strict):
    # Read reference transcripts
    # file 1 is a .pro file output by flux
    transcriptsTruth = parsePro(proFile)
    logging.info('Parsed %d transcripts from simulation .pro file' % len(transcriptsTruth))

    nexons = 0
    with open(truthGTF, 'r') as tsv:
        for line in tsv:
            row = line.strip().split('\t')
            if len(row) < 5:
                continue

            transcriptIdIndex = row[8].find('transcript_id')
            transcriptIdStart = row[8].find('"', transcriptIdIndex) + 1
            transcriptIdEnd = row[8].find('"', transcriptIdStart)
            transcriptId = row[8][transcriptIdStart : transcriptIdEnd]

            #if row[2] == 'transcript':
            #    transcriptsTruth[transcriptId] = Transcript(row[0], int(row[3]), int(row[4]), transcriptCovs[transcriptId])

            #if row[1] == 'protein_coding' and row[2] == 'exon' and transcriptId in transcriptsTruth:
            if row[2] == 'exon' and transcriptId in transcriptsTruth:
                transcriptsTruth[transcriptId].exons.append( (int(row[3]), int(row[4])) )
                nexons += 1


    for transcriptId, xscript in transcriptsTruth.items():
        assert len(xscript.exons) > 0
    transcriptsTruth = [v for v in transcriptsTruth.values()]

    #transcriptsTruth = []
    #for v in transcriptsTruth.values():
    #    if len(v.exons) > 0:
    #        transcriptsTruth.append(v)

    logging.info('Parsed %d protein-coding exons from simulation .gtf file' % nexons)

    # Read transcripts from first GTF
    transcriptsCompA = [v for v in read_xscripts(compGTF1).values()]
    logging.info('Parsed %d transcripts from first .gtf file' % len(transcriptsCompA))

    # Read transcripts from second GTF
    transcriptsCompB = [v for v in read_xscripts(compGTF2).values()]
    logging.info('Parsed %d transcripts from second .gtf file' % len(transcriptsCompB))

    connectionsA, connectionsB = buildGraph(transcriptsTruth, transcriptsCompA, transcriptsCompB)
    nconnected_a = sum(map(lambda x: len(x) > 0, connectionsA))
    nconnected_b = sum(map(lambda x: len(x) > 0, connectionsB))
    logging.info('%d of %d transcripts in set A have near neighbor in true transcript set' %
                 (nconnected_a, len(transcriptsCompA)))
    logging.info('%d of %d transcripts in set B have near neighbor in true transcript set' %
                 (nconnected_b, len(transcriptsCompB)))
    if nconnected_a == 0:
        raise RuntimeError('0 transcripts from first GTF file are similar to any true transcript!')
    if nconnected_b == 0:
        raise RuntimeError('0 transcripts from second GTF file are similar to any true transcript!')

    matches = 0
    score = 0.0
    threshold = 10

    xs = []
    ys = []

    for i in range(len(connectionsA)):
        if strict:
            # strict
            if len(connectionsA[i]) == 1 and len(connectionsB[i]) == 1:
                xs += [transcriptsCompA[connectionsA[i][0]].cov]
                ys += [transcriptsCompB[connectionsB[i][0]].cov]

                matches += 1
                score += transcriptsCompA[connectionsA[i][0]].scoreTranscript(transcriptsCompB[connectionsB[i][0]], threshold)
        
        else:
        # loose
            if len(connectionsA[i]) > 0 and len(connectionsB[i]) > 0:
                matches += 1
                bestA = None
                bestScoreA = 0
                for j in connectionsA[i]:
                    t = transcriptsCompA[j]
                    currScore = t.scoreTranscript(transcriptsTruth[i], threshold)
                    if currScore > bestScoreA:
                        bestScoreA = currScore
                        bestA = t

                bestB = None
                bestScoreB = 0
                for j in connectionsB[i]:
                    t = transcriptsCompB[j]
                    currScore = t.scoreTranscript(transcriptsTruth[i], threshold)
                    if currScore > bestScoreB:
                        bestScoreB = currScore
                        bestB = t

                xs += [bestA.cov]
                ys += [bestB.cov]

                score += bestA.scoreTranscript(bestB, threshold)
        
    print('%d / %d  = %0.4f transcripts from file 1' % (matches, len(transcriptsCompA), float(matches)/len(transcriptsCompA)))
    print('%d / %d  = %0.4f transcripts from file 2' % (matches, len(transcriptsCompB), float(matches)/len(transcriptsCompB)))
    print('Score: %f' % (score / matches))

    #plt.scatter(xs, ys)
    #plt.show()


def parsePro(filename):
    ''' Return a dictionary with transcript id (e.g. 0300689) pointing to coverage level
    '''
    threshold = 0.00005

    transcripts = dict()
    with open(filename, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            if len(row) < 8:
                continue

            tag = row[1]
            sep1 = row[0].find(':')
            sep2 = row[0].find('-', sep1)
            sep3 = row[0].find('W', sep2)
            chrom = row[0][:sep1]
            start = int(row[0][sep1+1:sep2])
            end = int(row[0][sep2+1:sep3])
            fraction = float(row[8])

            if fraction > threshold:
                transcripts[tag] = Transcript(chrom, start, end, fraction, tag)
    return transcripts


def make_connections(xscript_pred, xscript_true, threshold=10):
    """ Find true transcript most similar to each transcript in xscript_pred.
        Return a mapping.  If similarity doesn't exceed threshold, transcripts
        aren't considered similar. """
    connect = []
    for i in range(len(xscript_true)):
        connect += [[]]

    for i, t in enumerate(xscript_pred):
        best_score, best_xscript = 0, 0
        for j, truet in enumerate(xscript_true):
            score = t.scoreTranscript(truet, threshold)
            if score > best_score:
                best_score, best_xscript = score, j

        if best_score > 0:
            connect[best_xscript] += [i]
    return connect


def buildGraph(transcriptsTrue, transcriptsPredictedA, transcriptsPredictedB, threshold=10):
    """ Build mappings between true/A and true/B """
    refConnectionsA = make_connections(transcriptsPredictedA, transcriptsTrue, threshold)
    refConnectionsB = make_connections(transcriptsPredictedB, transcriptsTrue, threshold)
    return refConnectionsA, refConnectionsB


def scoreTranscripts(exons1, exons2, threshold):
    thisMatches = []
    scores = []
    for i in range(len(exons1)):
        e1 = exons1[i]
        closest = 0
        closestScore = 0

        for j in range(len(exons2)):
            e2 = exons2[j]
            currScore = scoreExons(e1, e2, threshold)

            if currScore > closestScore:
                closestScore = currScore
                closest = j

        thisMatches.append(closest)
        scores.append(closestScore)

    otherMatches = []
    for i in range(len(exons2)):
        e1 = exons2[i]
        closest = 0
        closestScore = 0

        for j in range(len(exons1)):
            e2 = exons1[j]
            currScore = scoreExons(e1, e2, threshold)

            if currScore > closestScore:
                closestScore = currScore
                closest = j

        otherMatches.append(closest)

    totalScore = 0.0
    count = len(thisMatches)
    for i in range(len(thisMatches)):
        if otherMatches[thisMatches[i]] == i:
            totalScore += scores[i]
    for i in range(len(otherMatches)):
        if not thisMatches[otherMatches[i]] == i:
            count += 1

    return totalScore / float(count)

def scoreExons(exon1, exon2, threshold):
    ''' Returns a value between 0 and 1, where 1 indicates that the 2 exons are a perfect match and 0 indicates no match
    '''

    dx = min(abs(exon1[0] - exon2[0]), threshold)
    dy = min(abs(exon1[1] - exon2[1]), threshold)

    return 1 - dx / (2.0*threshold) - dy / (2.0*threshold)


logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', datefmt='%m/%d/%y-%H:%M:%S', level=logging.INFO)

ref_pro = sys.argv[1]
ref_gtf = sys.argv[2]
alignments1 = sys.argv[3]
alignments2 = sys.argv[4]

strict = True
if int(sys.argv[5]) == 0:
    strict = False
    logging.info('Comparing tripartite loose')
else:
    logging.info('Comparing tripartite strict')
compareTripartite(ref_pro, ref_gtf, alignments1, alignments2, strict)
