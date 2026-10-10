### removes linker sequence from Illumina sequencing results written by Jean-benoit Lalanne ###



import numpy as np

import os, subprocess

import sys



def removeLinker(fqIN,matchThres,linker='CTGTAGGCACCATCAAT'):

    inFile = open(fqIN, 'r')

    outFile = open(fqIN[:-6] + '-stripped_v2.fastq', 'w')

    dumpFile = open(fqIN[:-6] + '-notpassed_v2.fastq', 'w')

    line = inFile.readline()

    processed = 0

    passed = 0

    noinsert = 0

    while line != '':

        seq = inFile.readline()[:-1]

        line3 = inFile.readline()

        line4 = inFile.readline()

        length = findLinker(seq,matchThres,linker)

        if length > 3:

            outFile.write(line)

            outFile.write(seq[:length] + '\n')

            outFile.write(line3)

            outFile.write(line4[:length] + '\n')

            passed += 1

        elif length == -1:

            dumpFile.write(line)

            dumpFile.write(seq + '\n')

            dumpFile.write(line3)

            dumpFile.write(line4)

        else:

            noinsert += 1

        line = inFile.readline()

        processed += 1

    inFile.close()

    outFile.close()

    dumpFile.close()

    print 'File: ' + fqIN

    print 'Total  reads: '+str(processed)+'. Total  passed: '+str(passed)+'.'

    print 'Percent passed: '+str(passed*100/processed)+'%'

    print 'Percent no insert: '+str(noinsert*100/processed)+'%'

    return fqIN[:-4] + '-stripped.fq'



def findLinker(seq, matchThres,linker):

    exactmatch = seq.find(linker)

    if exactmatch + 1:

        return exactmatch

    else:

        s, l = [len(seq), len(linker)]

        maxpercent, maxindex = [0, 0]

        for i in range(s-1): #added -1 on 11/2/11

            percentmatch = sum([seq[i+j]==linker[j] for j in range(min(l, s-i))])*100/min(l, s-i)

            if percentmatch > maxpercent:

                maxpercent = percentmatch

                maxindex = i

        if maxpercent > matchThres:

            return maxindex

        else:

            return -1

        



if __name__ == '__main__':

    linker='CTGTAGGCACCATCAAT'

    fqIN = sys.argv[1]

    matchThres = 80

    fqOUT = removeLinker(fqIN, matchThres, linker)




