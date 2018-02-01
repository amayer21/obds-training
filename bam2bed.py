#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 17:38:37 2018
This is our SAM_to_BED.py @improved by Andreas

Define a function main that contains the code and is run at the end.
Introduce log file properly. I think it's gone append that file each time
we run the script => keep trace of all runs in 1 file.

@author: group, edited by Andreas
"""

import logging as L
import sys


def main():

    L.basicConfig(filename='mylogfile.log', level=L.DEBUG)

    file = sys.stdin
    output = sys.stdout

    chrom = ''
    start = 0
    end = 0
    ID = ''
    Score = ''
    Strand = '+'
    TotalCount = 0
    PairsCount = 0
    RawPairedReads = 0
    SumIntervals = 0

    for line in file:
        if line.startswith('@'):
            pass
        else:
            name = line.split('\t')
            TotalCount += 1
            if name[6] == '=':
                RawPairedReads += 1
                if int(name[8]) > 0:
                    # take only the reads in on the positive strand
                    PairsCount += 1
                    chrom = name[2]
                    start = int(name[3])
                    end = int(start)+int(name[8])
                    ID = name[0]
                    Score = name[4]
                    output.write(chrom + '\t' + str(start) + '\t' + str(end)
                        + '\t' + ID + '\t' + Score + '\t' + Strand + '\n')
                    SumIntervals += int(name[8])

    AvgIntervals = SumIntervals/PairsCount

    L.info('The number of initial number of reads in the input file is:\t'
           + str(TotalCount))
    L.info('The number of paired-end reads is:\t'+str(PairsCount))
    L.info('The number of raw paired reads is:\t'+str(RawPairedReads))
    L.info('The average interval length is:\t'+str(AvgIntervals))


if __name__ == "__main__":
    main()
