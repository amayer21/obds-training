#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 17:38:37 2018
@author: group, edited by Andreas

Extract data from a SAM file and generate a BED file 
(using only basic Python commands, no SAM tools or pysam)

- read SAM files from std input
- separate columns
- does it map as a pair? it has an '=' when it is the same chromosome
- select only Fwd reads (col9 > 0)
- find the start position of interval = col4 (left most nucleotide of Fwd read)
- calculate the end position of interval
    start of interval (col4) + length (col9)
    NB: col8 = left most nucl of Rev read, ie end position of Rev read
- as ID for the interval in BED file, we'll use the name of the read (col1)
    
Added by Andreas:

- Define a function main that contains the code and is run at the end.
- Introduce log file properly (to output read counts). 
    I think it's going to append that file each time we run the script 
    => keep trace of all runs in 1 file.

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
    seq_name = ''
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
            read = line.split('\t')
            TotalCount += 1
            if read[6] == '=':
                RawPairedReads += 1
                if int(read[8]) > 0:
                    # take only the reads in on the positive strand
                    PairsCount += 1
                    chrom = read[2]
                    start = int(read[3])
                    end = int(start)+int(read[8])
                    seq_name = read[0]
                    Score = '.'
                    output.write(chrom + '\t' + str(start) + '\t' + str(end)
                        + '\t' + seq_name + '\t' + Score + '\t' + Strand + '\n')
                    SumIntervals += int(read[8])

    AvgIntervals = SumIntervals/PairsCount

    L.info('The number of initial number of reads in the input file is:\t'
           + str(TotalCount))
    L.info('The number of paired-end reads is:\t'+str(PairsCount))
    L.info('The number of raw paired reads is:\t'+str(RawPairedReads))
    L.info('The average interval length is:\t'+str(AvgIntervals))


if __name__ == "__main__":
    main()
