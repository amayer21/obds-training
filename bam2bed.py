#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thursday Feb 01 2018 11:22
Exercise: write script to convert BAM into BED using Pysam,
starting from sam2bed script we wrote (@improved by Andreas)


@author: group, edited by Andreas
"""

import logging as L
import sys
import pysam


def main():

    L.basicConfig(filename='mylogfile.log', level=L.DEBUG)

    bamfile = pysam.AlignmentFile("BEL033_1000.bam", "rb")
    outf = open("BEL033_1000_v2.bed", 'w')

    pair_count = 0
    sum_intervals = 0

    initial_count = bamfile.count()

    for aln in bamfile.fetch():
        if aln.is_paired and aln.is_read1:
            pair_count += 1
            fgt_end = str(aln.reference_start+aln.template_length)
            sum_intervals += aln.template_length
            outf.write(aln.reference_name + '\t' + str(aln.reference_start)
                       + '\t' + fgt_end + '\t' + aln.query_name + '\t.\t.\n')

    bamfile.close()
    outf.close()

    av_fgt_size = sum_intervals / pair_count
    L.info("Initial count is {}".format(initial_count))
    L.info("Number of paired reads is {}".format(pair_count))
    L.info("Average fragment size is {:4.2f}".format(av_fgt_size))


if __name__ == "__main__":
    main()
