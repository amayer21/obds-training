#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 14:00:34 2018

@author: amayer

docs removed
"""
import sys
import os
from ruffus import transform, follows, suffix
import CGATPipelines.Pipeline as P

# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

# ---------------------------------------------------
# Specific pipeline tasks


@transform("*.fastq.gz",
           suffix(".fastq.gz"),
           ".bam")
def run_bowtie(infile, outfile):

    statement = '''bowtie --trim5 27 --trim3 27 --sam geckob %(infile)s
                    2> %(infile)s.error | samtools view -b > %(outfile)s '''

    P.run()


@transform(run_bowtie,
           suffix(".bam"),
           "_sorted.bam")
def sort_bam(infile, outfile):

    statement = '''samtools sort %(infile)s > %(outfile)s '''

    P.run()


@transform(sort_bam,
           suffix(".bam"),
           ".bam.bai")
def index_bam(infile, outfile):

    statement = '''samtools index %(infile)s > %(outfile)s '''

    P.run()


@follows(index_bam)
@transform(sort_bam,
           suffix(".bam"),
           ".tsv")
def count_hits(infile, outfile):

    statement = '''samtools idxstats %(infile)s > %(outfile)s '''

    P.run()

# ---------------------------------------------------
# Generic pipeline tasks


@follows(count_hits)
def full():
    pass


def main(argv):
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
