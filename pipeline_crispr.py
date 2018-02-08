#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 14:00:34 2018
@author: amayer

Exercise: pipeline process fastq file from CRISPR screen to generate read count
This is our first pipeline, so to keep it easy we use the trim function from 
Bowtie. Therefore we can process only the file with a specific stagger (GFP-1)
(To process all files, we'll need to use a regular expression to extract 
sgRNA instead)

Library file had to be indexed before to start 
(from command line: $ bowtie-build FASTA_GeckoV2_HGLib_B.fasta geckob)

- bowtie: trim 5' and 3' end and align it to GecKOv2A library
- the result of bowtie is piped to samtools command to convert into BAM file
- BAM file is sorted (may want to pipe that to previous command to avoid 
                      keeping intermediate files
- sorted BAM file is indexed
- statistics are extracted from BAM file into tsv file

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


# trim 5' and 3' end of reads to keep only the sgRNA sequence,
# then map reads to indexed library, then convert sam file into bam file
# also generate a .error file for the alignement
@transform("*.fastq.gz",
           suffix(".fastq.gz"),
           ".bam")
def run_bowtie(infile, outfile):
    job_threads = PARAMS['bowtie_threads']
    bowtie_options = PARAMS['bowtie_options']
    bowtie_index = PARAMS['bowtie_index']
    samtools_vie = PARAMS['samtools_view']
    statement = '''/usr/bin/time -o bowtie.time -v bowtie %(bowtie_options)s
                    %(bowtie_index)s %(infile)s
                    2> %(infile)s.error
                    | samtools view %(samtools_view)s > %(outfile)s '''

    P.run()


# Sort the BAM file (could have been piped to previous command)
@transform(run_bowtie,
           suffix(".bam"),
           "_sorted.bam")
def sort_bam(infile, outfile):

    statement = '''samtools sort %(infile)s > %(outfile)s '''

    P.run()


# Index the sorted BAM file and create index file (.bam.bai)
@transform(sort_bam,
           suffix(".bam"),
           ".bam.bai")
def index_bam(infile, outfile):

    statement = '''samtools index %(infile)s > %(outfile)s '''

    P.run()


# follows: wait for the index file to be created before to run next command
# then compute statistics about the bam file and export them to tsv file
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
