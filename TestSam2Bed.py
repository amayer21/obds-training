#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 17:45:17 2018
@author: Andreas

Set of test for the script sam2bed.py

We are going to test on 2 data sets:
    paired.sam is a sam file with about 2000 entries
    simple.sam is a sam file with only 2 entries (paired reads)
    
On each file we'll preform 3 tests:
        check if there is any bugs (report if sth prevent script to run)
        check that number lines in files = nb paired reads in original file
        check lines splits correctly
        
Before to run: create the bam files to test

"""

import unittest
import pysam
import os
from sam2bed import split_line


class TestSam2Bed(unittest.TestCase):

    filename = "paired.sam"

    def build_statement(self):
        return "python sam2bed.py < {} > output.bed".format(self.filename)

    def test_sam2bed_exits_without_error(self):

        # execute statement, 0 is success
        return_value = os.system(self.build_statement())
        self.assertEqual(return_value, 0)

    def test_number_of_lines_is_equal(self):

        # execute statement, 0 is success
        return_value = os.system(self.build_statement())
        self.assertEqual(return_value, 0)

        with pysam.AlignmentFile(self.filename) as inf:
            npaired_reads = len(set([x.query_name for x in inf if x.is_paired]))

        with open("output.bed") as inf:
            nlines = len([x for x in inf])
        self.assertEqual(npaired_reads, nlines)

    def test_line_splits_correctly(self):
        self.assertEqual(split_line("a\tb"), ["a", "b"])


class TestSam2BedSimple(TestSam2Bed):
    filename = "simple.sam"


if __name__ == "__main__":
    unittest.main()
