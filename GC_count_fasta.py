#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 14:17:29 2018
@author: amayer

Compute the %GC in a fasta file containing several sequences and record
the output in a tsv file (2 columns: seq name and %GC)

When line is a header (starts with >),
    we calculate the %GC of previous entry (except if first line), save it in
    a dictionary (for very big files may be worth writing directly in file)
    then re-initialise the counts and save the name of next read.

When line is sequence, we loop through counting function

"""

GCcount = 0
ATcount = 0
other = 0
GCpercent = 0
name = ''
savedGC = {}


def basecount(line, GCcount, ATcount, other):
    for base in line.rstrip().upper():
        if base == 'C' or base == 'G':
            GCcount = GCcount + 1
        elif base == 'A' or base == 'T':
            ATcount = ATcount+1
        else:
            other = other+1
    return GCcount, ATcount, other


with open("Homo_sapiens_CEBPB_sequence2.fa.txt", 'r') as f:
    for line in f:
        if line.startswith('>'):
            if name != '':        # ie not first line of file
                GCpercent = GCcount/(GCcount + ATcount)
                savedGC[name] = GCpercent
                print('The GC content of CEBPB gene is:', name, "%.2f"
                      % GCpercent)
            GCcount = 0
            ATcount = 0
            other = 0
            name = line.split()[0].lstrip('>')
            # to get the first element from line and remove the >
        else:
            GCcount, ATcount, other = basecount(line, GCcount, ATcount, other)

    GCpercent = GCcount/(GCcount + ATcount)
    print('The GC content of CEBPB gene is:', name, "%.2f" % GCpercent)
    savedGC[name] = GCpercent
    print(GCcount, ATcount, other, GCpercent)


with open('Homo_sapiens_CEBPB_GCcontent.tsv', 'w') as f:
    for key in savedGC:
        f.write(key + '\t' + str(savedGC[key]) + '\n')
