#! /usr/local/python_anaconda/bin/python3.4

"""
@Author: odedkushnir

"""

import sys
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re
import glob
import Bio.Seq as Seq
import matplotlib.gridspec as gridspec
import time

def main():
    # parser = OptionParser("usage: %prog [options]")
    # parser.add_option("-f", "--file1", dest="file1", help="fastq1 file")
    # parser.add_option("-e", "--file2", dest="file2", help="fastq2 file")
    # parser.add_option("-o", "--output_file", dest="output_file", help="output fastq merged file")
    # (options, args) = parser.parse_args()
    # file1 = options.file1
    # file2 = options.file2
    # output_file = options.output_file

   #for Local

    suffix = "noGEL.freqs"
    freqs_file = "/Volumes/STERNADILABTEMP$/volume1/okushnir/AccuNGS/" + suffix

    parse_reads(freqs_file)

def parse_reads(freqs):

    ''' this method returns a vector of reads corresponding to genome positions.
input:
        freqs file
output:
        an integer vector containing for each position in the genome it's num of reads.
'''

    path = freqs
    df = pd.read_csv(path, sep='\t')

    # remove all duplicates from Pos except the first occurrence
    # remove all x.number duplicates
    df[["Pos"]] = df[["Pos"]].astype(int)
    df = df.drop_duplicates("Pos")

    pos = df["Pos"]  # a vector of all positions
    reads = df["Read_count"]
    med_reads = reads.median()
    print(med_reads)
    return pos, reads


if __name__ == "__main__":
    main()