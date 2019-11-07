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
import csv
from itertools import islice


def main():
    # parser = OptionParser("usage: %prog [options]")
    # parser.add_option("-f", "--file1", dest="file1", help="fastq1 file")
    # parser.add_option("-e", "--file2", dest="file2", help="fastq2 file")
    # parser.add_option("-o", "--output_file", dest="output_file", help="output fastq merged file")
    # (options, args) = parser.parse_args()
    # file1 = options.file1
    # file2 = options.file2
    # output_file = options.output_file

    csvfile = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/RNASeq/HeLaRVB14RNAseq.csv"
    txt = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/RNASeq/human_MT_genes_mod.txt"
    count_engs(csvfile, txt)


def count_engs(csvfile, txt):
    summ = 0
    df = pd.DataFrame.from_csv(csvfile)

    with open(txt, 'r') as file2:
        fle2 = list(file2.readlines())
        print(fle2)
        df["flag"] = df['# Gene ID'].isin(fle2)
    print(df.head())
if __name__ == "__main__":
    main()