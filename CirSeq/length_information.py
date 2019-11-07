#! /usr/local/python-anaconda-2.7//bin/python

import sys
sys.path.insert(0, '/sternadi/home/volume1/taliakustin/scripts/')
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import json, re, sys, getopt, math
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from scipy.stats import norm, beta
from optparse import OptionParser
from file_utilities import check_dirname
import re
import glob

sns.set_context("talk")


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--dir", dest="dir", help="dir of temp files of cirseq pipeline run")
    parser.add_option("-o", "--output", dest="output", help="output folder to save")

    (options, args) = parser.parse_args()
    in_dir = options.dir
    out_dir = options.output
    in_dir = check_dirname(in_dir)
    out_dir = check_dirname(out_dir)

    repeat_summmery = get_repeats_num(in_dir)
    make_graph(repeat_summmery, out_dir)



def get_repeats_num(in_dir):
    files = glob.glob(in_dir + "/*.fasta.blast.freqs.stats")

    repeat_summmery = {}

    for file in files:
        pattern = re.compile("(\d+\t{1}\d+\n{1})", re.MULTILINE)
        text = open(file, "rb").read()
        reads = pattern.findall(text)
        for r in reads:
            key = int(r.split("\t")[0])
            value = int(r.split("\t")[1].split("\n")[0])
            if key in repeat_summmery:
                repeat_summmery[key] += value
            else:
                repeat_summmery[key] = value

    return repeat_summmery


def make_graph(repeat_summmery, out_dir):
    x = repeat_summmery.keys()
    y = repeat_summmery.values()

    plt.bar(x,y, color="pink", align="center")

    ax = plt.gca()
    ax.set_xlabel("num of reads")
    ax.set_ylabel("count")
    ax.set_xticks(list(range(1,11)))
    plt.xlim(min(x)-0.5, max(x)+0.5)
    plt.yscale("log")
    sns.set_style("darkgrid")

    path = out_dir + "/repeat_num.png"
    plt.savefig(path, dpi=680)
    plt.close('all')











if __name__ == "__main__":
    main()