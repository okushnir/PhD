#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import sys, argparse
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

def freqs_position_plot(file_name):
    # kjkdd

def main(args):
    # khjdkj

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("virus", type=str, help="name of the virus, RVB14")
    parser.add_argument("passages", type=str, help="from which passages, p0-p12")
    parser.add_argument("without", type=int, help="Exclude passage no.")
    parser.add_argument("input_dir", type=str, help="the path to the directory that contains data_mutation.csv")
    parser.add_argument("quality", type=str, help="what is the prefix for the data_mutation.csv file; quality of the pipline ; for example: q38")
    parser.add_argument("mutation_type", type=str, help="all/syn")
    args = parser.parse_args(sys.argv[1:])
    main(args)