#!/powerapps/share/python-anaconda-3.2019.7/bin/python

"""
@Author: odedkushnir

"""

import sys, argparse
import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re
import glob
import Bio.Seq as Seq
import matplotlib.gridspec as gridspec
import time

def plot_generation_time(file_name):
    """
    :param file_name:
    :return: plot and data frame
    """
    data = pd.read_csv(file_name)
    data = data.loc[data.Sample != 14]
    g1 = sns.lineplot(x="Sample", y="Copy No.", data=data)
    g1.set(xlabel="Time[hr]", ylabel="Copy number of RNA")
    g1.set_yscale("log")
    output_dir = (file_name.split("/")[0:-1])
    output_dir = '/'.join(output_dir)
    plt.savefig(output_dir + "/Copy_number.png", dpi=300)
    return data

# def main(args):
def main():
    file_name = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/RVB14/RealTime Results xlsx/admin_2016-06-14 13-05-44_BR003827 -  Quantification Cq Results.csv"
    plot_generation_time(file_name)

if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    # parser.add_argument("virus", type=str, help="name of the virus, RVB14")
    # parser.add_argument("passages", type=str, help="from which passages, p0-p12")
    # parser.add_argument("without", type=int, help="Exclude passage no.")
    # parser.add_argument("input_dir", type=str, help="the path to the directory that contains data_mutation.csv")
    # parser.add_argument("quality", type=str, help="what is the prefix for the data_mutation.csv file; quality of the pipline ; for example: q38")
    # parser.add_argument("mutation_type", type=str, help="all/syn")
    # args = parser.parse_args(sys.argv[1:])
    # main(args)
    main()