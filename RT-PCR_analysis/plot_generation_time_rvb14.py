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
    rvb14_weight = (6.022*10**23)/(7212*1*10**9*660)
    data = data.loc[(data.Sample != "IVT") & (data.Sample != "STD1") & (data.Sample != "STD2") & (data.Sample != "STD3")
                    & (data.Sample != "STD4") & (data.Sample != "STD5") & (data.Sample != "NTC") & (data.Sample != "13hr")
                    & (data.Sample != "7hr")& (data.Sample != "11hr")]
    pattern = re.compile("(\d\d|\d)", re.MULTILINE)
    data["Sample"] = data["Sample"].str.extract(pattern, expand=True)
    data["Sample"] = data["Sample"].astype(int)
    data["Copy_No"] = data["Starting Quantity (SQ)"]*rvb14_weight
    brust_size = 40
    limit = 7.56*10**5 * brust_size
    g1 = sns.lineplot(x="Sample", y="Copy_No", data=data)
    g1.set(xlabel="Time[hr]", ylabel="vRNA Copy number")
    g1.set_yscale("log")
    g1.set(xticks=(range(1, 15)))
    g1.set(title="RVB14")
    plt.axhline(limit, color='red', ls='--')

    output_dir = (file_name.split("/")[0:-1])
    output_dir = '/'.join(output_dir)
    plt.savefig(output_dir + "/Generation_Time_Curve.png", dpi=300)
    return data

# def main(args):
def main():
    file_name = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/RV/RVB14/RealTime Results xlsx/20200303_RVB14_Time/Oded_2020-03-03 13-51-03_BR003827 -  Quantification Cq Results_0.csv"
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