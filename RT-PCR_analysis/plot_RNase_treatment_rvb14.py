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
                    & (data.Sample != "STD4") & (data.Sample != "STD5") & (data.Sample != "STD6") & (data.Sample != "STD7")& (data.Sample != "NTC")]
    # data["Treatment"] = np.where(data["Sample"] == "treatment #1", "True", np.where(data["Sample"] == "treatment #2",
    #                                                                                 "True", np.where(data["Sample"] ==
    #                                                                                                  "treatment #3",
    #                                                                                                  "True", "False")))
    # pattern = re.compile("(\d\d|\d)", re.MULTILINE)
    # data["Sample"] = data["Sample"].str.extract(pattern, expand=True)
    # data["Sample"] = data["Sample"].astype(int)
    data["Copy_No"] = data["Starting Quantity (SQ)"]*rvb14_weight
    brust_size = 50
    limit = 7.56*10**5 * brust_size

    g1 = sns.lineplot(x="Sample", y="Cq", data=data, legend=False) #, hue="Treatment", hue_order=("False", "True")
    g1.set(xlabel="Sample", ylabel="Cq")
    g1.set(yticks=(range(0, 40, 2)))
    # g1.set_yscale("log")
    # g1.set_xticklabels(g1.get_xticklabels(), rotation=30)
    g1.set(title="RVB14")
    plt.xticks(rotation=45)
    plt.tight_layout()
    # plt.axhline(limit, color='red', ls='--')

    output_dir = (file_name.split("/")[0:-1])
    output_dir = '/'.join(output_dir)
    plt.savefig(output_dir + "/Cq.png", dpi=300)
    plt.close()

    g2 = sns.scatterplot(x="Sample", y="Log Starting Quantity", data=data, legend=False) # hue="Treatment", hue_order=("False", "True")
    g2.set(xlabel="Sample", ylabel="Log Starting Quantity (SQ) [ng]")
    g2.set(yticks=(range(-4, 8, 2)))
    # g1.set_yscale("log")
    # g1.set_xticklabels(g1.get_xticklabels(), rotation=30)
    g2.set(title="RVB14")
    plt.xticks(rotation=45)
    plt.tight_layout()
    # plt.axhline(limit, color='red', ls='--')

    output_dir = (file_name.split("/")[0:-1])
    output_dir = '/'.join(output_dir)
    plt.savefig(output_dir + "/SQ.png", dpi=300)
    plt.close()

    g3 = sns.scatterplot(x="Sample", y="Copy_No", data=data,
                         legend=False)  # hue="Treatment", hue_order=("False", "True")
    g3.set(xlabel="Sample", ylabel="Copy Number")
    # g3.set(yticks=(range(0, 14, 2)))
    g3.set_yscale("log")
    # g1.set_xticklabels(g1.get_xticklabels(), rotation=30)
    g3.set(title="RVB14")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.axhline(limit, color='red', ls='--')

    output_dir = (file_name.split("/")[0:-1])
    output_dir = '/'.join(output_dir)
    plt.savefig(output_dir + "/Copy_Number.png", dpi=300)

    return data

# def main(args):
def main():
    file_name = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/RV/RVB14/RealTime Results xlsx/20200811_RVB14_Amicon_Ultra_Capsid_Free/Oded_2020-08-11 22-57-27_BR003827 -  Quantification Cq Results_0.csv"
    data = plot_generation_time(file_name)
    print(data)

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