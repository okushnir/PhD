#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from AccuNGS_analysis import add_Protein_to_pd_df
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.colors as color

# sns.set(font_scale=2)
sns.set_style("ticks")
# sns.set_context("poster")
sns.despine()

def find_haplotype(pdser):
    global result
    for x in pdser:
        y = x.split(", ")
        if len(y)>2:
            continue
        result = pdser.str.findall(x)
    return result


def compare_2_strings(string1, string2):
    listA = string1.split(", ")
    listB = string2.split(", ")
    result = list(set(listA) & set(listB))
    return result

def find_max_per_column(list1, list2):
    L = [list1, list2]
    len_l = 0
    max = 0
    for l in L:
        len_l = len(l)
        if len_l > max:
            max = len_l
            max_list = l
    return max_list


def main():
    # input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages"
    input_dir = "/Users/odedkushnir/PhD_Projects/After_review/AccuNGS/RV/passages"
    output_dir = input_dir + "/20220101_co_occur_plots"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)
    replica_lst = ["1", "2", "3"]
    data_replica1 = pd.read_pickle(output_dir + "/all_co_occur_raw_replica1.pkl")
    data_replica2 = pd.read_pickle(output_dir + "/all_co_occur_raw_replica2.pkl")
    data_replica3 = pd.read_pickle(output_dir + "/all_co_occur_raw_replica3.pkl")
    data_all = pd.concat([data_replica1, data_replica2, data_replica3])

    data = data_all[data_all["Stretch_new"] == True]
    # data = data[data["label"] != "RNA Control\nPrimer ID"]
    data = data[data["numric_ref"] == 1]
    data["replica"] = data["replica"].astype(int)
    data = data.rename(columns={"replica": "Replica", "Pos": "Position"})
    data.to_csv(output_dir + "/data_all.csv")

    """removed by hand the groups with primerID Control in the Stretches 20220101"""

    data = pd.read_csv(output_dir + "/data_all.csv")
    data["Stretch"] = data["Stretch"].apply(lambda x: x.split("s")[1]).astype(int)
    data = data.rename(columns={"Position": "Genome Position"})
    kde_plot_stretch = sns.kdeplot(data=data, x="Genome Position", hue="Replica", bw_adjust=0.2, palette="Dark2")
    # kde_plot_stretch.set_title("ADAR-like stretches density across 12 passages of RV\nin 3 biological replicates")
    # kde_plot_stretch.set_yscale("log")
    # kde_plot_stretch.set_ylim(10**-5, 10**-1)
    plt.savefig(output_dir + "/kde_plot_withot_title.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()