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
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages"
    output_dir = input_dir + "/20201228_co_occur_plots"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)
    # replica_lst = ["1", "2", "3"]
    data_replica1 = pd.read_pickle(output_dir + "/all_co_occur_raw_replica1.pkl")
    data_replica2 = pd.read_pickle(output_dir + "/all_co_occur_raw_replica2.pkl")
    data_replica3 = pd.read_pickle(output_dir + "/all_co_occur_raw_replica3.pkl")
    data_all = pd.concat([data_replica1, data_replica2, data_replica3])

    data = data_all[data_all["Stretch_new"] == True]
    data = data[data["numric_ref"] == 1]
    data["replica"] = data["replica"].astype(int)
    data = data.rename(columns={"replica":"Replica", "Pos":"Position"})
    kde_plot_stretch = sns.kdeplot(data=data, x="Position", hue="Replica", bw_adjust=0.2, palette="Dark2")
    # kde_plot_stretch.set_title("ADAR-like stretches density across 12 passages of RV\nin 3 biological replicates")
    plt.savefig(output_dir + "/kde_plot_withot_title.png", dpi=300)
    plt.close()


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