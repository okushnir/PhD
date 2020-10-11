#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import sys, argparse
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

def load_file(path, label):
    df = pd.read_csv(path, sep="\t", names=["start", "end", "fisher_stat", "pval", "variant_freq"])
    df["sample"] = label
    df = df[df["pval"] < 1]
    return df

def main():
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/190627_RV_CV/merged/RVB14"
    all_p2 = load_file(input_dir + "/RVB14_p2/20191029_q38/accungs_associations/all.txt", "02")
    all_p5 = load_file(input_dir + "/RVB14_p5/20191029_q38/accungs_associations/all.txt", "05")
    all_p8 = load_file(input_dir + "/RVB14_p8/20191029_q38/accungs_associations/all.txt","08")
    all_p10 = load_file(input_dir + "/RVB14_p10/20191029_q38/accungs_associations/all.txt","10")
    all_p12 = load_file(input_dir + "/RVB14_p12/20191029_q38/accungs_associations/all.txt", "12")
    all_df = pd.concat([all_p2, all_p5, all_p8, all_p10, all_p12])
    lst = [6019, 6043, 6049, 6067, 6076, 6078, 6079, 6082, 6097]
    all_df["filter"] = all_df["start"].isin(lst)
    all_df = all_df[all_df["filter"] == True]

    all_df["start"] = all_df["start"].astype(str)
    all_df["end"] = all_df["end"].astype(str)
    # all_df["paires"] = all_df["start"] + ", " + all_df["end"]
    all_df = (all_df.groupby(["sample", "start"]).agg(lambda x: x.mean() if np.issubdtype(x.dtype, np.number)
    else ', '.join(x))).reset_index()
    print(all_df.to_string())
    all_df["start"] = all_df["start"].astype(int)
    all_df["end"] = all_df["end"].astype(int)



    # all_df.to_csv(input_dir + "/all_co_occur_grouped_paires.csv", sep=",", encoding='utf-8')

    # g2 = sns.swarmplot(x="sample", y="variant_freq", data=all_df, hue="paires")
    # g2.set(yscale="log")
    # plt.show()
    # g2.savefig(output_dir + "/co_occur_protein.png", dpi=300)
    # plt.close()

if __name__ == "__main__":
    main()