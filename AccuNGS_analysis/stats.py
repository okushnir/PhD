#! /Users/odedku/Anaconda3/python
"""
@Author: odedkushnir
"""
import datetime
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from natsort import index_natsorted
import glob
import re


def sts(file_path, new_file, virus, skiprows=6):
    with open(file_path, "r") as stsfile, open(new_file, 'w') as outfile:
        for line in stsfile:
            new_line = line.replace(": Kruskal-Wallis paired samples with Bonferroni correction", "").replace("P_val:", "").replace("Stat=", ", ")
            outfile.write(new_line)
    outfile.close()
    sts_df = pd.read_csv(new_file, skiprows=skiprows, header=None)
    sts_df = sts_df.rename({0: "Passage_Mutation vs. Passage_Mutation", 1: "p-value", 2: "stat"}, axis="columns")
    sts_df["Significant"] = np.where((sts_df["p-value"] >= 0.01) & (sts_df["p-value"] <= 0.05), "*",
                                     np.where((sts_df["p-value"] >= 0.001) & (sts_df["p-value"] <= 0.01), "**",
                                              np.where((sts_df["p-value"] >= 0.0001) & (sts_df["p-value"] <= 0.001), "***",
                                                       np.where(sts_df["p-value"] <= 0.0001, "****", "ns"))))
    sts_df["Test"] = "Kruskal-Wallis with Bonferroni correction"
    sts_df["Virus"] = virus
    sts_df.drop("stat", axis=1, inplace=True)
    sts_df = sts_df[["Virus", "Passage_Mutation vs. Passage_Mutation", "p-value", "Significant", "Test"]]
    sts_df = sts_df.sort_values(by="Passage_Mutation vs. Passage_Mutation",  key=lambda x: np.argsort(index_natsorted(sts_df["Passage_Mutation vs. Passage_Mutation"])))
    return sts_df

def iter_stat(inputdir_dir, rv_sheet, sheet):
    sts_all_lst = []
    for key, value in inputdir_dir.items():
        if key == "RV":
            data_sheet = rv_sheet
        else:
            data_sheet = sheet
        file_path = value + data_sheet + ".csv"
        new_file = value + data_sheet + "_mod.csv"
        sts_df = sts(file_path, new_file, key)
        sts_all_lst.append(sts_df)
    sts_all_df = pd.concat(sts_all_lst)
    return sts_all_df


def main():
    input_dir = "C:/Users/odedku/PhD_Projects/After_review/"
    inputdir_dir = {"RV": input_dir + "AccuNGS/RV/passages/20220314_inosine_predict_context/",
                    "CV": input_dir + "AccuNGS/CVB3/20220314_plots/", "OPV2": input_dir + "Cirseq/PV/OPV/20220314_plots/",
                    "PV1": input_dir + "Cirseq/PV/Mahoney/20220314_plots/"}
    sts_all_df = iter_stat(inputdir_dir, rv_sheet="sts_trans2", sheet="sts_trans")
    sts_all_df.to_csv("C:/Users/odedku/PhD_Projects/After_review/Stats/sts_all_trans_versions.csv")

    sts_all_df = iter_stat(inputdir_dir, rv_sheet="sts2", sheet="sts")
    sts_all_df.to_csv("C:/Users/odedku/PhD_Projects/After_review/Stats/sts_all_transitions.csv")

    sts_all_df = iter_stat(inputdir_dir, rv_sheet="sts_adar_2", sheet="sts_adar")
    sts_all_df.to_csv("C:/Users/odedku/PhD_Projects/After_review/Stats/sts_adar_all.csv")



if __name__ == "__main__":
    main()