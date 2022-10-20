#! /Users/odedku/Anaconda3/python
"""
@Author: odedkushnir
"""
import pandas as pd
import numpy as np
from natsort import index_natsorted


def sts(file_path, new_file, virus, skiprows=7):
    with open(file_path, "r") as stsfile, open(new_file, 'w') as outfile:
        for line in stsfile:
            new_line = line.replace(": Mann-Whitney-Wilcoxon test greater with Bonferroni correction", "").replace("P_val:", "").replace("U_stat=", ", ")
            outfile.write(new_line)
    outfile.close()
    sts_df = pd.read_csv(new_file, skiprows=skiprows, header=None)
    sts_df = sts_df.rename({0: "Passages and mutations compared", 1: "p-value", 2: "stat"}, axis="columns")
    sts_df["Significant"] = np.where((sts_df["p-value"] >= 0.01) & (sts_df["p-value"] <= 0.05), "*",
                                     np.where((sts_df["p-value"] >= 0.001) & (sts_df["p-value"] <= 0.01), "**",
                                              np.where((sts_df["p-value"] >= 0.0001) & (sts_df["p-value"] <= 0.001), "***",
                                                       np.where(sts_df["p-value"] <= 0.0001, "****", "ns"))))
    sts_df["Test"] = "Mann-Whitney-Wilcoxon test greater with Bonferroni correction"
    sts_df["Virus"] = virus
    sts_df.drop("stat", axis=1, inplace=True)
    sts_df = sts_df[["Virus", "Passages and mutations compared", "p-value", "Significant", "Test"]]
    sts_df = sts_df.sort_values(by="Passages and mutations compared",  key=lambda x: np.argsort(index_natsorted(sts_df["Passages and mutations compared"])))
    return sts_df


def iter_stat(inputdir_dir, sheet, rv_sheet=None):
    sts_all_lst = []
    for key, value in inputdir_dir.items():
        if key == "RV":
            data_sheet = rv_sheet
        elif key == "Non-Synonymous":
            data_sheet = rv_sheet
        else:
            data_sheet = sheet
        file_path = value + data_sheet + ".csv"
        new_file = value + data_sheet + "_mod.csv"
        sts_df = sts(file_path, new_file, key)
        sts_all_lst.append(sts_df)
    sts_all_df = pd.concat(sts_all_lst)
    return sts_all_df


def sts_fits(file_path, new_file, skiprows=6):
    with open(file_path, "r") as stsfile, open(new_file, 'w') as outfile:
        for line in stsfile:
            new_line = line.replace(": Welch's t-test independent samples with Bonferroni correction", "").replace("P_val:", "").replace("t=", ", ")
            outfile.write(new_line)
    outfile.close()
    sts_df = pd.read_csv(new_file, skiprows=skiprows, header=None)
    sts_df = sts_df.rename({0: "Types of mutations compared", 1: "p-value", 2: "stat"}, axis="columns")
    # sts_df["Significant"] = np.where((sts_df["p-value"] >= 0.01) & (sts_df["p-value"] <= 0.05), "*",
    #                                  np.where((sts_df["p-value"] >= 0.001) & (sts_df["p-value"] <= 0.01), "**",
    #                                           np.where((sts_df["p-value"] >= 0.0001) & (sts_df["p-value"] <= 0.001), "***",
    #                                                    np.where(sts_df["p-value"] <= 0.0001, "****", "ns"))))
    sts_df["Test"] = "Welch's t-test independent samples with Bonferroni correction"
    sts_df.drop("stat", axis=1, inplace=True)
    sts_df = sts_df[["Types of mutations compared", "p-value", "Test"]]
    sts_df = sts_df.sort_values(by="Types of mutations compared",  key=lambda x: np.argsort(index_natsorted(sts_df["Types of mutations compared"])))
    return sts_df

def main():
    # input_dir = "C:/Users/odedku/PhD_Projects/After_review/"
    # inputdir_dir = {"RV": input_dir + "AccuNGS/RV/passages/20221018_inosine_predict_context/",
    #                 "CV": input_dir + "AccuNGS/CVB3/20221018_plots/", "OPV2": input_dir + "Cirseq/PV/OPV/20221018_plots/",
    #                 "PV1": input_dir + "Cirseq/PV/Mahoney/20221018_plots/"}
    # sts_all_df = iter_stat(inputdir_dir, rv_sheet="sts_trans2", sheet="sts_trans")
    # sts_all_df.to_csv("C:/Users/odedku/PhD_Projects/After_review/Stats/sts_all_trans_versions.csv")

    # sts_all_df = iter_stat(inputdir_dir, rv_sheet="sts2", sheet="sts")
    # sts_all_df.to_csv("C:/Users/odedku/PhD_Projects/After_review/Stats/sts_all_transitions.csv")
    #
    # sts_all_df = iter_stat(inputdir_dir, rv_sheet="sts_adar_2", sheet="sts_adar")
    # sts_all_df.to_csv("C:/Users/odedku/PhD_Projects/After_review/Stats/sts_adar_all.csv")

    # capsid_dir = {"Synonymous": input_dir + "AccuNGS/RV/capsid/20221013_inosine_predict_context/",
    #               "Non-Synonymous": input_dir + "AccuNGS/RV/capsid/20221013_inosine_predict_context/"}
    # sts_all_df = iter_stat(capsid_dir, sheet="sts_adar_synon", rv_sheet="sts_adar_nonsynon")
    # sts_all_df.to_csv("C:/Users/odedku/PhD_Projects/After_review/Stats/sts_capsid_adar_all.csv")

    input_dir = "D:/My Drive/Studies/PhD/Thesis/FinalVer/AfterReview/20221018_fits_syn_plots"
    sts_file = input_dir + "/sts.csv"
    new_file = input_dir + "/sts_edited.csv"
    sts_fits_df = sts_fits(sts_file, new_file)
    sts_fits_df.to_csv("D:/My Drive/Studies/PhD/Thesis/FinalVer/AfterReview/20221018_fits_syn_plots/sts_fits_df.csv")


if __name__ == "__main__":
    main()