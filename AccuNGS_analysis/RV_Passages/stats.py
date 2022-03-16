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
from scipy import stats as sts
import glob
import re


def sts(file_path, new_file, skiprows=6):
    with open(file_path, "r") as stsfile, open(new_file, 'w') as outfile:
        for line in stsfile:
            new_line = line.replace(": Kruskal-Wallis paired samples with Bonferroni correction", "").replace("P_val:", "").replace("Stat=", ", ")
            print(new_line)
            outfile.write(new_line)
    outfile.close()
    sts_df = pd.read_csv(new_file, skiprows=skiprows, header=None)
    sts_df.rename({0: "Passage_Mutation vs. Passage_Mutation", 1: "p-value", 2: "stat"}, axis="columns")
    sts_df["test"] = "Kruskal-Wallis with Bonferroni correction"
    sts_df["significant"] = np.where((sts_df[1] >= 0.01) & (sts_df[1] <= 0.05), "*",
                                     np.where(sts_df[1] >= 0.001) & (sts_df[1] <= 0.01), "**",
                                     np.where(sts_df[1] >= 0.01) & (sts_df[1] <= 0.05), "***") #C2>=0.001,C2<=0.01),"**",IF(AND(C2>=0.0001,C2<=0.001),"***",IF(AND(C2<=0.0001),"****","ns"))))
    print(sts_df.to_string())

def main():
    input_dir = "C:/Users/odedku/PhD_Projects/After_review/AccuNGS/RV/passages/20220314_inosine_predict_context/"
    data_sheet = "sts_adar_2"
    file_path = input_dir + data_sheet + ".csv"
    new_file = input_dir + data_sheet + "_mod.csv"
    sts(file_path, new_file)

if __name__ == "__main__":
    main()