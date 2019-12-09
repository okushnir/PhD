#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import sys, argparse
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from AccuNGS_analysis import add_Protein_to_pd_df
from AccuNGS_analysis import new_analysis_RV
import numpy as np

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
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/190627_RV_CV/merged/RVB14"
    output_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/20191209_plots"
    data_p2 = pd.read_csv(input_dir + "/RVB14_p2/20191029_q38/co_occur_non_adar.csv")
    data_p5 = pd.read_csv(input_dir + "/RVB14_p5/20191029_q38/co_occur_non_adar.csv")
    data_p8 = pd.read_csv(input_dir + "/RVB14_p8/20191029_q38/co_occur_non_adar.csv")
    data_p10 = pd.read_csv(input_dir + "/RVB14_p10/20191029_q38/co_occur_non_adar.csv")
    data_p12 = pd.read_csv(input_dir + "/RVB14_p12/20191029_q38/co_occur_non_adar.csv")
    data_p0 = pd.read_csv(input_dir + "/RVB14_p0/20191029_q38/co_occur_non_adar.csv")
    df = pd.concat([data_p0, data_p2, data_p5, data_p8, data_p10, data_p12])
    # df = df.loc[df.Mutation == "C>U"]
    df.to_csv(input_dir + "/all_co_occur_nonadar_raw.csv", sep=",", encoding='utf-8')

    reg_lst = [629, 835, 1621, 2329, 3196, 3634, 3925, 4915, 5170, 5239, 5785, 7165]
    df = add_Protein_to_pd_df.add_Protein_to_pd_df(df, reg_lst)
    df["Pos"] = df["Pos"].astype(str)
    df_co_occur_new = (df.groupby(["label", "Stretch"]).agg(lambda x: x.mean() if np.issubdtype(x.dtype, np.number)
    else ', '.join(x))).reset_index()
    df_co_occur_new = df_co_occur_new[["label", "Stretch", "Pos", "Frequency"]]
    df_co_occur_new = df_co_occur_new.sort_values(by=["Pos", "label"])
    df_co_occur_new["filter"] = df_co_occur_new["Pos"].str.contains(",")
    df_co_occur_new = df_co_occur_new[df_co_occur_new["filter"] == True]
    df_co_occur_new["Next_1line_Pos"] = df_co_occur_new["Pos"].shift(periods=-1)
    df_co_occur_new["Next_2line_Pos"] = df_co_occur_new["Pos"].shift(periods=-2)
    df_co_occur_new["Next_3line_Pos"] = df_co_occur_new["Pos"].shift(periods=-3)
    df_co_occur_new["Next_4line_Pos"] = df_co_occur_new["Pos"].shift(periods=-4)
    df_co_occur_new["Next_5line_Pos"] = df_co_occur_new["Pos"].shift(periods=-5)
    df_co_occur_new = df_co_occur_new.fillna("0, ")

    # finds the intersection group
    df_co_occur_new["Next_1line_Pos_filter"] = df_co_occur_new.apply(
        lambda x: compare_2_strings(x["Pos"], x["Next_1line_Pos"]), axis=1)
    df_co_occur_new["Next_2line_Pos_filter"] = df_co_occur_new.apply(
        lambda x: compare_2_strings(x["Pos"], x["Next_2line_Pos"]), axis=1)
    df_co_occur_new["Next_3line_Pos_filter"] = df_co_occur_new.apply(
        lambda x: compare_2_strings(x["Pos"], x["Next_3line_Pos"]), axis=1)
    df_co_occur_new["Next_4line_Pos_filter"] = df_co_occur_new.apply(
        lambda x: compare_2_strings(x["Pos"], x["Next_4line_Pos"]), axis=1)
    df_co_occur_new["Next_5line_Pos_filter"] = df_co_occur_new.apply(
        lambda x: compare_2_strings(x["Pos"], x["Next_5line_Pos"]), axis=1)
                    # df_co_occur_new["Next_line_Pos_filter"] = np.where((df_co_occur_new["Pos"] >= df_co_occur_new["Next_line_Pos"]), 1, 0)
    df_co_occur_new.to_csv(input_dir + "/all_co_occur_nonadar_grouped.csv", sep=",", encoding='utf-8')

    df_co_occur_hap = df_co_occur_new[(df_co_occur_new["Next_1line_Pos_filter"].map(len) >= 2) |
                                      (df_co_occur_new["Next_2line_Pos_filter"].map(len) >= 2) |
                                      (df_co_occur_new["Next_3line_Pos_filter"].map(len) >= 2) |
                                      (df_co_occur_new["Next_4line_Pos_filter"].map(len) >= 2) |
                                      (df_co_occur_new["Next_5line_Pos_filter"].map(len) >= 2)]
    df_co_occur_hap["No_of_lines_to_jump"] = np.where(
        df_co_occur_hap.Next_5line_Pos_filter.map(lambda x: len(x)) >= 2, 5, np.where(
            df_co_occur_hap.Next_4line_Pos_filter.map(lambda x: len(x)) >= 2, 4, np.where(
            df_co_occur_hap.Next_3line_Pos_filter.map(lambda x: len(x)) >= 2, 3, np.where(
            df_co_occur_hap.Next_2line_Pos_filter.map(lambda x: len(x)) >= 2, 2, np.where(df_co_occur_hap.Next_1line_Pos_filter.map(lambda x: len(x)) >= 1, 1, 0)))))
    df_jump = df_co_occur_hap[["Pos", "No_of_lines_to_jump"]]

    # Concate df_co_occur_hap and df_co_occur_new base on Pos
    df_co_occur_merge = df_co_occur_new.merge(df_jump, "left",  on="Pos")
    #  TODO:find the bug for the last line in the group and group the groups

    df_co_occur_merge.to_csv(input_dir + "/all_co_occur_nonadar_merge.csv", sep=",", encoding='utf-8')
    # print(df_co_occur_hap.to_string())

    # Plots
    data = pd.read_csv(input_dir + "/all_co_occur_nonadar_merge_edited.csv", sep=",", encoding='utf-8')
    data["group"] = data["group"].fillna(0)
    data = data.loc[data.group > 0]
    data["group"] = data["group"].astype(int)
    data["group"] = data["group"].astype(str)
    data["group"] = data["group"].apply(lambda x: "Group " + x)
    data["passage"] = data["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    data["passage"] = data["passage"].astype(int)
    label_order = ["RVB14-p2", "RVB14-p5", "RVB14-p8", "RVB14-p10", "RVB14-p12"]

    g4 = sns.lineplot(x="passage", y="Frequency", data=data, hue="group", palette="tab20")# ,kind="point", order=label_order
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(output_dir + "/co_occur_groups.png", dpi=300)

    g5 = sns.regplot(x="passage", y="Frequency", data=data)
    plt.tight_layout()
    g5.set(ylim=(0.0001, 0.0045))
    plt.savefig(output_dir + "/co_occur_reg.png", dpi=300)


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