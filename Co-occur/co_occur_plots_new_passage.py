#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import sys, argparse
import matplotlib
import os
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
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages"
    output_dir = input_dir + "/20201027_co_occur_plots"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    data_p2_1 = pd.read_csv(input_dir + "/p2_1/20201012_q38/co_occur_all.csv")
    data_p2_2 = pd.read_csv(input_dir + "/p2_2/20201012_q38/co_occur_all.csv")
    data_p2_3 = pd.read_csv(input_dir + "/p2_3/20201012_q38/co_occur_all.csv")
    data_p5_1 = pd.read_csv(input_dir + "/p5_1/20201012_q38/co_occur_all.csv")
    data_p5_2 = pd.read_csv(input_dir + "/p5_2/20201012_q38/co_occur_all.csv")
    data_p5_3 = pd.read_csv(input_dir + "/p5_3/20201012_q38/co_occur_all.csv")
    data_p8_1 = pd.read_csv(input_dir + "/p8_1/20201012_q38/co_occur_all.csv")
    data_p8_2 = pd.read_csv(input_dir + "/p8_2/20201012_q38/co_occur_all.csv")
    data_p8_3 = pd.read_csv(input_dir + "/p8_3/20201012_q38/co_occur_all.csv")
    data_p10_1 = pd.read_csv(input_dir + "/p10_1/20201012_q38/co_occur_all.csv")
    data_p10_2 = pd.read_csv(input_dir + "/p10_2/20201012_q38/co_occur_all.csv")
    data_p10_3 = pd.read_csv(input_dir + "/p10_3/20201012_q38/co_occur_all.csv")
    data_p12_1 = pd.read_csv(input_dir + "/p12_1/20201012_q38/co_occur_all.csv")
    data_p12_2 = pd.read_csv(input_dir + "/p12_2/20201012_q38/co_occur_all.csv")
    data_p12_3 = pd.read_csv(input_dir + "/p12_3/20201012_q38/co_occur_all.csv")

    data_p0 = pd.read_csv(input_dir + "/RVB14_p0/20191029_q38/co_occur.csv")
    df = pd.concat([data_p0, data_p2, data_p5, data_p8, data_p10, data_p12])
    df.to_csv(output_dir + "/all_co_occur_raw.csv", sep=",", encoding='utf-8')

    reg_lst = [629, 835, 1621, 2329, 3196, 3634, 3925, 4915, 5170, 5239, 5785, 7165]
    df = add_Protein_to_pd_df.add_Protein_to_pd_df_func(df, reg_lst)

    grouped = df.groupby(["Pos"])
    df_co_occur = pd.DataFrame(grouped.size().reset_index(name="Group_Count"))
    df_co_occur = df_co_occur.merge(df, how="outer", on="Pos")
    df_co_occur = df_co_occur.loc[df_co_occur.Group_Count > 1]


    # df_co_occur["no_variants"] = df_co_occur["Frequency"] * df_co_occur["Read_count"]
    # df_co_occur["frac_and_weight"] = list(zip(df_co_occur.no_variants, df_co_occur.Read_count))
    # df_co_occur = df_co_occur[df_co_occur["Protein"] != "3'UTR"]
    df_co_occur["Passage"] = df_co_occur["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    df_co_occur["Passage"] = np.where(df_co_occur["Passage"] == "RNA Control", 0, df_co_occur["Passage"])
    df_co_occur["Passage"] = df_co_occur["Passage"].astype(int)
    
    df_co_occur.to_csv(output_dir + "/all_co_occur_protein.csv", sep=",", encoding='utf-8')
    sample_order = ["RVB14-p0", "RVB14-p2", "RVB14-p5", "RVB14-p8", "RVB14-p10", "RVB14-p12"]

    #
    # # Plots
    # # g1 = sns.relplot(x="Pos", y="Frequency", data=df_co_occur, hue="label", col="Group_Count", col_wrap=2,
    # #                  hue_order=sample_order)  # , style="Stretch")
    # # g1.set(yscale="log")
    # # # plt.show()
    # # g1.savefig(output_dir + "/co_occur.png", dpi=300)
    # # plt.close()
    #
    # g2 = sns.relplot(x="Pos", y="Frequency", data=df_co_occur, hue="Protein", style="Group_Count", palette="tab10",
    #                  col="passage", col_wrap=3)
    # g2.set(yscale="log")
    # g2.set(xlim=(3500, 7500))
    # # plt.show()
    # g2.savefig(output_dir + "/co_occur_protein.png", dpi=300)
    # plt.close()

    g2_passage = sns.relplot(x="Pos", y="Frequency", data=df_co_occur, hue="Protein", style="Group_Count",
                             palette="tab10", col="Passage", col_wrap=3)#, style="Stretch")
    g2_passage.set(yscale="log")
    g2_passage.set(xlim=(3500, 7500))
    g2_passage.set_xticklabels(fontsize=10, rotation=45)
    # plt.show()
    g2_passage.savefig(output_dir + "/co_occur_protein_passage.png", dpi=300)
    plt.close()

    # g3 = sns.catplot("Pos", "frac_and_weight", data=df_co_occur.loc[df_co_occur.Group_Count == 5], hue="Protein",
    #                  kind="point", dodge=True, join=False, estimator=new_analysis_RV.weighted_varaint, orient="v")
    # g3.set_axis_labels("", "Variant Frequency")
    # # g1.set_xticklabels(fontsize=9, rotation=45)
    # g3.set(yscale='log')
    # # g1.set(ylim=(10**-7, 10**-3))
    # g3.savefig(output_dir + "/co_occur_protein_variant.png", dpi=300)

    # df["Pos"] = df["Pos"].astype(str)
    # df_co_occur_new = (df.groupby(["label", "Stretch"]).agg(lambda x: x.mean() if np.issubdtype(x.dtype, np.number)
    # else ', '.join(x))).reset_index()
    # df_co_occur_new = df_co_occur_new[["label", "Stretch", "Pos", "Frequency"]]
    # df_co_occur_new = df_co_occur_new.sort_values(by=["Pos", "label"])
    # df_co_occur_new["filter"] = df_co_occur_new["Pos"].str.contains(",")
    # df_co_occur_new = df_co_occur_new[df_co_occur_new["filter"] == True]
    # df_co_occur_new["Next_1line_Pos"] = df_co_occur_new["Pos"].shift(periods=-1)
    # df_co_occur_new["Next_2line_Pos"] = df_co_occur_new["Pos"].shift(periods=-2)
    # df_co_occur_new["Next_3line_Pos"] = df_co_occur_new["Pos"].shift(periods=-3)
    # df_co_occur_new["Next_4line_Pos"] = df_co_occur_new["Pos"].shift(periods=-4)
    # df_co_occur_new["Next_5line_Pos"] = df_co_occur_new["Pos"].shift(periods=-5)
    # df_co_occur_new = df_co_occur_new.fillna("0, ")
    #
    # # finds the intersection group
    # df_co_occur_new["Next_1line_Pos_filter"] = df_co_occur_new.apply(
    #     lambda x: compare_2_strings(x["Pos"], x["Next_1line_Pos"]), axis=1)
    # df_co_occur_new["Next_2line_Pos_filter"] = df_co_occur_new.apply(
    #     lambda x: compare_2_strings(x["Pos"], x["Next_2line_Pos"]), axis=1)
    # df_co_occur_new["Next_3line_Pos_filter"] = df_co_occur_new.apply(
    #     lambda x: compare_2_strings(x["Pos"], x["Next_3line_Pos"]), axis=1)
    # df_co_occur_new["Next_4line_Pos_filter"] = df_co_occur_new.apply(
    #     lambda x: compare_2_strings(x["Pos"], x["Next_4line_Pos"]), axis=1)
    # df_co_occur_new["Next_5line_Pos_filter"] = df_co_occur_new.apply(
    #     lambda x: compare_2_strings(x["Pos"], x["Next_5line_Pos"]), axis=1)
    #                 # df_co_occur_new["Next_line_Pos_filter"] = np.where((df_co_occur_new["Pos"] >= df_co_occur_new["Next_line_Pos"]), 1, 0)
    # df_co_occur_new.to_csv(input_dir + "/all_co_occur_grouped.csv", sep=",", encoding='utf-8')
    #
    # df_co_occur_hap = df_co_occur_new[(df_co_occur_new["Next_1line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_2line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_3line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_4line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_5line_Pos_filter"].map(len) >= 2)]
    #                 # df_co_occur_hap["Next_1line_Pos_filter"] = np.where(
    #                 #     df_co_occur_hap.Next_1line_Pos_filter.map(lambda x: len(x)) >= 2, df_co_occur_hap["Next_1line_Pos_filter"], "[]")
    #                 # df_co_occur_hap["Next_2line_Pos_filter"] = np.where(
    #                 #     df_co_occur_hap.Next_2line_Pos_filter.map(lambda x: len(x)) >= 2, df_co_occur_hap["Next_2line_Pos_filter"], "[]")
    #                 # df_co_occur_hap["Next_3line_Pos_filter"] = np.where(
    #                 #     df_co_occur_hap.Next_3line_Pos_filter.map(lambda x: len(x)) >= 2, df_co_occur_hap["Next_3line_Pos_filter"], "[]")
    #                 # df_co_occur_hap["Next_4line_Pos_filter"] = np.where(
    #                 #     df_co_occur_hap.Next_4line_Pos_filter.map(lambda x: len(x)) >= 2, df_co_occur_hap["Next_4line_Pos_filter"], "[]")
    # df_co_occur_hap["No_of_lines_to_jump"] = np.where(
    #     df_co_occur_hap.Next_5line_Pos_filter.map(lambda x: len(x)) >= 2, 5, np.where(
    #         df_co_occur_hap.Next_4line_Pos_filter.map(lambda x: len(x)) >= 2, 4, np.where(
    #         df_co_occur_hap.Next_3line_Pos_filter.map(lambda x: len(x)) >= 2, 3, np.where(
    #         df_co_occur_hap.Next_2line_Pos_filter.map(lambda x: len(x)) >= 2, 2, np.where(df_co_occur_hap.Next_1line_Pos_filter.map(lambda x: len(x)) >= 1, 1, 0)))))
    # df_jump = df_co_occur_hap[["Pos", "No_of_lines_to_jump"]]
    #
    # # Concate df_co_occur_hap and df_co_occur_new base on Pos
    # df_co_occur_merge = df_co_occur_new.merge(df_jump, "left",  on="Pos")
    #
    #                 # df_co_occur_hap["Max"] = df_co_occur_hap.apply(
    #                 #     lambda x: find_max_per_column(x["Next_2line_Pos_filter"], x["Next_1line_Pos_filter"]), axis=1)
    #                 # df_co_occur_hap["Max"] = np.where(
    #                 #     df_co_occur_hap.Max.map(lambda x: len(x)) > df_co_occur_hap.Next_3line_Pos_filter.map(
    #                 #         lambda x: len(x)), df_co_occur_hap["Max"], df_co_occur_hap["Next_3line_Pos_filter"])
    #                 # df_co_occur_hap["Max"] = np.where(
    #                 #     df_co_occur_hap.Max.map(lambda x: len(x)) > df_co_occur_hap.Next_4line_Pos_filter.map(
    #                 #         lambda x: len(x)), df_co_occur_hap["Max"], df_co_occur_hap["Next_4line_Pos_filter"])
    #                 # df_co_occur_hap["Max"] = np.where(
    #                 #     df_co_occur_hap.Max.map(lambda x: len(x)) > df_co_occur_hap.Next_5line_Pos_filter.map(
    #                 #         lambda x: len(x)), df_co_occur_hap["Max"], df_co_occur_hap["Next_5line_Pos_filter"])
    #                 # df_co_occur_hap["Max"] = df_co_occur_hap.apply(
    #                 #     lambda x: find_max_per_column(x["Next_3line_Pos_filter"], x["Next_2line_Pos_filter"]), axis=1)
    #                 # df_co_occur_hap["Max"] = df_co_occur_hap.apply(
    #                 #     lambda x: find_max_per_column(x["Next_4line_Pos_filter"], x["Next_3line_Pos_filter"]), axis=1)
    #                 # df_co_occur_hap["Max"] = df_co_occur_hap.Next_3line_Pos_filter.map(lambda x: len(x))#.max()
    #  TODO:find the bug for the last line in the group and group the groups

    # df_co_occur_merge.to_csv(input_dir + "/all_co_occur_merge.csv", sep=",", encoding='utf-8')
    # print(df_co_occur_hap.to_string())

    data = pd.read_csv(input_dir + "/all_co_occur_merge_edited.csv", sep=",", encoding='utf-8')
    data["group"] = data["group"].fillna(0)
    data = data.loc[data.group > 0]
    data["group"] = data["group"].astype(int)
    data["group"] = data["group"].astype(str)
    data["group"] = data["group"].apply(lambda x: "Group " + x)
    data["Passage"] = data["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    data["Passage"] = data["Passage"].astype(int)
    merge_df = df_co_occur.merge(data, how="outer", on="Stretch")
    columns = ["Pos_x", "Base", "Frequency_x", "Ref", "Read_count", "Rank", "Prob", "Mutation", "Stretch", "meandist",
               "Co-occurrences_identified", "ADAR_context", "ADAR_reverse_context", "Editing_context", "ADAR", "label_x",
               "Protein", "Passage_x", "group"]
    merge_df = merge_df[columns]
    merge_df = merge_df.rename(columns={"Pos_x": "Pos"})
    merge_df = merge_df.rename(columns={"Frequency_x": "Frequency"})
    merge_df = merge_df.rename(columns={"label_x": "label"})
    merge_df = merge_df.rename(columns={"Passage_x": "Passage"})
    merge_df["group"] = merge_df["group"].fillna("out")
    merge_df = merge_df[merge_df["group"] != "out"]
    merge_df = merge_df[merge_df["Protein"] != "3'UTR"]
    merge_df = merge_df[merge_df["Passage"] != 0]
    merge_df = merge_df[merge_df["group"] != "Group 13"]

    group_order = ["Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7", "Group 8", "Group 9",
                   "Group 10", "Group 11", "Group 12"]

    # data = data.rename(columns={"passage": "Passage"})
    g4 = sns.lineplot(x="Passage", y="Frequency", data=data, hue="group", palette="tab10")# ,kind="point", order=label_order

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(output_dir + "/co_occur_groups.png", dpi=300)
    #
    # g5 = sns.regplot(x="Passage", y="Frequency", data=data)
    # plt.tight_layout()
    # g5.set(ylim=(0.0001, 0.0045))
    # plt.savefig(output_dir + "/co_occur_reg.png", dpi=300)

    g6_passage = sns.relplot(x="Pos", y="Frequency", data=merge_df, hue="group", hue_order=group_order,
                             palette="tab20", col="Passage", col_wrap=3, legend="brief")#, style="Stretch")
    g6_passage.set(yscale="log")
    g6_passage.set(ylim=(10 ** -6, 10 ** -1))
    g6_passage.set(xlim=(3500, 7500))
    g6_passage.set_xticklabels(fontsize=10, rotation=45)
    # plt.show()
    g6_passage.savefig(output_dir + "/co_occur_protein_passage_new.png", dpi=300)
    plt.close()

    g7_passage = sns.relplot(x="Passage", y="Frequency", data=merge_df, kind="line", palette="tab20", col="group",
                             col_wrap=4, col_order=group_order, hue="group", legend=False)
    g7_passage.set(yscale="log")
    g7_passage.set(ylim=(10 ** -6, 10 ** -1))
    # g7_passage.set(xlim=(3500, 7500))
    # g7_passage.set_xticklabels(fontsize=10, rotation=45)
    # plt.show()
    g7_passage.savefig(output_dir + "/co_occur_group_new.png", dpi=300)
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