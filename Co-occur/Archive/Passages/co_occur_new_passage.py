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
    replica_lst = ["1", "2", "3"]
    for replica in replica_lst:
        data_p2 = pd.read_csv(input_dir + "/p2_%s/20201012_q38/co_occur_all.csv" % replica)
        data_p5 = pd.read_csv(input_dir + "/p5_%s/20201012_q38/co_occur_all.csv" % replica)
        data_p8 = pd.read_csv(input_dir + "/p8_%s/20201012_q38/co_occur_all.csv" % replica)
        data_p10 = pd.read_csv(input_dir + "/p10_%s/20201012_q38/co_occur_all.csv" % replica)
        data_p12 = pd.read_csv(input_dir + "/p12_%s/20201012_q38/co_occur_all.csv" % replica)

        data_p0 = pd.read_csv("/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/"
                                          "merged/controls/IVT_3_Control/20201012_q38/co_occur_all.csv")
        data_rnd = pd.read_csv("/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/190627_RV_CV/merged/RVB14/RVB14_p0/"
                               "20191029_q38/co_occur_all.csv")
        if replica == "3":
            df = pd.concat([data_rnd, data_p2, data_p5, data_p8, data_p10, data_p12])
        else:
            df = pd.concat([data_p0, data_p2, data_p5, data_p8, data_p10, data_p12])
        df = df[["Pos", "Rank", "Stretch", "meandist", "Co-occurrences_identified",	"APOBEC3G_context", "APOBEC3F_context",
                 "ADAR_context", "ADAR_reverse_context", "Editing_context", "ADAR", "label"]]
        data_mutation = pd.read_csv(input_dir + "/inosine_predict_context" + "/data_filter.csv")
        data_mutation[data_mutation["replica"] == int(replica)]
        data_mutation["Pos"] = data_mutation["Pos"].astype(int)
        df["Pos"] = df["Pos"].astype(int)
        df = df.merge(data_mutation, how="right", on=["Pos", "label", "Rank"])
        df = df.sort_values(by=["Pos", "label", "Stretch"])
        df["numric_ref"] = np.where(df["Ref"] == "A", 1, 0)
        df["numric_ref"] = np.where(df["Ref"] == "C", 2, df["numric_ref"])
        df["numric_ref"] = np.where(df["Ref"] == "U", 3, df["numric_ref"])
        df["numric_ref"] = np.where(df["Ref"] == "G", 4, df["numric_ref"])
        df["Stretch"] = df["Stretch"].fillna("None")
        df["Stretch_new"] = np.where(df["Stretch"] != "None", True, False)
        df.to_csv(output_dir + "/all_co_occur_raw_replica%s.csv" % replica, sep=",", encoding='utf-8')
        df.to_pickle(output_dir + "/all_co_occur_raw_replica%s.pkl" % replica)

        # df["Pos"] = df["Pos"].astype(str)
        # df_co_occur_new = (df.groupby(["label", "Stretch"]).agg(lambda x: x.mean() if np.issubdtype(x.dtype, np.number)
        # else ', '.join(x))).reset_index()
        # df_co_occur_new = df_co_occur_new.rename(columns={"Freq": "Frequency", "Base_x": "Base", "Ref_x": "Ref"})
        # df_co_occur_new = df_co_occur_new[["label", "Stretch", "Pos", "Frequency", "Type", "Base", "Ref"]]
        # df_co_occur_new = df_co_occur_new.sort_values(by=["Pos", "label"])
        # df_co_occur_new["filter"] = df_co_occur_new["Pos"].str.contains(",")
        # df_co_occur_new = df_co_occur_new[df_co_occur_new["filter"] == True]
        # df_co_occur_new["Next_line_Pos"] = df_co_occur_new["Pos"].shift(periods=-1)
        # df_co_occur_new = df_co_occur_new.fillna("0, ")
        #
        # grouped2 = df_co_occur_new.groupby(["Pos"])
        # df_grouped = pd.DataFrame(grouped2.size().reset_index(name="Group_Count"))
        # df_grouped["Next_line_Pos"] = df_grouped["Pos"].shift(periods=-1)
        # df_grouped = df_grouped.fillna("0, ")
        # df_grouped["Next_line_Pos_intercept"] = df_grouped.apply(
        #     lambda x: compare_2_strings(x["Pos"], x["Next_line_Pos"]), axis=1)
        # df_grouped = df_grouped.rename_axis('index1').reset_index()
        # df_grouped["index1"] = df_grouped.apply(lambda x: np.where(x.Next_line_Pos_intercept != [], x.index1 + 1, x.index1), axis=1)
        # df_grouped["Group"] = df_grouped["index1"].apply(lambda x: str(x))
        #
        # df_co_occur_merge = df_co_occur_new.merge(df_grouped, on="Pos", how="left")
        # grouped3 = df_co_occur_merge.groupby(["Group"])
        # df_co_occur_group = pd.DataFrame(grouped3.size().reset_index(name="Count"))
        # df_co_occur_merge2 = df_co_occur_merge.merge(df_co_occur_group, how="left", on="Group")
        # df_co_occur_merge2 = df_co_occur_merge2.loc[df_co_occur_merge2.Count > 1]
        # df_co_occur_merge2["Stretch_with_Non-Synonymous"] = df_co_occur_merge2["Type"].str.contains("Non-Synonymous")
        # df_co_occur_merge2["Passage"] = df_co_occur_merge2["label"].apply(lambda x: x.split("-")[0].split("p")[-1])
        # df_co_occur_merge2["Replica"] = df_co_occur_merge2["label"].apply(lambda x: x.split("-")[-1])
        # df_co_occur_merge2["Passage"] = np.where(df_co_occur_merge2["label"] == "RNA Control\nPrimer ID", 0, df_co_occur_merge2["Passage"])
        # df_co_occur_merge2["Replica"] = np.where(df_co_occur_merge2["label"] == "RNA Control\nPrimer ID", 1, df_co_occur_merge2["Replica"])
        # df_co_occur_merge2["Passage"] = np.where(df_co_occur_merge2["label"] == "RNA Control_RND", 0, df_co_occur_merge2["Passage"])
        # df_co_occur_merge2["Replica"] = np.where(df_co_occur_merge2["label"] == "RNA Control_RND", 3, df_co_occur_merge2["Replica"])
        # df_co_occur_merge2["Passage"] = df_co_occur_merge2["Passage"].astype(int)
        # df_co_occur_merge2["Replica"] = df_co_occur_merge2["Replica"].astype(int)
        # print(df_co_occur_merge2.to_string())
        # df_co_occur_merge2.to_csv(output_dir + "/all_co_occur_merge_replica%s.csv" % replica, sep=",", encoding='utf-8')


        # # Plots


        # sample_order = ["RVB14-p0", "RVB14-p2", "RVB14-p5", "RVB14-p8", "RVB14-p10", "RVB14-p12"]
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

        # g2_passage = sns.relplot(x="Pos", y="Frequency", data=df_co_occur, hue="Protein",
        #                          palette="tab10", col="Passage", col_wrap=3)# , style="Group_Count"
        # g2_passage.set(yscale="log")
        # g2_passage.set(xlim=(3500, 7500))
        # g2_passage.set_xticklabels(fontsize=10, rotation=45)
        # # plt.show()
        # g2_passage.savefig(output_dir + "/co_occur_protein_passage.png", dpi=300)
        # plt.close()

        # g3 = sns.catplot("Pos", "frac_and_weight", data=df_co_occur.loc[df_co_occur.Group_Count == 5], hue="Protein",
        #                  kind="point", dodge=True, join=False, estimator=new_analysis_RV.weighted_varaint, orient="v")
        # g3.set_axis_labels("", "Variant Frequency")
        # # g1.set_xticklabels(fontsize=9, rotation=45)
        # g3.set(yscale='log')
        # # g1.set(ylim=(10**-7, 10**-3))
        # g3.savefig(output_dir + "/co_occur_protein_variant.png", dpi=300)


        # data = pd.read_csv(input_dir + "/all_co_occur_merge_edited.csv", sep=",", encoding='utf-8')
        # data["group"] = data["group"].fillna(0)
        # data = data.loc[data.group > 0]
        # data["group"] = data["group"].astype(int)
        # data["group"] = data["group"].astype(str)
        # data["group"] = data["group"].apply(lambda x: "Group " + x)
        # data = df_co_occur_merge2

        # merge_df = df_co_occur.merge(data, how="outer", on="Stretch")
        # columns = ["Pos_x", "Base", "Frequency_x", "Ref", "Read_count", "Rank", "Prob", "Mutation", "Stretch", "meandist",
        #            "Co-occurrences_identified", "ADAR_context", "ADAR_reverse_context", "Editing_context", "ADAR", "label_x",
        #            "Protein", "Passage_x", "group"]
        # merge_df = merge_df[columns]
        # merge_df = merge_df.rename(columns={"Pos_x": "Pos"})
        # merge_df = merge_df.rename(columns={"Frequency_x": "Frequency"})
        # merge_df = merge_df.rename(columns={"label_x": "label"})
        # merge_df = merge_df.rename(columns={"Passage_x": "Passage"})
        # merge_df["group"] = merge_df["group"].fillna("out")
        # merge_df = merge_df[merge_df["group"] != "out"]
        # merge_df = merge_df[merge_df["Protein"] != "3'UTR"]
        # merge_df = merge_df[merge_df["Passage"] != 0]
        # merge_df = merge_df[merge_df["group"] != "Group 13"]

        # group_order = ["Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7", "Group 8", "Group 9",
        #                "Group 10", "Group 11", "Group 12"]

        # data = data.rename(columns={"passage": "Passage"})
        # g4 = sns.lineplot(x="Passage", y="Frequency", data=data, hue="Group", palette="tab10")# ,kind="point", order=label_order
        #
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        # plt.tight_layout()
        # plt.savefig(output_dir + "/co_occur_groups_replica%s.png" % replica, dpi=300)
        # plt.close()
        # # g5 = sns.regplot(x="Passage", y="Frequency", data=data)
        # # plt.tight_layout()
        # # g5.set(ylim=(0.0001, 0.0045))
        # # plt.savefig(output_dir + "/co_occur_reg.png", dpi=300)
        #
        # g6_passage = sns.relplot(x="Passage", y="Frequency", data=data, hue="Stretch_with_Non-Synonymous",
        #                          palette="Set2", col="Group", col_wrap=4, legend="brief", height=5)#, style="Stretch")
        # g6_passage.set(yscale="log")
        # g6_passage.set(ylim=(10 ** -6, 10 ** -1))
        # # g6_passage.set(xlim=(3500, 7500))
        # # g6_passage.set_xticklabels(fontsize=10, rotation=45)
        # # plt.show()
        # g6_passage.savefig(output_dir + "/co_occur_replica%s.png" % replica, dpi=600 )
        # plt.close()
        #
        # g7_passage = sns.relplot(x="Passage", y="Frequency", data=merge_df, kind="line", palette="tab20", col="group",
        #                          col_wrap=4, col_order=group_order, hue="group", legend=False)
        # g7_passage.set(yscale="log")
        # g7_passage.set(ylim=(10 ** -6, 10 ** -1))
        # # g7_passage.set(xlim=(3500, 7500))
        # # g7_passage.set_xticklabels(fontsize=10, rotation=45)
        # # plt.show()
        # g7_passage.savefig(output_dir + "/co_occur_group_new.png", dpi=300)
        # plt.close()

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