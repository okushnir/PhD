#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
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
    output_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/20201009_plots"
    data_p2 = pd.read_csv(input_dir + "/RVB14_p2/20191029_q38/co_occur.csv")
    data_p5 = pd.read_csv(input_dir + "/RVB14_p5/20191029_q38/co_occur.csv")
    data_p8 = pd.read_csv(input_dir + "/RVB14_p8/20191029_q38/co_occur.csv")
    data_p10 = pd.read_csv(input_dir + "/RVB14_p10/20191029_q38/co_occur.csv")
    data_p12 = pd.read_csv(input_dir + "/RVB14_p12/20191029_q38/co_occur.csv")
    data_p0 = pd.read_csv(input_dir + "/RVB14_p0/20191029_q38/co_occur.csv")
    df = pd.concat([data_p0, data_p2, data_p5, data_p8, data_p10, data_p12])
    df.to_csv(output_dir + "/all_co_occur_raw.csv", sep=",", encoding='utf-8')

    df["Protein"] = np.where(df["Pos"] <= 503, "5'UTR",
                                np.where(df["Pos"] <= 709, "1A",
                                    np.where(df["Pos"] <= 1495, "1B",
                                        np.where(df["Pos"] <= 2197, "1C",
                                            np.where(df["Pos"] <= 3094, "1D",
                                                np.where(df["Pos"] <= 3523, "2A",
                                                    np.where(df["Pos"] <= 3808, "2B",
                                                        np.where(df["Pos"] <= 4771, "2C",
                                                            np.where(df["Pos"] <= 5005, "3A",
                                                                np.where(df["Pos"] <= 5068, "3B",
                                                                   np.where(df["Pos"] <= 5617, "3C",
                                                                    np.where(df["Pos"] <= 6728, "3D", "3'UTR"))))))))))))

    grouped = df.groupby(["Pos"])
    df_co_occur = pd.DataFrame(grouped.size().reset_index(name="Group_Count"))
    df_co_occur = df_co_occur.merge(df, how="outer", on="Pos")
    df_co_occur = df_co_occur.loc[df_co_occur.Group_Count > 1]


    df_co_occur["no_variants"] = df_co_occur["Frequency"] * df_co_occur["Read_count"]
    df_co_occur["frac_and_weight"] = list(zip(df_co_occur.no_variants, df_co_occur.Read_count))
    # df_co_occur = df_co_occur[df_co_occur["Protein"] != "3'UTR"]
    df_co_occur.to_csv(output_dir + "/all_co_occur_protein.csv", sep=",", encoding='utf-8')
    sample_order = ["RVB14-p0", "RVB14-p2", "RVB14-p5", "RVB14-p8", "RVB14-p10", "RVB14-p12"]
    #
    #
    # # Plots
    # # g1 = sns.relplot(x="Pos", y="Frequency", data=df_co_occur, hue="label", col="Group_Count", col_wrap=2,
    # #                  hue_order=sample_order)  # , style="Stretch")
    # # g1.set(yscale="log")
    # # # plt.show()
    # # g1.savefig(output_dir + "/co_occur.png", dpi=300)
    # # plt.close()
    #
    g2 = sns.relplot(x="Pos", y="Frequency", data=df_co_occur, hue="Protein", size="Group_Count", palette="tab10")#, col="Group_Count", col_wrap=2, style="Stretch")
    g2.set(yscale="log")
    g2.set(xlim=(3500, 7500))
    # plt.show()
    g2.savefig(output_dir + "/co_occur_protein.png", dpi=300)
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
    data["passage"] = data["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    data["passage"] = data["passage"].astype(int)
    label_order = ["RVB14-p2", "RVB14-p5", "RVB14-p8", "RVB14-p10", "RVB14-p12"]

    g4 = sns.lineplot(x="passage", y="Frequency", data=data, hue="group", palette="tab10")# ,kind="point", order=label_order

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