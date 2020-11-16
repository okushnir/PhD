#! /usr/local/python/anaconda_python-3.6.1

"""
@Author: odedkushnir

"""


import os.path
import glob
import pandas as pd
import numpy as np
from AccuNGS_analysis.add_Protein_to_pd_df import add_Protein_to_pd_df_func
import seaborn as sns
from AccuNGS_analysis import old_statannot
import matplotlib.pyplot as plt

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

def creating_co_occur_passages_df(input_dir, prefix, date, q, region_lst):

    # input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14"

    output_dir = input_dir + "/Co_occur_all_passages"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    dirs = glob.glob(input_dir + prefix)
    lst_srr = []
    for passage in dirs:
        file_path = glob.glob(passage + "/%s_%s/co_occur_all.csv" % (date, q))
        lst_srr.append(file_path[0])
    print(lst_srr)

    columns = ["Pos", "Base", "Frequency", "Ref", "Read_count",	"Rank",	"Prob",	"Mutation",	"Stretch", "meandist",
               "Co-occurrences_identified",	"ADAR_context",	"ADAR_reverse_context",	"APOBEC3G_context",
               "APOBEC3F_context", "Editing_context", "ADAR", "label"]
    data = pd.DataFrame(columns=columns)
    for i in range(len(lst_srr)):
        sample_file = lst_srr[i]
        print("loading " + sample_file + " as sample")
        data_co_occur_all = pd.read_csv(sample_file, sep=",")
        data = data.append(data_co_occur_all, ignore_index=True)

    data = add_Protein_to_pd_df_func(data, region_lst)
    grouped = data.groupby(["Pos"])
    df_co_occur = pd.DataFrame(grouped.size().reset_index(name="Group_Count"))
    df_co_occur = df_co_occur.merge(data, how="outer", on="Pos")
    df_co_occur = df_co_occur.loc[df_co_occur.Group_Count > 1]
    df_co_occur["Passage"] = df_co_occur["label"].apply(lambda x: x.split("-")[0][1:])
    # df_co_occur["Passage"] = np.where(df_co_occur["Passage"] == "RNA Control", 0, df_co_occur["Passage"])
    df_co_occur["replica"] = df_co_occur["label"].apply(lambda x: x.split("-")[-1])


    df_co_occur["Passage"] = df_co_occur["Passage"].astype(str)
    df_co_occur["replica"] = df_co_occur["replica"].astype(str)
    df_co_occur["New_Stretch"] = df_co_occur["Stretch"] + "_" + df_co_occur["Passage"] + "_" + df_co_occur["replica"]
    df_co_occur["Passage"] = df_co_occur["Passage"].astype(int)
    df_co_occur["replica"] = df_co_occur["replica"].astype(int)

    df_co_occur.to_csv(output_dir + "/all_co_occur_protein.csv", sep=",", encoding='utf-8')
    df_co_occur.to_pickle(output_dir + "/all_co_occur_protein.pkl")
    return df_co_occur

def grouped_co_occur(df, input_dir,output_dir):
    data_mutation = pd.read_pickle(input_dir +"/Rank0_data_mutation" + "/q38_data_mutation.pkl")
    data_mutation["Pos"] = data_mutation["Pos"].astype(int)
    df = df.merge(data_mutation, "inner", on=["Pos", "label", "Rank"])
    df["Pos"] = df["Pos"].astype(str)
    df["Frequency_x"] = df["Frequency_x"].astype(str)
    df["Frequency_y"] = df["Frequency_y"].astype(str)
    df_co_occur_new = (df.groupby(["label", "Stretch"]).agg(lambda x: x.mean() if np.issubdtype(x.dtype, np.number)
    else ', '.join(x))).reset_index()
    df_co_occur_new = df_co_occur_new[["label", "Stretch", "Pos", "Frequency_x", "meandist", "Type", "Mutation_x"]]
    df_co_occur_new = df_co_occur_new.rename(columns={"Frequency_x": "Frequency", "Mutation_x": "Mutation"})
    df_co_occur_new = df_co_occur_new.sort_values(by=["Pos", "label"])
    df_co_occur_new["filter"] = df_co_occur_new["Pos"].str.contains(",")
    df_co_occur_new = df_co_occur_new[df_co_occur_new["filter"] == True]
    df_co_occur_new["Stretch_Freq"] = df_co_occur_new["Frequency"].apply(lambda x: np.mean(list(map(float, x.split(",")))))
    df_co_occur_new["Stretch_Type_Non-Synonymous"] = df_co_occur_new["Type"].str.contains("Non-Synonymous")

    df_co_occur_new.to_csv(output_dir + "/all_co_occur_grouped.csv", sep=",", encoding='utf-8')
    df_co_occur_new.to_pickle(output_dir + "/all_co_occur_grouped.pkl")

    g1 = sns.boxplot(x="Stretch_Type_Non-Synonymous", y="Stretch_Freq", data=df_co_occur_new,
                     order=[False, True])
    old_statannot.add_stat_annotation(g1, data=df_co_occur_new, x="Stretch_Type_Non-Synonymous", y="Stretch_Freq",
                        boxPairList=[(False, True)], test='Mann-Whitney', textFormat='star',
                                      loc='outside', verbose=2,
                                      order=[False, True])
    g1.set(xticklabels=("Synonymous\nStretches", "Stretches\nWith Non-Synonymous"))
    g1.set(xlabel="")
    plt.tight_layout()
    plt.savefig(output_dir + "/Frequencies_of_Stertch_Type_Non-Synonymous", dpi=300)

    return df_co_occur_new


def regression_stretches(df_co_occur_new, output_dir, input_dir=None):
    df_pos = pd.read_pickle(output_dir + "/all_co_occur_protein.pkl")
    df_adar = pd.read_csv("/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/"
                             "inosine_predict_context/data_filter.csv", sep=",")
    df_co_occur_new = df_co_occur_new.merge(df_pos, "outer", on=["label", "Stretch"])
    columns = ["label", "Stretch", "meandist_x", "Type", "Mutation_x", "Stretch_Freq", "Stretch_Type_Non-Synonymous",
               "Pos_y", "Base", "Frequency_y", "Ref", "Read_count", "Rank", "Prob", "Mutation_y",
               "Co-occurrences_identified",	"ADAR_context",	"ADAR_reverse_context",	"APOBEC3G_context",
               "APOBEC3F_context", "Editing_context", "ADAR", "Protein", "Passage", "replica"]
    df_co_occur_new = df_co_occur_new[columns]
    df_co_occur_new.rename(columns={"Pos_y": "Pos", "meandist_x": "meandist", "Mutation_x": "Mutation_stretch",
                                    "Frequency_y": "Frequency", "Mutation_y": "Mutation", "Type": "Type_Stretch",
                                    "Passage": "passage"}, inplace=True)
    df_co_occur_new = df_co_occur_new.merge(df_adar, "inner", on=["Pos", "label", "passage", "replica"])
    print(df_co_occur_new.to_string())


def main():
    # for Cluster
    # parser = OptionParser("usage: %prog [options]")
    # parser.add_option("-s", "--SRR_dir", dest="SRR_dir", help="the directory of all the SRR's")
    # parser.add_option("-n", "--SRR_no", dest="SRR_no", help="SRR Number")
    # (options, args) = parser.parse_args()
    #
    # SRR_dir = options.SRR_dir
    # SRR_No = options.SRR_no
    # freqs_file = check_dirname(freqs_file)


    # for Local
    """RV"""
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages"
    output_dir = input_dir + "/Co_occur_all_passages"
    prefix = "/p*"
    min_coverage = 5000
    virus = "RVB14"
    date = "20201012"
    q = "q38"
    region_lst = [629, 835, 1621, 2329, 3196, 3634, 3925, 4915, 5170, 5239, 5785, 7165]
    # data_all_passages = creating_co_occur_passages_df(input_dir, prefix, date, q, region_lst)

    # data_all_passages = pd.read_pickle(output_dir + "/all_co_occur_protein.pkl")
    # df_co_occur_new = grouped_co_occur(data_all_passages, input_dir, output_dir)

    df_co_occur_new = pd.read_pickle(output_dir + "/all_co_occur_grouped.pkl")
    df_regression = regression_stretches(df_co_occur_new, output_dir)


    # """RV-Capsid_Free"""
    # input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/capsid"
    # prefix = "/*_3*"
    # min_coverage = 5000
    # virus = "RVB14"
    # date = "20201012"
    # q = "q38"
    #
    # control_file_id = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/controls/" \
    #                   "IVT_3_Control/20201012_q38/IVT-3-Control.merged.with.mutation.type.freqs"
    # label_control1 = "RNA Control\nPrimer ID"
    # control_file_mix = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/p8_2/" \
    #                    "20201012_q38/p8-2.merged.with.mutation.type.freqs"
    # label_control2 = "p8 Mixed Population"
    # control_dict = {label_control1: control_file_id, label_control2: control_file_mix}
    # creating_data_mutation_df(input_dir, prefix, min_coverage, virus, date, q, control_dict)
    #
    # """RV-Patients"""
    # input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/patients"
    # prefix = "/*"
    # min_coverage = 5000
    # virus = "RVB14"
    # date = "20201017"
    # q = "q30_consensusX5"
    #
    # control_file_id = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/controls/IVT_5_Control/20201012_q38/IVT-5-Control.merged.with.mutation.type.freqs"
    # label_control1 = "RNA Control\nPrimer ID"
    # control_file_cell = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/controls/p3_Control/20201012_q38/p3-Control.merged.with.mutation.type.freqs"
    # label_control2 = "p3 Cell Culture\nControl"
    # control_dict = {label_control1: control_file_id, label_control2: control_file_cell}
    # creating_data_mutation_df(input_dir, prefix, min_coverage, virus, date, q, control_dict)
    #
    # """CV"""
    # input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/CVB3"
    # prefix = "/CVB3_p*"
    # min_coverage = 5000
    # virus = "CVB3"
    # date = "q38"
    # q ="3UTR"
    #
    # control_file = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/CVB3/CVB3_RNA_Control/q38_3UTR/" \
    #                "CVB3-RNA-Control.merged.with.mutation.type.freqs"
    # label_control = "CVB3-RNA Control"
    # control_dict = {label_control: control_file}
    # creating_data_mutation_df(input_dir, prefix, min_coverage, virus, date, q, control_dict)
    #
    # """PV1"""
    # input_dir = "/Users/odedkushnir/Projects/fitness/CirSeq/PV/Mahoney"
    # prefix = "/p*"
    # min_coverage = 10000
    # virus = "Human poliovirus 1"
    # date = "20181210"
    # q = "q30"
    # creating_data_mutation_df(input_dir, prefix, min_coverage, virus, date, q)
    #
    # """OPV2"""
    # input_dir = "/Users/odedkushnir/Projects/fitness/CirSeq/PV/OPV"
    # prefix = "/p*"
    # min_coverage = 10000
    # virus = "OPV"
    # date = "20190226"
    # q = "q23"
    # creating_data_mutation_df(input_dir, prefix, min_coverage, virus, date, q)


if __name__ == "__main__":
    main()

    # df_co_occur_new["Next_1line_Pos"] = df_co_occur_new["Pos"].shift(periods=-1)
    # df_co_occur_new["Next_2line_Pos"] = df_co_occur_new["Pos"].shift(periods=-2)
    # df_co_occur_new["Next_3line_Pos"] = df_co_occur_new["Pos"].shift(periods=-3)
    # df_co_occur_new["Next_4line_Pos"] = df_co_occur_new["Pos"].shift(periods=-4)
    # df_co_occur_new["Next_5line_Pos"] = df_co_occur_new["Pos"].shift(periods=-5)
    # df_co_occur_new["Next_6line_Pos"] = df_co_occur_new["Pos"].shift(periods=-6)
    # df_co_occur_new["Next_7line_Pos"] = df_co_occur_new["Pos"].shift(periods=-7)
    # df_co_occur_new["Next_8line_Pos"] = df_co_occur_new["Pos"].shift(periods=-8)
    # df_co_occur_new["Next_9line_Pos"] = df_co_occur_new["Pos"].shift(periods=-9)
    # df_co_occur_new["Next_10line_Pos"] = df_co_occur_new["Pos"].shift(periods=-10)
    # df_co_occur_new["Next_11line_Pos"] = df_co_occur_new["Pos"].shift(periods=-11)
    # df_co_occur_new["Next_12line_Pos"] = df_co_occur_new["Pos"].shift(periods=-12)
    # df_co_occur_new["Next_13line_Pos"] = df_co_occur_new["Pos"].shift(periods=-13)
    # df_co_occur_new["Next_14line_Pos"] = df_co_occur_new["Pos"].shift(periods=-14)
    # df_co_occur_new["Next_15line_Pos"] = df_co_occur_new["Pos"].shift(periods=-15)
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
    #
    # df_co_occur_new["Next_6line_Pos_filter"] = df_co_occur_new.apply(
    #     lambda x: compare_2_strings(x["Pos"], x["Next_6line_Pos"]), axis=1)
    # df_co_occur_new["Next_7line_Pos_filter"] = df_co_occur_new.apply(
    #     lambda x: compare_2_strings(x["Pos"], x["Next_7line_Pos"]), axis=1)
    # df_co_occur_new["Next_8line_Pos_filter"] = df_co_occur_new.apply(
    #     lambda x: compare_2_strings(x["Pos"], x["Next_8line_Pos"]), axis=1)
    # df_co_occur_new["Next_9line_Pos_filter"] = df_co_occur_new.apply(
    #     lambda x: compare_2_strings(x["Pos"], x["Next_9line_Pos"]), axis=1)
    # df_co_occur_new["Next_10line_Pos_filter"] = df_co_occur_new.apply(
    #     lambda x: compare_2_strings(x["Pos"], x["Next_10line_Pos"]), axis=1)
    #
    # df_co_occur_new["Next_11line_Pos_filter"] = df_co_occur_new.apply(
    #     lambda x: compare_2_strings(x["Pos"], x["Next_11line_Pos"]), axis=1)
    # df_co_occur_new["Next_12line_Pos_filter"] = df_co_occur_new.apply(
    #     lambda x: compare_2_strings(x["Pos"], x["Next_12line_Pos"]), axis=1)
    # df_co_occur_new["Next_13line_Pos_filter"] = df_co_occur_new.apply(
    #     lambda x: compare_2_strings(x["Pos"], x["Next_13line_Pos"]), axis=1)
    # df_co_occur_new["Next_14line_Pos_filter"] = df_co_occur_new.apply(
    #     lambda x: compare_2_strings(x["Pos"], x["Next_14line_Pos"]), axis=1)
    # df_co_occur_new["Next_15line_Pos_filter"] = df_co_occur_new.apply(
    #     lambda x: compare_2_strings(x["Pos"], x["Next_15line_Pos"]), axis=1)
    #
    #                 # df_co_occur_new["Next_line_Pos_filter"] = np.where((df_co_occur_new["Pos"] >= df_co_occur_new["Next_line_Pos"]), 1, 0)
    # df_co_occur_new.to_csv(output_dir + "/all_co_occur_grouped_with_lines.csv", sep=",", encoding='utf-8')
    #
    # df_co_occur_hap = df_co_occur_new[(df_co_occur_new["Next_1line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_2line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_3line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_4line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_5line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_6line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_7line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_8line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_9line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_10line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_11line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_12line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_13line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_14line_Pos_filter"].map(len) >= 2) |
    #                                   (df_co_occur_new["Next_15line_Pos_filter"].map(len) >= 2)]
    #                 # df_co_occur_hap["Next_1line_Pos_filter"] = np.where(
    #                 #     df_co_occur_hap.Next_1line_Pos_filter.map(lambda x: len(x)) >= 2, df_co_occur_hap["Next_1line_Pos_filter"], "[]")
    #                 # df_co_occur_hap["Next_2line_Pos_filter"] = np.where(
    #                 #     df_co_occur_hap.Next_2line_Pos_filter.map(lambda x: len(x)) >= 2, df_co_occur_hap["Next_2line_Pos_filter"], "[]")
    #                 # df_co_occur_hap["Next_3line_Pos_filter"] = np.where(
    #                 #     df_co_occur_hap.Next_3line_Pos_filter.map(lambda x: len(x)) >= 2, df_co_occur_hap["Next_3line_Pos_filter"], "[]")
    #                 # df_co_occur_hap["Next_4line_Pos_filter"] = np.where(
    #                 #     df_co_occur_hap.Next_4line_Pos_filter.map(lambda x: len(x)) >= 2, df_co_occur_hap["Next_4line_Pos_filter"], "[]")
    # df_co_occur_hap["No_of_lines_to_jump"] = np.where(
    #     df_co_occur_hap.Next_15line_Pos_filter.map(lambda x: len(x)) >= 2, 15, np.where(
    #         df_co_occur_hap.Next_14line_Pos_filter.map(lambda x: len(x)) >= 2, 14, np.where(
    #             df_co_occur_hap.Next_13line_Pos_filter.map(lambda x: len(x)) >= 2, 13, np.where(
    #                 df_co_occur_hap.Next_12line_Pos_filter.map(lambda x: len(x)) >= 2, 12, np.where(
    #                     df_co_occur_hap.Next_11line_Pos_filter.map(lambda x: len(x)) >= 2, 11, np.where(
    #                         df_co_occur_hap.Next_10line_Pos_filter.map(lambda x: len(x)) >= 2, 10, np.where(
    #                             df_co_occur_hap.Next_9line_Pos_filter.map(lambda x: len(x)) >= 2, 9, np.where(
    #                                 df_co_occur_hap.Next_8line_Pos_filter.map(lambda x: len(x)) >= 2, 8, np.where(
    #                                     df_co_occur_hap.Next_7line_Pos_filter.map(lambda x: len(x)) >= 2, 7, np.where(
    #                                         df_co_occur_hap.Next_6line_Pos_filter.map(lambda x: len(x)) >= 2, 6, np.where(
    #                                             df_co_occur_hap.Next_5line_Pos_filter.map(lambda x: len(x)) >= 2, 5, np.where(
    #                                              df_co_occur_hap.Next_4line_Pos_filter.map(lambda x: len(x)) >= 2, 4, np.where(
    #                                                 df_co_occur_hap.Next_3line_Pos_filter.map(lambda x: len(x)) >= 2, 3, np.where(
    #                                                      df_co_occur_hap.Next_2line_Pos_filter.map(lambda x: len(x)) >= 2, 2, np.where(
    #                                                         df_co_occur_hap.Next_1line_Pos_filter.map(lambda x: len(x)) >= 1, 1, 0)))))))))))))))
    # df_jump = df_co_occur_hap[["Pos", "No_of_lines_to_jump", "label"]]
    #
    # # Concate df_co_occur_hap and df_co_occur_new base on Pos
    #
    #
    # df_co_occur_hap["Max"] = df_co_occur_hap.apply(
    #     lambda x: find_max_per_column(x["Next_2line_Pos_filter"], x["Next_1line_Pos_filter"]), axis=1)
    # # df_co_occur_hap["Max"] = np.where(
    # #     df_co_occur_hap.Max.map(lambda x: len(x)) > df_co_occur_hap.Next_3line_Pos_filter.map(
    # #         lambda x: len(x)), df_co_occur_hap["Max"], df_co_occur_hap["Next_3line_Pos_filter"])
    # # df_co_occur_hap["Max"] = np.where(
    # #     df_co_occur_hap.Max.map(lambda x: len(x)) > df_co_occur_hap.Next_4line_Pos_filter.map(
    # #         lambda x: len(x)), df_co_occur_hap["Max"], df_co_occur_hap["Next_4line_Pos_filter"])
    # # df_co_occur_hap["Max"] = np.where(
    # #     df_co_occur_hap.Max.map(lambda x: len(x)) > df_co_occur_hap.Next_5line_Pos_filter.map(
    # #         lambda x: len(x)), df_co_occur_hap["Max"], df_co_occur_hap["Next_5line_Pos_filter"])
    # # df_co_occur_hap["Max"] = df_co_occur_hap.apply(
    # #     lambda x: find_max_per_column(x["Next_3line_Pos_filter"], x["Next_2line_Pos_filter"]), axis=1)
    # # df_co_occur_hap["Max"] = df_co_occur_hap.apply(
    # #     lambda x: find_max_per_column(x["Next_4line_Pos_filter"], x["Next_3line_Pos_filter"]), axis=1)
    # # df_co_occur_hap["Max"] = df_co_occur_hap.Next_3line_Pos_filter.map(lambda x: len(x))#.max()
    #  # TODO:find the bug for the last line in the group and group the groups
    #
    # df_co_occur_merge = df_co_occur_new.merge(df_jump, "left", on=["Pos", "label"])
    #
    # df_co_occur_merge.to_csv(output_dir + "/all_co_occur_merge.csv", sep=",", encoding='utf-8')
    # print(df_co_occur_hap.to_string())