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
from statannot import statannot
import matplotlib.pyplot as plt
from scipy import stats

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

def weighted_varaint(x, **kws):
    var, count = map(np.asarray, zip(*x))
    return var.sum() / count.sum()

def creating_co_occur_passages_df(input_dir, experiment, prefix, date, q, region_lst):
    """
    concatenate all co_occur files from co_occur pipeline
    :param input_dir: input directory path, where all the sub directory are located
    :param prefix: the prefix of the sub directory, for example "/p*" in the passage directory
    :param date: the prefix of the directory contains the freqs
    :param q: which q-score was used in the pipeline
    :param region_lst: the protein list
    :return: pd.DataFrame with the stretches
    """

    output_dir = input_dir + "/Co_occur_all_%s" % experiment
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
    df_co_occur = df_co_occur.loc[df_co_occur.Group_Count >= 1]
    if experiment == "capsid":
        df_co_occur["Passage"] = 9
    else:
        df_co_occur["Passage"] = df_co_occur["label"].apply(lambda x: x.split("-")[0][1:])
    # df_co_occur["Passage"] = np.where(df_co_occur["Passage"] == "RNA Control", 0, df_co_occur["Passage"])
    if experiment == "capsid":
        df_co_occur["replica"] = df_co_occur["label"].apply(lambda x: x.split("-")[1][1])
        df_co_occur["replica"] = np.where(df_co_occur["replica"] == "3", "2", df_co_occur["replica"])
    else:
        df_co_occur["replica"] = df_co_occur["label"].apply(lambda x: x.split("-")[-1])


    df_co_occur["Passage"] = df_co_occur["Passage"].astype(str)
    df_co_occur["replica"] = df_co_occur["replica"].astype(str)
    df_co_occur["New_Stretch"] = df_co_occur["Stretch"] + "_" + df_co_occur["Passage"] + "_" + df_co_occur["replica"]
    df_co_occur["Passage"] = df_co_occur["Passage"].astype(int)
    df_co_occur["replica"] = df_co_occur["replica"].astype(int)

    df_co_occur.to_csv(output_dir + "/all_co_occur_protein.csv", sep=",", encoding='utf-8')
    df_co_occur.to_pickle(output_dir + "/all_co_occur_protein.pkl")
    return df_co_occur


def grouped_co_occur(df, input_dir, experiment, output_dir, q):
    """
    :param df: pd.DataFrame with the stretches
    :param input_dir: input directory path, where all the sub directory are located
    :param output_dir: output directory path, all files will be stored there
    :return: pd.DataFrame with aggregated stretches
    """
    data_mutation = pd.read_pickle(input_dir +"/Rank0_data_mutation" + "/%s_data_mutation.pkl" % q)
    data_mutation["Pos"] = data_mutation["Pos"].astype(int)
    df = df.merge(data_mutation, "inner", on=["Pos", "label", "Rank"])
    df["Pos"] = df["Pos"].astype(str)
    df["Frequency_x"] = df["Frequency_x"].astype(str) #Maoz Script
    df["Frequency_y"] = df["Frequency_y"].astype(str) #freqs file
    df["Freq"] = df["Freq"].astype(str)
    df_co_occur_new = (df.groupby(["label", "Stretch"]).agg(lambda x: x.mean() if np.issubdtype(x.dtype, np.number)
    else ', '.join(x))).reset_index()
    df_co_occur_new = df_co_occur_new[["label", "Stretch", "Pos", "Freq", "meandist", "Type", "Mutation_x"]]
    df_co_occur_new = df_co_occur_new.rename(columns={"Freq": "Frequency", "Mutation_x": "Mutation"})
    if experiment == "capsid":
        # df_co_occur_new = df_co_occur_new.loc[df_co_occur_new.label != "Capsid-32-Ultra"]
        # df_co_occur_new = df_co_occur_new.loc[df_co_occur_new.label != "Free-32-Ultra"]
        df_co_occur_new = df_co_occur_new.loc[df_co_occur_new.label != "Free-33-Amicon"]
    df_co_occur_new = df_co_occur_new.sort_values(by=["Pos", "label"])
    df_co_occur_new["filter"] = df_co_occur_new["Pos"].str.contains(",")
    df_co_occur_new = df_co_occur_new[df_co_occur_new["filter"] == True]
    df_co_occur_new["Stretch_Freq"] = df_co_occur_new["Frequency"].apply(lambda x: np.mean(list(map(float, x.split(",")))))
    df_co_occur_new["Stretch_Type_Non-Synonymous"] = df_co_occur_new["Type"].str.contains("Non-Synonymous")
    if experiment == "capsid":
        df_co_occur_new["RNA"] = df_co_occur_new["label"].apply(lambda x: x.split("-")[0])
        df_co_occur_new["replica"] = df_co_occur_new["label"].apply(lambda x: x.split("-")[1][1])
        # df_co_occur_new["replica"] = np.where(df_co_occur_new["replica"] == "3", "2", df_co_occur_new["replica"])
    df_co_occur_new.to_csv(output_dir + "/all_co_occur_grouped.csv", sep=",", encoding='utf-8')
    df_co_occur_new.to_pickle(output_dir + "/all_co_occur_grouped.pkl")
    return df_co_occur_new


def regression_stretches(df_co_occur_new, output_dir, df_adar_path, experiment, min_coverage):
    """
    :param df_co_occur_new: pd.DataFrame with aggregated stretches
    :param output_dir: output directory path, all files will be stored there
    :param df_adar_path: pd.DataFrame contained the adar preference grades
    :param input_dir: input directory path
    :return: mereged pd.DataFrame with stretches and ADAR preferences and for each position
    """
    df_pos = pd.read_pickle(output_dir + "/all_co_occur_protein.pkl")
    df_adar = pd.read_csv(df_adar_path, sep=",")

    df_co_occur_new = df_co_occur_new.merge(df_pos, "outer", on=["label", "Stretch"])
    if experiment == "capsid":
        df_co_occur_new["replica"] = df_co_occur_new["label"].apply(lambda x: x.split("-")[1][1])
        # df_co_occur_new["replica"] = np.where(df_co_occur_new["replica"] == "3", "2", df_co_occur_new["replica"])
        df_co_occur_new["Passage"] = 9
    columns = ["label", "Stretch", "meandist_x", "Type", "Mutation_x", "Stretch_Freq", "Stretch_Type_Non-Synonymous",
               "Pos_y", "Mutation_y", "Co-occurrences_identified", "ADAR_context", "ADAR_reverse_context", "APOBEC3G_context",
               "APOBEC3F_context", "Editing_context", "ADAR", "Passage", "replica"]
    df_co_occur_new = df_co_occur_new[columns]
    df_co_occur_new.rename(columns={"Pos_y": "Pos", "meandist_x": "meandist", "Mutation_x": "Mutation_stretch",
                                    "Type": "Type_Stretch", "Passage": "passage", "Mutation_y": "Mutation"}, inplace=True)
    if experiment == "capsid":
        df_co_occur_new["Pos"] = df_co_occur_new["Pos"].astype(int)
        df_adar["Pos"] = df_adar["Pos"].astype(int)
        df_co_occur_new["label"] = df_co_occur_new["label"].astype(str)
        df_adar["label"] = df_adar["label"].astype(str)
        df_co_occur_new["Mutation"] = df_co_occur_new["Mutation"].astype(str)
        df_adar["Mutation"] = df_adar["Mutation"].astype(str)
        df_co_occur_new["replica"] = df_co_occur_new["replica"].astype(int)
        df_adar["replica"] = df_adar["replica"].astype(int)
        df_co_occur_new = df_co_occur_new.merge(df_adar, "outer", on=["Pos", "label", "Mutation", "replica"])
        df_co_occur_new = df_co_occur_new[df_co_occur_new["label"] != "p8 Mixed Population"]
    else:
        df_co_occur_new = df_co_occur_new.merge(df_adar, "outer", on=["Pos", "label", "Mutation", "passage", "replica"])
    df_co_occur_new["Mutation"] = df_co_occur_new["Mutation"].astype(str)
    df_co_occur_new["fiveGrade"] = df_co_occur_new["fiveGrade"].astype(float)
    df_co_occur_new["5`_ADAR_Preference"] = np.where(df_co_occur_new["fiveGrade"] == 0, "None",
                                                     df_co_occur_new["5`_ADAR_Preference"])
    df_co_occur_new["Stretch_Type_Non-Synonymous"] = np.where((df_co_occur_new["Stretch_Type_Non-Synonymous"] != True) &
                                                              (df_co_occur_new["Stretch_Type_Non-Synonymous"] != False) &
                                                              (df_co_occur_new["Type"] == "Synonymous"),
                                                              "Non-Stretch_Synonymous",
                                                              df_co_occur_new["Stretch_Type_Non-Synonymous"])
    df_co_occur_new["Stretch_Type_Non-Synonymous"] = np.where((df_co_occur_new["Stretch_Type_Non-Synonymous"] != True) &
                                                              (df_co_occur_new["Stretch_Type_Non-Synonymous"] != False) &
                                                              (df_co_occur_new["Type"] == "Synonymous") &
                                                              (df_co_occur_new["Mutation"] == "A>G") &
                                                              (df_co_occur_new["5`_ADAR_Preference"] == "Low"),
                                                              "non-ADAR A>G Synonymous",
                                                              df_co_occur_new["Stretch_Type_Non-Synonymous"])
    df_co_occur_new["Stretch_Type_Non-Synonymous"] = np.where(df_co_occur_new["Stretch_Type_Non-Synonymous"] == False,
                                                               "Stretch_Synonymous",
                                                              df_co_occur_new["Stretch_Type_Non-Synonymous"])
    df_co_occur_new["Stretch_Type_Non-Synonymous"] = np.where(df_co_occur_new["Stretch_Type_Non-Synonymous"] == True,
                                                               "Stretch_syn_non-syn",
                                                              df_co_occur_new["Stretch_Type_Non-Synonymous"])
    df_co_occur_new["Frequency"] = df_co_occur_new["Frequency"].fillna("0")
    df_co_occur_new["Stretch_Type_Non-Synonymous"] = df_co_occur_new["Stretch_Type_Non-Synonymous"].fillna("Else")
    df_co_occur_new = df_co_occur_new[df_co_occur_new["label"] != "RNA Control\nPrimer ID"]
    df_co_occur_new = df_co_occur_new[df_co_occur_new["label"] != "RNA Control\nRND"]
    df_co_occur_new = df_co_occur_new[df_co_occur_new["Read_count"] > min_coverage]
    df_co_occur_new["Frequency"] = df_co_occur_new["Frequency"].astype(float)
    df_co_occur_new.to_csv(output_dir + "/all_co_occur_grouped_adar_preferences.csv", sep=",", encoding='utf-8')
    df_co_occur_new.to_pickle(output_dir + "/all_co_occur_grouped_adar_preferences.pkl")
    return df_co_occur_new


def plot_regression(df_co_occur_new, output_dir):
    replica_lst = [1, 2, 3]
    for replica in replica_lst:
        df_co_occur_new_replica = df_co_occur_new[df_co_occur_new["replica"] == replica]

        df_co_occur_new_replica = df_co_occur_new_replica[df_co_occur_new_replica["Stretch_Type_Non-Synonymous"] != "Else"]
        data_filter_grouped_high = df_co_occur_new_replica[df_co_occur_new_replica["5`_ADAR_Preference"] == "High"]
        data_filter_grouped_high.to_csv(output_dir + "/data_filter_grouped_high_%s.csv" % str(replica), sep=",",
                                        encoding='utf-8')

        data_reg_syn = data_filter_grouped_high[
            data_filter_grouped_high["Stretch_Type_Non-Synonymous"] == "Stretch_Synonymous"]
        data_reg_non_stretch = data_filter_grouped_high[
            data_filter_grouped_high["Stretch_Type_Non-Synonymous"] == "Non-Stretch_Synonymous"]
        data_reg_non_syn = data_filter_grouped_high[
            data_filter_grouped_high["Stretch_Type_Non-Synonymous"] == "Stretch_syn_non-syn"]

        slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(data_reg_syn['passage'],
                                                                            data_reg_syn
                                                                            ['Frequency'])
        slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(data_reg_non_stretch['passage'],
                                                                            data_reg_non_stretch[
                                                                                'Frequency'])
        slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(data_reg_non_syn['passage'],
                                                                            data_reg_non_syn[
                                                                                'Frequency'])
        stretch_order = ["Stretch_Synonymous", "Non-Stretch_Synonymous", "Stretch_syn_non-syn"]

        reg_plot = sns.lmplot(x="passage", y="Frequency", data=data_filter_grouped_high,
                              col="Stretch_Type_Non-Synonymous",
                              fit_reg=True, legend=True, line_kws={'label': "Linear Reg"},
                              hue="Stretch_Type_Non-Synonymous", hue_order=stretch_order, col_order=stretch_order,
                              robust=True, truncate=True)
        label_line_1 = "y={0:.3g}x+{1:.3g}  pval={2:.3g} r-squared={3:.3g}".format(slope1, intercept1, p_value1,
                                                                                   (r_value1 ** 2))
        label_line_2 = "y={0:.3g}x+{1:.3g}  pval={2:.3g} r-squared={3:.3g}".format(slope2, intercept2, p_value2,
                                                                                   (r_value2 ** 2))
        label_line_3 = "y={0:.3g}x+{1:.3g}  pval={2:.3g} r-squared={3:.3g}".format(slope3, intercept3, p_value3,
                                                                                   (r_value3 ** 2))
        reg_plot.fig.subplots_adjust(wspace=.02)
        ax = reg_plot.axes[0, 0]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_1)

        ax = reg_plot.axes[0, 1]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_2)

        ax = reg_plot.axes[0, 2]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_3)

        reg_plot.set(ylim=(0, 0.007))

        reg_plot.savefig(output_dir + "/lmplot_replica#%s.png" % str(replica), dpi=300)
        plt.close()

def plot_regression_stretches_adar_preferences(df_co_occur_new, output_dir):
    """

    :param df_co_occur_new: mereged pd.DataFrame with stretches and ADAR preferences and for each position
    :param output_dir: output directory path, all files will be stored there
    :return: plot of ADAR preferences and stretches
    """
    replica_lst = [1, 2, 3]
    for replica in replica_lst:
        df_co_occur_new_replica = df_co_occur_new[df_co_occur_new["replica"] == replica]
        df_co_occur_new_replica["no_variants"] = df_co_occur_new_replica["Frequency"] * df_co_occur_new["Read_count"]
        df_co_occur_new_replica["freq_and_weight"] = list(zip(df_co_occur_new_replica.no_variants, df_co_occur_new_replica.Read_count))
        data_filter_grouped = df_co_occur_new_replica.groupby(["passage", "Stretch_Type_Non-Synonymous", "5`_ADAR_Preference"])[
            "freq_and_weight"].agg(
            lambda x: weighted_varaint(x))
        data_filter_grouped = data_filter_grouped.reset_index()
        data_filter_grouped = data_filter_grouped.rename(columns={"freq_and_weight": "Frequency"})
        data_filter_grouped["Frequency"] = data_filter_grouped["Frequency"].astype(float)

        data_reg_syn = data_filter_grouped[data_filter_grouped["Stretch_Type_Non-Synonymous"] == "Stretch_Synonymous"]
        data_reg_non_stretch = data_filter_grouped[data_filter_grouped["Stretch_Type_Non-Synonymous"] == "Non-Stretch_Synonymous"]
        data_reg_non_syn = data_filter_grouped[data_filter_grouped["Stretch_Type_Non-Synonymous"] == "Stretch_syn_non-syn"]
        data_reg_non_adar_ag = data_filter_grouped[data_filter_grouped["Stretch_Type_Non-Synonymous"] ==
                                                   "non-ADAR A>G Synonymous"]

        data_reg_else = data_filter_grouped[data_filter_grouped["Stretch_Type_Non-Synonymous"] ==
                                                   "Else"]

        slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(data_reg_syn['passage'],
                                                                            data_reg_syn
                                                                            ['Frequency'])
        slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(data_reg_non_stretch['passage'],
                                                                            data_reg_non_stretch[
                                                                                'Frequency'])
        slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(data_reg_non_syn['passage'],
                                                                            data_reg_non_syn[
                                                                                'Frequency'])
        slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(data_reg_non_adar_ag['passage'],
                                                                                      data_reg_non_adar_ag[
                                                                                          'Frequency'])

        stretch_order_adar = ["Stretch_Synonymous", "Non-Stretch_Synonymous", "Stretch_syn_non-syn", "non-ADAR A>G Synonymous"]
        adar_plot = sns.lmplot(x="passage", y="Frequency", data=data_filter_grouped, col="Stretch_Type_Non-Synonymous",
                                fit_reg=True, legend=True, line_kws={'label': "Linear Reg"},
                              hue="5`_ADAR_Preference", hue_order=["High", "Intermediate", "Low", "None"],
                               col_order=stretch_order_adar)
        label_line_1 = "y={0:.3g}x+{1:.3g} pval={2:.3g} r-squared={3:.3g}".format(slope1, intercept1, p_value1,
                                                                                 r_value1**2)
        label_line_2 = "y={0:.3g}x+{1:.3g} pval={2:.3g} r-squared={3:.3g}".format(slope2, intercept2, p_value2,
                                                                                 r_value2**2)
        label_line_3 = "y={0:.3g}x+{1:.3g} pval={2:.3g} r-squared={3:.3g}".format(slope3, intercept3, p_value3,
                                                                                 r_value3**2)
        label_line_4 = "y={0:.3g}x+{1:.3g} pval={2:.3g} r-squared={3:.3g}".format(slope4, intercept4, p_value4,
                                                                                 r_value4**2)

        adar_plot.fig.subplots_adjust(wspace=.02)
        ax = adar_plot.axes[0, 0]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_1)

        ax = adar_plot.axes[0, 1]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_2)

        ax = adar_plot.axes[0, 2]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_3)

        ax = adar_plot.axes[0, 3]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_4)
        adar_plot.savefig(output_dir + "/new_adar_prefe_lmplot%s.png" % str(replica), dpi=300)
        plt.close()


def plot_capsid_free_plot(df_regression, output_dir):
    g1 = sns.catplot(x="RNA", y="Stretch_Freq", data=df_regression, hue="5`_ADAR_Preference", order=["Capsid", "Free"],
                     hue_order=["High", "Intermidate", "Low", "None"], kind="strip", col="replica")
    g1.set(yscale='log')
    g1.savefig(output_dir + "/capsid_adar_pref.png", dpi=300)
    plt.close()


def main():
    """RV"""
    # input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages"
    # experiment = "passages"
    # output_dir = input_dir + "/Co_occur_all_%s" % experiment
    # prefix = "/p*"
    # min_coverage = 5000
    # virus = "RVB14"
    # date = "20201012"
    # q = "q38"
    # """1"""
    # region_lst = [629, 835, 1621, 2329, 3196, 3634, 3925, 4915, 5170, 5239, 5785, 7165]
    # data_all_passages = creating_co_occur_passages_df(input_dir, prefix, date, q, region_lst)
    #
    # """2"""
    # data_all_passages = pd.read_pickle(output_dir + "/all_co_occur_protein.pkl")
    # df_co_occur_new = grouped_co_occur(data_all_passages, input_dir, experiment, output_dir)
    #
    # g1 = sns.boxplot(x="Stretch_Type_Non-Synonymous", y="Stretch_Freq", data=df_co_occur_new,
    #                  order=[False, True])
    # old_statannot.add_stat_annotation(g1, data=df_co_occur_new, x="Stretch_Type_Non-Synonymous", y="Stretch_Freq",
    #                     boxPairList=[(False, True)], test='Mann-Whitney', textFormat='star',
    #                                   loc='outside', verbose=2,
    #                                   order=[False, True])
    # g1.set(xticklabels=("Synonymous\nStretches", "Stretches\nWith Non-Synonymous"))
    # g1.set(xlabel="")
    # plt.tight_layout()
    # plt.savefig(output_dir + "/Frequencies_of_Stertch_Type_Non-Synonymous", dpi=300)
    # plt.close()
    # df_adar_path = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/" \
    #                "inosine_predict_context/data_filter.csv"
    # """3"""
    # df_co_occur_new = pd.read_pickle(output_dir + "/all_co_occur_grouped.pkl")
    # df_regression = regression_stretches(df_co_occur_new, output_dir, df_adar_path, min_coverage)
    # plot_regression(df_regression, output_dir)
    # """4"""
    # df_regression = pd.read_pickle(output_dir + "/all_co_occur_grouped_adar_preferences.pkl")
    # plot_regression_stretches_adar_preferences(df_regression, output_dir)


    """RV-Capsid_Free"""
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/capsid"
    experiment = "capsid"
    output_dir = input_dir + "/Co_occur_all_%s" % experiment
    prefix = "/*_3*"
    min_coverage = 5000
    virus = "RVB14"
    date = "20201012"
    q = "q38"
    """1"""
    region_lst = [629, 835, 1621, 2329, 3196, 3634, 3925, 4915, 5170, 5239, 5785, 7165]
    data_all_passages = creating_co_occur_passages_df(input_dir, experiment, prefix, date, q, region_lst)

    """2"""
    # data_all_passages = pd.read_pickle(output_dir + "/all_co_occur_protein.pkl")
    df_co_occur_new = grouped_co_occur(data_all_passages, input_dir, experiment, output_dir, q)
    #
    # Plots
    df_co_occur_new = pd.read_pickle(output_dir + "/all_co_occur_grouped.pkl")
    g1 = sns.catplot(x="RNA", y="Stretch_Freq", data=df_co_occur_new, order=["Capsid", "Free"], kind="strip",
                     col="replica")
    g1.set(yscale='log')
    g1.set(ylim=(10**-5, 10**-1))
    plt.tight_layout()
    plt.savefig(output_dir + "/Frequencies_of_Stertch", dpi=300)
    plt.close()
    stretch_g1 = sns.catplot(x="Stretch_Type_Non-Synonymous", y="Stretch_Freq", data=df_co_occur_new,
                             order=[False, True], hue="RNA", hue_order=["Capsid", "Free"], kind="strip", col="replica")
    stretch_g1.set(xticklabels=("Synonymous\nStretches", "Stretches\nWith Non-Synonymous"))
    stretch_g1.set(xlabel="")
    plt.tight_layout()
    plt.savefig(output_dir + "/Frequencies_of_Stertch_Type_Non-Synonymous_replicas", dpi=300)

    """3"""
    df_adar_path = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/capsid/" \
                   "inosine_predict_context/data_filter.csv"

    df_co_occur_new = pd.read_pickle(output_dir + "/all_co_occur_grouped.pkl")
    df_regression = regression_stretches(df_co_occur_new, output_dir, df_adar_path, experiment, min_coverage)
    """4"""
    # df_regression = pd.read_pickle(output_dir + "/all_co_occur_grouped_adar_preferences.pkl")
    plot_capsid_free_plot(df_regression, output_dir)

    """RV-Patients"""
    # input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/patients"
    # output_dir = input_dir + "/Co_occur_patients"
    # prefix = "/*"
    # min_coverage = 5000
    # virus = "RVB14"
    # date = "20201017"
    # q = "q30_consensusX5"
    #
    # region_lst = [503, 709, 1495, 2197, 3094, 3523, 3808, 4771, 5005, 5068, 5617, 6728]
    # data_all_passages = creating_co_occur_passages_df(input_dir, prefix, date, q, region_lst)
    # # data_all_passages = pd.read_pickle(output_dir + "/all_co_occur_protein.pkl")
    # df_co_occur_new = grouped_co_occur(data_all_passages, input_dir, output_dir, q=q.split("_")[0])
    # df_adar_path = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/" \
    #                "inosine_predict_context/data_filter.csv"
    #
    # df_co_occur_new = pd.read_pickle(output_dir + "/all_co_occur_grouped.pkl")
    # df_regression = regression_stretches(df_co_occur_new, output_dir, df_adar_path, min_coverage)
    # # df_regression = pd.read_pickle(output_dir + "/all_co_occur_grouped_adar_preferences.pkl")
    # df_regression_adar = regression_stretches_adar_preferences(df_regression, output_dir)
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