#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""
import glob
import sys, argparse
import pandas as pd
import numpy as np
from textwrap import wrap
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(font_scale=1.2)
sns.set_style("ticks")
sns.despine()

def insert_to_df(df, row):
    insert_loc = df.index.max()
    if pd.isna(insert_loc):
        df.loc[0] = row
    else:
        df.loc[insert_loc + 1] = row
    return df


def string_delimiter(string):
    string_lst = wrap(string, 2)
    return string_lst


def AG_read_counter(data):
    data = pd.DataFrame(data.columns.values[None, :], columns=data.columns).append(data).reset_index(drop=True)
    data.columns = range(data.shape[1])
    columns = ["read", "start_pos", "end_pos", "start_read", "end_read", "direction", "read_len", "blast_out"]
    df = pd.DataFrame(data.values, columns=columns)
    df1 = df[["read", "blast_out"]]
    df1["new_blast_out"] = df1["blast_out"].str.replace('\d+', '')
    df1["new_blast_out"] = df1["new_blast_out"].apply(lambda x: string_delimiter(x))
    df1["new_blast_out_count"] = df1["new_blast_out"].apply(lambda x: x.count("AG"))
    df1["in_stretch"] = np.where(df1["new_blast_out_count"] >= 3, True, False)
    df_grouped = df1.groupby("in_stretch", as_index=False).count()
    df_grouped = df_grouped[["in_stretch", "read"]]
    df_grouped = df_grouped.rename(columns={"read": "count"})
    df_grouped = df_grouped.set_index("in_stretch")
    return df1, df_grouped


def my_crosstab(df_control_grouped, df_grouped, key, value):
    crosstab_df = pd.merge(df_control_grouped, df_grouped, on=df_grouped.index)
    crosstab_df = crosstab_df.set_index("key_0")
    crosstab_df = crosstab_df.rename(columns={"count_x": "Control", "count_y": key})
    crosstab_df = crosstab_df.apply(lambda x: x/value)
    oddsratio, pval = stats.fisher_exact(crosstab_df, alternative="two-sided")
    insert_to_df(crosstab_df, [None, oddsratio])
    insert_to_df(crosstab_df, [None, pval])
    crosstab_df = crosstab_df.rename(index={2: "odssratio", 3: "pval", True: "No._of_reads_with_stretch_A>G", False: "No._of_reads_without_stretch_A>G"})
    # print("data of:{0}".format(key))
    # print(crosstab_df.to_string())
    # print("oddsratio={0}    pval={1}".format(oddsratio, pval))
    return crosstab_df


def create_crosstab_df(input_dir, prefix, data_dict):
    data_control = pd.read_table(input_dir + "/IVT_3_Control/{0}".format(prefix), sep="\t")
    df_control, df_control_grouped = AG_read_counter(data_control)
    crosstab_lst = []
    for key, value in data_dict.items():
        df, df_grouped = AG_read_counter(value[0])
        df_grouped.to_pickle(input_dir + "/{0}/20201012_q38/grouped.pkl".format(key))
        crosstab_df = my_crosstab(df_control_grouped, df_grouped, key, value[1])
        crosstab_df.to_pickle(input_dir + "/{0}/20201012_q38/corsstab_df.pkl".format(key))
        crosstab_df.to_csv(input_dir + "/{0}/20201012_q38/corsstab_df.csv".format(key), sep=",")
        crosstab_lst.append(crosstab_df)
    return crosstab_lst

def main():
    # input_dir = "/Users/odedkushnir/Google Drive/Studies/PhD/Stretch_analysis"
    input_dir = "C:/Users/odedku/Stretch_analysis"
    prefix = "20201012_q38/all_parts.blast"
    p2_1 = pd.read_table(input_dir + "/p2_1/{0}".format(prefix), sep="\t")
    p2_2 = pd.read_table(input_dir + "/p2_2/{0}".format(prefix), sep="\t")
    p5_1 = pd.read_table(input_dir + "/p5_1/{0}".format(prefix), sep="\t")
    p5_2 = pd.read_table(input_dir + "/p5_2/{0}".format(prefix), sep="\t")
    p8_1 = pd.read_table(input_dir + "/p8_1/{0}".format(prefix), sep="\t")
    p8_2 = pd.read_table(input_dir + "/p8_2/{0}".format(prefix), sep="\t")
    p10_1 = pd.read_table(input_dir + "/p10_1/{0}".format(prefix), sep="\t")
    p10_2 = pd.read_table(input_dir + "/p10_2/{0}".format(prefix), sep="\t")
    p12_1 = pd.read_table(input_dir + "/p12_1/{0}".format(prefix), sep="\t")
    p12_2 = pd.read_table(input_dir + "/p12_2/{0}".format(prefix), sep="\t")
    # Dictionary of passage and number of PrimerID
    data_dict = {"p2_1": [p2_1, 10000], "p2_2": [p2_2, 10000], "p5_1": [p5_1, 10000], "p5_2": [p5_2, 10000],
                 "p8_1": [p8_1, 10000], "p8_2": [p8_2, 10000], "p10_1": [p10_1, 10000], "p10_2": [p10_2, 10000],
                 "p12_1": [p12_1, 10000], "p12_2": [p12_2, 10000]}

    """NOT from memory"""
    create_crosstab_df(input_dir, prefix, data_dict)

    """from memory"""
    passage_lst = glob.glob(input_dir + "/p*")
    crosstab_lst = []
    for passage in passage_lst:
        passage_num = passage.split("\\")[-1]
        crosstab_df = pd.read_pickle(input_dir + "/{0}/20201012_q38/corsstab_df.pkl".format(passage_num))
        crosstab_lst.append(crosstab_df)
    crosstab_df_all = pd.concat(crosstab_lst, axis=1)
    crosstab_df_all = crosstab_df_all[["Control", "p2_1", "p2_2", "p5_1", "p5_2", "p8_1", "p8_2", "p10_1", "p10_2", "p12_1", "p12_2"]]
    crosstab_df_all = crosstab_df_all.iloc[0:4, 9:]
    crosstab_df_all = crosstab_df_all.transpose()
    crosstab_df_all["Stretch_percentage"] = crosstab_df_all["No._of_reads_with_stretch_A>G"] / \
                                           (crosstab_df_all["No._of_reads_with_stretch_A>G"]+crosstab_df_all["No._of_reads_without_stretch_A>G"])*100
    crosstab_df_all.reset_index(inplace=True, drop=False)
    crosstab_df_all["passage"] = np.where(crosstab_df_all["index"] != "Control", crosstab_df_all.apply(lambda x: str(x["index"]).split("_")[0].split("p")[-1], axis=1), 0)
    crosstab_df_all["replica"] = np.where(crosstab_df_all["index"] != "Control", crosstab_df_all.apply(lambda x: str(x["index"]).split("_")[-1], axis=1), 1)
    crosstab_df_all["passage"] = crosstab_df_all["passage"].astype(int)
    crosstab_df_all.to_csv(input_dir + "/crosstab_df_all.csv", sep=",")
    mean_crosstab_df_all = crosstab_df_all.groupby("passage", as_index=False).mean()
    mean_crosstab_df_all["sem"] = crosstab_df_all.groupby("passage", as_index=False).sem()["Stretch_percentage"]
    mean_crosstab_df_all.to_csv(input_dir + "/mean_crosstab_df_all.csv", sep=",")

    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(crosstab_df_all['passage'],
                                                                        crosstab_df_all[
                                                                            'Stretch_percentage'])
    fig1 = sns.lmplot(x="passage", y="Stretch_percentage", data=crosstab_df_all, fit_reg=True, line_kws={'label': "Linear Reg"},)
    fig1.set(xlabel="Passage", ylabel="Stretch Percentage [%]", xlim=(0, 12))
    ax = fig1.axes[0, 0]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    label_line_1 = "y={0:.3g}x+{1:.3g}\nstderr={2:.3g} Rsq={3:.3g}".format(slope1, intercept1, std_err1, r_value1**2)
    L_labels[0].set_text(label_line_1)
    plt.show()

    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(mean_crosstab_df_all['passage'],
                                                                        mean_crosstab_df_all[
                                                                            'Stretch_percentage'])
    fig2 = sns.lmplot(x="passage", y="Stretch_percentage", data=mean_crosstab_df_all, fit_reg=True, line_kws={'label': "Linear Reg"},)
    fig2.set(xlabel="Passage", ylabel="Stretch Percentage [%]", xlim=(0, 12))
    ax = fig2.axes[0, 0]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    label_line_2 = "y={0:.3g}x+{1:.3g}\nstderr={2:.3g} Rsq={3:.3g}".format(slope2, intercept2, std_err2, r_value2**2)
    L_labels[0].set_text(label_line_2)
    plt.show()


if __name__ == "__main__":
    main()

#     data = data.rename(columns={"ref_pos": "Pos", "base": "Base"})
# data = data[data["Pos"] != "ref_pos"]
# data["Pos"] = data["Pos"].astype(int)
#
# consensus_data = pd.read_table(input_dir + "/IVT_3_Control/20201012_q38/IVT-3-Control.merged.freqs", sep="\t")
# consensus_data = consensus_data[consensus_data["Rank"] == 0]
# consensus_data = consensus_data[["Pos", "Ref"]]
# consensus_data = consensus_data[consensus_data["Ref"] != "-"]
# consensus_data["Pos"] = consensus_data["Pos"].astype(int)
# df = pd.merge(data, consensus_data, on="Pos", how="left")
# df["Mutation"] = df["Ref"] + ">" + df["Base"]
# df_ag = df[df["Mutation"] == "A>G"]
# df_ag["next_pos"] = df_ag.shift(-1)["Pos"]
# df_ag["next_read"] = df_ag.shift(-1)["read"]
# df_ag["delta"] = df_ag["next_pos"] - df["Pos"]
#
# df_ag["in_stertch"] = np.where(df_ag["delta"] < 4, True, False)
# df_ag["in_stertch"] = np.where(df_ag["delta"] <= 0, False, df_ag["in_stertch"])
# df_ag["same_read"] = np.where(df_ag["read"] == df_ag["next_read"], True, False)
# df_grouped = df_ag.groupby("in_stertch", as_index=False).count()