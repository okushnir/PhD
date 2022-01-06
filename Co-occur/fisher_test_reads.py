#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import sys, argparse
import pandas as pd
import numpy as np
from textwrap import wrap
import scipy.stats as stats


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

def main():
    # input_dir = "/Users/odedkushnir/Google Drive/Studies/PhD/test"
    input_dir = "C:/Users/odedku/test"
    data_control = pd.read_table(input_dir + "/IVT_3_Control/20201012_q38/all_parts.blast", sep="\t")
    data_p2_1 = pd.read_table(input_dir + "/p2_1/20201012_q38/all_parts.blast", sep="\t")
    df_control, df_control_grouped = AG_read_counter(data_control)
    df_p2_1, df_p2_1_grouped = AG_read_counter(data_p2_1)
    # print(df_control.to_string())
    # print(df_p2_1.to_string())
    # crossta_df = pd.crosstab(index=df_control["in_stretch"], columns=df_p2_1["in_stretch"])
    crosstab_df = pd.merge(df_control_grouped, df_p2_1_grouped, on=df_p2_1_grouped.index)
    crosstab_df = crosstab_df.set_index("key_0")
    crosstab_df = crosstab_df.rename(columns={"count_x": "Control", "count_y": "Test"})
    oddsratio, pval = stats.fisher_exact(crosstab_df, alternative="two-sided")
    print(crosstab_df.to_string())
    print("oddsratio={0}    pval={1}".format(oddsratio, pval))

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