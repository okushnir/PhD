#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import sys, argparse
import pandas as pd
import numpy as np
from string import digits

def pos_next_line(data, i):
    pos_next_line = data["Pos"].iloc[i+1]
    return pos_next_line

def main():

    input_dir = "/Users/odedkushnir/Google Drive/Studies/PhD/test"
    data = pd.read_table(input_dir + "/p2_1/20201012_q38/all_parts.blast", sep="\t")

    data = pd.DataFrame(data.columns.values[None, :], columns=data.columns).append(data).reset_index(drop=True)
    data.columns = range(data.shape[1])
    columns = ["read", "start_pos", "end_pos", "start_read", "end_read", "direction", "read_len", "blast_out"]
    df = pd.DataFrame(data.values, columns=columns)
    df1 = df[["read", "blast_out"]]
    df1["new_blast_out"] = df1["blast_out"].str.replace('\d+', '')
    # TODO:
    # count the AG in new_blast_out
    # if >3 True
    # count True/False"""

    # df_grouped = df1.groupby("new_blast_out", as_index=False).count()


    print(df1.to_string())

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