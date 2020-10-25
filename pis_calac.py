import os
import warnings
import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns;
sns.set_context("poster")
import sys


def pis_calc(data, pivot_cols=[], interval = (0, sys.maxsize)): # start_pos=0, end_pos= sys.maxsize
    """
    Calculates PI diversity per position, than calculates mean per group according to pivot_vols. Assumes data is not indexed.
    :param data:
    :param pivot_cols:
    :param min_read_count:
    :return:
    """

    def pairwise_differences_proportion(row):
        if row["Minor"] == 0:
            return 0
        total_part = (row["Total"] * (row["Total"] - 1))
        numerator = total_part - ((row["Major"] * (row["Major"] - 1)) + (row["Minor"] * (row["Minor"] - 1)))
        denominator = total_part
        return numerator * 1.0 / denominator
    # Filters
    # TODO- extract this filter too\ insert all others here
    # choose interval
    filtered_data = data[(data["Pos"] >= interval[0]) & (data["Pos"] <= interval[1])]
    if filtered_data.empty:
        warnings.warn("No relevant data after filtering. Skipping")
        return None
    filtered_data['counts_for_position'] = np.round(filtered_data['Read_count'] * filtered_data['Frequency'])
    selecting_cols = pivot_cols[:]
    selecting_cols.extend(["Pos", "Base", "counts_for_position", "Rank"])
    filtered_data = filtered_data[selecting_cols]
    filtered_data["Rank"] = np.where(filtered_data["Rank"] == 0, "Major", "Minor")
    group_cols = pivot_cols[:]
    group_cols.extend(["Pos", "Rank"])
    # selecting max on cfp, per Pos\Rank (Major\minor?)- than will be graded per Pos, and summed
    # TODO: use all minor variants (not only max)
    filtered_data = filtered_data.groupby(group_cols)['counts_for_position'].aggregate(max).unstack().reset_index()
    if 'Minor' not in filtered_data.columns:
        return 0
    filtered_data["Total"] = filtered_data["Major"] + filtered_data["Minor"]
    filtered_data["pdp"] = filtered_data.apply(lambda row: pairwise_differences_proportion(row), axis=1)
    if any(pivot_cols):
        bysample_diversity = filtered_data.groupby(pivot_cols)['pdp'].agg(['count', 'sum']).reset_index()
        bysample_diversity["Pi"] = bysample_diversity["sum"] * 1.0 / bysample_diversity["count"]
        output_cols = pivot_cols[:]
        output_cols.append("Pi")
        pis = bysample_diversity[output_cols]
    else:
        pis = filtered_data['pdp'].mean()
    return pis