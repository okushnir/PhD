"""
@Author: odedkushnir

"""

import sys
import glob
import pandas as pd
import numpy as np

import shutil
import os


def add_mutation(freqs_file, out_dir):
    """
    :param freqs_file: pandas DataFrame after find_mutation_type function
    :return: csv dataframe
    """
    data = pd.read_table(freqs_file)
    data.reset_index(drop=True, inplace=True)
    flag = '-' in data.Base.values
    if flag is True:
        data = data[data.Ref != '-']
        data = data[data.Base != '-']
        data.reset_index(drop=True, inplace=True)
        # data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
        # data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
        # raise Exception("This script does not support freqs file with deletions, for support please contact Maoz ;)"
    data['Base'].replace('T', 'U', inplace=True)
    data['Ref'].replace('T', 'U', inplace=True)
    min_read_count = 10000
    data = data[data['Pos'] == np.round(data['Pos'])]  # remove insertions
    data['Pos'] = data[['Pos']].apply(pd.to_numeric)
    data = data[data['Read_count'] > min_read_count]
    data['mutation_type'] = data['Ref'] + data['Base']
    data = data[data['Ref'] != data['Base']]
    data = data[data["Base"] != "-"]
    data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
    data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
    data["Mutation"] = data["Ref"] + "->" + data["Base"]
    data.to_csv(out_dir + "RVB14p2.csv", sep=",", encoding='utf-8')
    return data


def context_mutation(freqs_file, out_dir):

    # out_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/Cirseq/RV/20190115_q30rRNA/"
    # data = pd.read_csv(out_dir + "RVB14p2.csv")

    data = pd.read_table(freqs_file)
    data['Base'].replace('T', 'U', inplace=True)
    data['Ref'].replace('T', 'U', inplace=True)

    # generate consensus of both sample and control
    consensus_data = data[['Pos', 'Base', "Read_count", "Freq"]][data['Rank'] == 0]
    consensus_data = consensus_data[consensus_data['Pos'] == np.round(consensus_data['Pos'])]  # removes insertions
    consensus_data['Next'] = consensus_data['Base'] + consensus_data.shift(-1)['Base']
    consensus_data['Prev'] = consensus_data.shift(1)['Base'] + consensus_data['Base']
    consensus_data['Pos'] = consensus_data[['Pos']].apply(pd.to_numeric)
    consensus_data["Context"] = consensus_data.shift(1)['Base'] + consensus_data['Base'] + consensus_data.shift(-1)['Base']
    consensus_data["major_read_count"] = np.round(consensus_data['Read_count'] * consensus_data['Freq'])
    consensus_data = consensus_data.rename(columns={"Base": "Consensus"})
    data = pd.merge(data, consensus_data[['Pos', 'Prev', 'Next', 'Context', 'Consensus', "major_read_count"]],
                    on=["Pos"])
    data["mutation_type"] = data["Consensus"] + data["Base"]
    data["Mutation"] = data["Ref"] + "->" + data["Base"]
    data.to_csv(out_dir + "/data_mutation.csv", sep=',', encoding='utf-8')
    mutation_for_rna = ["AG"]
    dataForPlotting = data[(data["mutation_type"].isin(mutation_for_rna)) & (data["Rank"] != 0)]
    dataForPlotting.to_csv(out_dir + "/data_XpA_by_mutation.csv", sep=',', encoding='utf-8')
    return dataForPlotting

def main():
    freqs_file = "/Volumes/STERNADILABHOME$/volume3/okushnir/Cirseq/RV/20171220_q23r2_blastn/RVB14-p2.with.mutation.type.freqs"
    out_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/Cirseq/RV/20171220_q23r2_blastn/"
    data = add_mutation(freqs_file, out_dir)
    context_data = context_mutation(freqs_file, out_dir)
if __name__ == "__main__":
    main()