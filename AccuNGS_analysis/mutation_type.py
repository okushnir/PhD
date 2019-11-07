
#! /usr/local/python/anaconda_python-3.6.1

"""
@Author: odedkushnir

"""

import sys
import os
import re
import glob
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
import cirseq_utilities
import subprocess
from time import sleep
from os import listdir
from os.path import isfile, join
from optparse import OptionParser
import numpy as np
import pbs_jobs
import make_reference_from_consensus

def add_mutation_type_to_df(data, ncbi_id):
#this part creates required data to work with (basically the consensus sequence that will be the template for the "wildtype" amino acid)
    min_read_count = 10000

    data['counts_for_position'] = np.round(data['Read_count']*data['Freq'])
    data = data[data['Pos'] == np.round(data['Pos'])] #remove insertions
    data['Pos'] = data[['Pos']].apply(pd.to_numeric)
    data = data[(data['Read_count'] > min_read_count)] #| (data["Read_count"] > 5000)] &(data["control"]==label40806)


    # data['Base'].replace('T', 'U', inplace=True)
    # data['Ref'].replace('T', 'U', inplace=True)
    """Make Consensus column in dataframe"""
    consensus_data = data[['Pos', 'Base', "Read_count", "Freq", 'counts_for_position']][data['Rank'] == 0]
    consensus_data = consensus_data[consensus_data['Pos'] == np.round(consensus_data['Pos'])]  # removes insertions
    consensus_data['Pos'] = consensus_data[['Pos']].apply(pd.to_numeric)
    consensus_data["Context"] = consensus_data.shift(1)['Base'] + consensus_data['Base'] + consensus_data.shift(-1)['Base']
    consensus_data = consensus_data.rename(columns={"Base": "Consensus"})

    data = pd.merge(data, consensus_data[['Pos', 'Consensus', 'Context']], on=["Pos"])
    data['BaseToBase'] = data['Consensus'] + data['Base']
    data['Mutation_class'] = np.where(data["Rank"] == 0, "self",
        np.where(data['BaseToBase'] == 'GA', 'transition', np.where(data['BaseToBase'] == 'AG', 'transition',
                                                                       np.where(data['BaseToBase'] == 'CT',
                                                                                'transition',
                                                                                np.where(data['BaseToBase'] == 'TC',
                                                                                         'transition', 'transversion')))
                 ))
    
    # data["adj_freq"] = np.where(1>data["counts_for_position"], 1, data["counts_for_position"]) / data["Read_count"]
    #
    # dataToMutate = data[data["Rank"] == 0]
    # dataToMutate["adj_freq"] = 1-dataToMutate["adj_freq"]
    # dataToMutate["Mutation_type"] = "All"

    
    # transitions = data
    # #transitions = transitions[transitions["Base"] != '-'] #remove deletions
    # transitions = transitions[transitions["Mutation_class"] != "transversion"]
    # grouped = transitions.groupby(["Pos"], as_index=False)["counts_for_position"].agg("sum")
    # transitions = pd.merge(transitions,grouped, on=["Pos"])
    # transitions["adjusted_freq"] = (transitions["counts_for_position_x"] / transitions["counts_for_position_y"]).apply(lambda x : 5*10**-6 if x == 0 else x)

#this part creates the actual translation
    start_pos, end_pos = cirseq_utilities.find_coding_region(ncbi_id)
    start_pos_mod3 = start_pos % 3

    noncoding_data = data[(data["Base"]!="-") & (data["Pos"] < start_pos) & (data["Pos"] > end_pos)]

    data = data[(data["Base"] != "-") & (data["Pos"] >= start_pos) & (data["Pos"] <= end_pos)]
    data["prevBase"] = data.shift(4)['Consensus']
    data["nextBase"] = data.shift(-4)['Consensus']
    data["prev2bases"] = data.shift(8)['Consensus'] + data.shift(4)['Consensus']
    data["next2bases"] = data.shift(-4)['Consensus'] + data.shift(-8)['Consensus']

    data.to_csv("/Users/odedkushnir/Projects/fitness/test/test.csv", sep=',', encoding='utf-8')
    data["Consensus_codon"] = np.where(data["Pos"] % 3 == start_pos_mod3, data["Consensus"] + data["next2bases"],
                                                np.where(data["Pos"] % 3 == start_pos_mod3 + 1, data["Context"],
                                                         data["prev2bases"] + data["Consensus"]))
    data["Mutated_codon"] = np.where(data["Pos"] % 3 == start_pos_mod3, data["Base"] + data["next2bases"]
                                     , np.where(data["Pos"]% 3 == start_pos_mod3 + 1, data["prevBase"] + data["Base"] +
                                                data["nextBase"], data["prev2bases"] + data["Base"]))
    data.dropna(subset=["Consensus_codon", "Mutated_codon"], inplace=True)
    # data.to_csv("/Users/odedkushnir/Projects/fitness/test/test2.csv", sep=',', encoding='utf-8')
    data["Consensus_aa"] =data["Consensus_codon"].apply(lambda x: Seq(x).translate()[0] if "-" not in x else "-")
    data["Mutated_aa"] = data["Mutated_codon"].apply(lambda x: Seq(x).translate()[0] if "-" not in x else "-")
    data["Type"] = np.where(data["Mutated_aa"] == "*", "Premature Stop Codon", np.where(data["Mutated_aa"] == data["Consensus_aa"],"Synonymous","Non-Synonymous"))

    return data


def main():
    data = "/Users/odedkushnir/Projects/fitness/test/ERR2352296.freqs"
    data = cirseq_utilities.freqs_to_dataframe(data)
    data = add_mutation_type_to_df(data, "MH732737")
    data.to_csv("/Users/odedkushnir/Projects/fitness/test/data_aa.csv", sep=',', encoding='utf-8')
    print(data)


if __name__ == "__main__":
    main()