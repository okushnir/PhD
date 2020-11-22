#! /usr/local/python/anaconda_python-3.6.1


import pandas as pd
import numpy as np
import time
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO
from Utilities import cirseq_utilities



def find_mutation_type(freqs_file, ncbi_id, min_read_count=None, seq_method="AccuNGS"):
    """
    This function adds Mutation type to the freqs file
    :param freqs_file:  The path of the relevant freqs file
    :param ncbi_id: NCBI number
    :param min_read_count: Default 10,0000
    :return:DataFrame of the freqs file with mutation type column, saves it into txt file
    """

    #Debug
    # freqs_file = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/RV-p11/q30_3UTR_new/RV-p11.freqs"
    # ncbi_id = "NC_001490"
    # start_pos, end_pos = (629, 7168)

    start_time = time.time()
    file_name = freqs_file
    data = freqs_to_dataframe(freqs_file)
    data.reset_index(drop=True, inplace=True)
    orign_data = data
    start_pos, end_pos = find_coding_region(ncbi_id)

    if min_read_count is None:
        min_read_count = 10000
    data = add_mutation_type_to_df(data, ncbi_id, min_read_count)

    top_data = orign_data.loc[orign_data['Pos'] < start_pos]
    top_data["counts_for_position"] = ""
    top_data["Consensus"] = ""
    top_data["Context"] = ""
    top_data["BaseToBase"] = ""
    top_data["Mutation_class"] = ""
    top_data["adj_freq"] = ""
    top_data["prevBase"] = ""
    top_data["nextBase"] = ""
    top_data["prev2bases"] = ""
    top_data["next2bases"] = ""
    top_data["Consensus_codon"] = ""
    top_data["Mutated_codon"] = ""
    top_data["Consensus_aa"] = ""
    top_data["Mutated_aa"] = ""
    top_data["Type"] = ""
    top_data['Pos'] = top_data['Pos'].astype(int)
    lower_data = orign_data.loc[orign_data['Pos'] > end_pos]
    lower_data["counts_for_position"] = ""
    lower_data["Consensus"] = ""
    lower_data["Context"] = ""
    lower_data["BaseToBase"] = ""
    lower_data["Mutation_class"] = ""
    lower_data["adj_freq"] = ""
    lower_data["prevBase"] = ""
    lower_data["nextBase"] = ""
    lower_data["prev2bases"] = ""
    lower_data["next2bases"] = ""
    lower_data["Consensus_codon"] = ""
    lower_data["Mutated_codon"] = ""
    lower_data["Consensus_aa"] = ""
    lower_data["Mutated_aa"] = ""
    lower_data["Type"] = ""
    lower_data['Pos'] = lower_data['Pos'].astype(int)
    frames = [top_data, data, lower_data]
    data_final = pd.concat(frames)
    if seq_method == "CirSeq":
        data_final = data_final[["Pos", "Base", "Freq", "Ref", "Read_count", "Rank", "Prob",
                                    "counts_for_position", "Consensus", "BaseToBase", "Mutation_class", "adj_freq",
                                    "prevBase", "nextBase", "prev2bases", "next2bases", "Consensus_codon",
                                    "Mutated_codon", "Consensus_aa", "Mutated_aa", "Type"]]
    else:
        data_final = data_final[["Pos", "Base", "Freq", "Ref", "Read_count", "Rank", "Prob",
                                        "counts_for_position", "Consensus", "BaseToBase", "Mutation_class", "adj_freq",
                                        "prevBase", "nextBase", "prev2bases", "next2bases", "Consensus_codon",
                                        "Mutated_codon", "Consensus_aa", "Mutated_aa", "Type", "pval", "Var_perc", "SNP_Profile"]]

    file_name = file_name[0:-5]
    file_name += "with.mutation.type.freqs"
    data_final.to_csv(file_name, sep='\t', encoding='utf-8')
    print("The File is ready in the folder")
    print("--- %s sec ---" % (time.time() - start_time))
    print("start_pos:%i" % start_pos)
    print("end_pos:%i" % end_pos)
    return data_final


def add_mutation_type_to_df(data, ncbi_id, min_read_count=None):
    """
    This function adds mutation types just for to the Coding Region of freqs file
    :param data: DataFrame of the freqs file
    :param ncbi_id: NCBI number
    :param min_read_count: Default 10,000
    :return:
    """


    if min_read_count is None:
        min_read_count = 10000
    # this part creates required data to work with (basically the consensus sequence that will be the template for the "wildtype" amino acid)

    data['counts_for_position'] = np.round(data['Read_count'] * data['Freq'])
    data = data[data['Pos'] == np.round(data['Pos'])]  # remove insertions
    data['Pos'] = data[['Pos']].apply(pd.to_numeric)
    data = data[(data['Read_count'] > min_read_count)]  # | (data["Read_count"] > 5000)] &(data["control"]==label40806)

    """Make Consensus column in the dataframe"""
    consensus_data = data[['Pos', 'Base', "Read_count", "Freq", 'counts_for_position']][data['Rank'] == 0]
    consensus_data = consensus_data[consensus_data['Pos'] == np.round(consensus_data['Pos'])]  # removes insertions
    consensus_data['Pos'] = consensus_data[['Pos']].apply(pd.to_numeric)

    #TO DO: write a function that checks continuity in Pos column - if the Pos is NOT continuously then ther si a frame shift!!!

    consensus_data["Context"] = consensus_data.shift(1)['Base'] + consensus_data['Base'] + consensus_data.shift(-1)[
        'Base']
    consensus_data = consensus_data.rename(columns={"Base": "Consensus"})

    data = pd.merge(data, consensus_data[['Pos', 'Consensus', 'Context']], on=["Pos"])
    data['BaseToBase'] = data['Consensus'] + data['Base']
    data['Mutation_class'] = np.where(data["Rank"] == 0, "self",
                                      np.where(data['BaseToBase'] == 'GA', 'transition',
                                               np.where(data['BaseToBase'] == 'AG', 'transition',
                                                        np.where(data['BaseToBase'] == 'CT',
                                                                 'transition',
                                                                 np.where(data['BaseToBase'] == 'TC',
                                                                          'transition', 'transversion')))
                                               ))
    # this part creates the actual translation
    start_pos, end_pos = find_coding_region(ncbi_id)
    start_pos_mod3 = start_pos % 3

    noncoding_data = data[(data["Base"] != "-") & (data["Pos"] < start_pos) & (data["Pos"] > end_pos)]

    data = data[(data["Base"] != "-") & (data["Pos"] >= start_pos) & (data["Pos"] <= end_pos)]
    data["prevBase"] = data.shift(4)['Consensus']
    data["nextBase"] = data.shift(-4)['Consensus']
    data["prev2bases"] = data.shift(8)['Consensus'] + data.shift(4)['Consensus']
    data["next2bases"] = data.shift(-4)['Consensus'] + data.shift(-8)['Consensus']

    # data.to_csv("/Users/odedkushnir/Projects/fitness/test/SRR6955546/test.csv", sep=',', encoding='utf-8')
    if (start_pos_mod3 == 0) | (start_pos_mod3 == 1):
        data["Consensus_codon"] = np.where(data["Pos"] % 3 == start_pos_mod3, data["Consensus"] + data["next2bases"],
                                           np.where(data["Pos"] % 3 == start_pos_mod3 + 1, data["Context"], data["prev2bases"] + data["Consensus"]))

        data["Mutated_codon"] = np.where(data["Pos"] % 3 == start_pos_mod3, data["Base"] + data["next2bases"],
                                           np.where(data["Pos"] % 3 == start_pos_mod3 + 1, data["prevBase"] + data["Base"] + data["nextBase"], data["prev2bases"] + data["Base"]))
    elif start_pos_mod3 == 2:
        data["Consensus_codon"] = np.where(data["Pos"] % 3 == start_pos_mod3, data["Consensus"] + data["next2bases"],
                                           np.where(data["Pos"] % 3 == start_pos_mod3 + 1,
                                                    data["prev2bases"] + data["Consensus"],
                                                    np.where(data["Pos"] % 3 == 0, data["Context"],
                                                             np.where(data["Pos"] % 3 == 1,
                                                                      data["prev2bases"] + data["Consensus"],
                                                                      data["Context"]))))
        data["Mutated_codon"] = np.where(data["Pos"] % 3 == start_pos_mod3, data["Base"] + data["next2bases"],
                                           np.where(data["Pos"] % 3 == start_pos_mod3 + 1,
                                                    data["prev2bases"] + data["Base"],
                                                    np.where(data["Pos"] % 3 == 0, data["prevBase"] + data["Base"] + data["nextBase"],
                                                             np.where(data["Pos"] % 3 == 1,
                                                                      data["prev2bases"] + data["Base"],
                                                                      data["Context"]))))
    data.dropna(subset=["Consensus_codon", "Mutated_codon"], inplace=True)
    data["Consensus_aa"] = data["Consensus_codon"].apply(lambda x: Seq(x).translate()[0] if "-" not in x else "-")
    data["Mutated_aa"] = data["Mutated_codon"].apply(lambda x: Seq(x).translate()[0] if "-" not in x else "-")
    data["Type"] = np.where(data["Mutated_aa"] == data["Consensus_aa"], "Synonymous", np.where(data["Mutated_aa"] == "*"
                                                                                               , "Premature Stop Codon",
                                                                                               "Non-Synonymous"))
    return data


def freqs_to_dataframe(freqs_file):
    """
    This function returns arranged DataFrame without deletions
    :param freqs_file: The path of the relevant freqs file
    :return: DataFrame without deletions
    """
    data = pd.read_table(freqs_file)
    data = data[data.Ref != '-']
    data = data[data.Base != '-']
    data.reset_index(drop=True, inplace=True)
    return data


def find_coding_region(ncbi_id):
    """
    :param ncbi_id: NCBI_ID number
    :return:the start and the end positions of the Coding region
    """
    try:
        Entrez.email = "A.N.Other@example.com"
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=ncbi_id)
        ncbi_gb = SeqIO.read(handle, "gb")
        handle.close()
        start_pos = 0
        end_pos = 0
        if ncbi_id == "M33854":
            start_pos = 739
            end_pos = 7299
            return start_pos, end_pos
        for feature in ncbi_gb.features:
            if feature.type == 'CDS':
                if end_pos < feature.location.end:
                    start_pos = feature.location.start + 1
                    end_pos = feature.location.end + 0
                return start_pos, end_pos
    except Exception as e:
        print(e)


def filter_by_coverage_mutation(type_file, output_dir, min_read_count=None):
    """
    Creates DataFrame of frequencies with Mutations X->Y column
    :param type_file: freqs file after it has ben processed with find_mutation_type func
    :param output_dir: where to save the returned DF
    :param min_read_count: Default 10,0000
    :return: DF with Mutation and mutation_type column with minimal coverage
    """
    data = type_file
    data.reset_index(drop=True, inplace=True)
    flag = '-' in data.Base.values
    if flag is True:
        data = data[data.Ref != '-']
        data = data[data.Base != '-']
        data.reset_index(drop=True, inplace=True)
    data['Base'].replace('T', 'U', inplace=True)
    data['Ref'].replace('T', 'U', inplace=True)
    if min_read_count is None:
        min_read_count = 100000
    data = data[data['Pos'] == np.round(data['Pos'])]  # remove insertions
    data['Pos'] = data[['Pos']].apply(pd.to_numeric)
    data = data[data['Read_count'] > min_read_count]
    data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
    data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
    data["Mutation"] = data["Ref"] + ">" + data["Base"]
    data = data[["Pos", "Base", "Freq", "Ref", "Read_count", "Rank", "Prob",
                                    "counts_for_position", "Consensus", "BaseToBase", "Mutation_class", "adj_freq",
                                    "prevBase", "nextBase", "prev2bases", "next2bases", "Consensus_codon",
                                    "Mutated_codon", "Consensus_aa", "Mutated_aa", "Type", "abs_counts", "Frequency",
                                    "Mutation", "label", "passage", "replica"]]
    data.to_csv(output_dir + "/data_mutation.csv", sep=',', encoding='utf-8')
    return data



def main():
    sample1 = "ERR2352296"
    ncbi_id1 = "MH732737"

    sample2 = "SRR6955546"
    ncbi_id2 = "V01149"

    sample3 = "SRR8196193"
    ncbi_id3 = "NC_001430"

    sample4 = "RV-p11"
    ncbi_id4 = "NC_001490"
    freqs_file = "/Users/odedkushnir/Projects/fitness/test/%s/%s.freqs" % (sample4, sample4)
    data = find_mutation_type(freqs_file, ncbi_id4)

    data.to_csv("/Users/odedkushnir/Projects/fitness/test/%s/data_cirseq.csv" % (sample4), sep=',', encoding='utf-8' )
    print(data)

if __name__ == "__main__":
    main()