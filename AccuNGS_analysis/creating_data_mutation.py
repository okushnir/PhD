#! /usr/local/python/anaconda_python-3.6.1

"""
@Author: odedkushnir

"""


import os.path
from Utilities import sequnce_utilities
import glob
import pandas as pd
import numpy as np
from AccuNGS_analysis.add_Protein_to_pd_df import add_Protein_to_pd_df_func
import urllib



def checkKey(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        raise Exception()


def creting_data_mutation_df(input_dir, min_coverage, virus, date, q, control_dict={}):

    # input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14"

    output_dir = input_dir + "/test"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    org_dic = {"CVB3": "M33854", "RVB14": "NC_001490", "RV": "NC_001490", "Echovirus E7": "MH732737",
               "Coxsackievirus A16": "NC_001612", "Enterovirus A": "NC_001612", "Echovirus E3": "KX808644",
               "Coxsackievirus B2": "AF081485", "Echovirus E25": "KJ957190", "Human poliovirus 1": "V01149",
               "Human poliovirus 2": "V01149", "Human poliovirus 3": "V01149", "Enterovirus C": "V01149",
               "Enterovirus D68": "NC_001430", "Enterovirus D": "NC_001430", "Rhinovirus B": "NC_001490",
               "Coxsackievirus B3 (strain Nancy)": "JN048468", "Rhinovirus C": "LC428177",
               "Echovirus E6": "JX976771", "RVA": "JX025555"}

    dirs = glob.glob(input_dir + "/p*")
    lst_srr = []
    for passage in dirs:
        # Checks if the file .merged.with.mutation.type.freqs file exists
        file_path = glob.glob(passage + "/%s_%s/*.merged*freqs" % (date, q))
        if len(file_path) >= 1:
            if "merged.with.mutation.type.freqs" in str(file_path):
                for file in file_path:
                    if (file.split(".")[-2] == "type") & (len(file.split(".")) > 3):
                        lst_srr.append(file)
                        ncbi_id = checkKey(org_dic, virus)
            else:
                print(file_path[0].split('/')[-1].split(".")[0].split("-")[0])
                # virus = os.path.basename(file_path[0].split('/')[-1].split(".")[0].split("-")[0])
                try:
                    ncbi_id = checkKey(org_dic, virus)
                    print("Adding Mutation type to:%s" % (file_path[0].split('/')[-1].split(".")[0]))
                except Exception as e:
                    print("type error: " + str(e))
                try:
                    append_mutation = sequnce_utilities.find_mutation_type(file_path[0], ncbi_id, min_coverage)
                    lst_srr.append(file_path[0].split('freqs')[0] + "with.mutation.type.freqs")
                except Exception as e:
                    print(e)
                continue
    print(lst_srr)

    columns = ["Pos", "Base", "Freq", "Ref", "Read_count", "Rank", "Prob", "counts_for_position", "Consensus",
               "BaseToBase", "Mutation_class", "adj_freq", "prevBase", "nextBase", "prev2bases", "next2bases",
               "Consensus_codon", "Mutated_codon", "Consensus_aa", "Mutated_aa", "Type", "pval", "Var_perc",
               "SNP_Profile"]
    data = pd.DataFrame(columns=columns)
    for i in range(len(lst_srr)):
        sample_file = lst_srr[i]
        label_sample = sample_file.split("/")[-1].split(".")[0]
        print("loading " + sample_file + " as sample")
        data_mutations = pd.read_table(sample_file)
        data_mutations["label"] = label_sample
        data_mutations["passage"] = label_sample.split("-")[0].split("p")[-1]
        data_mutations["replica"] = label_sample.split("-")[-1]
        data = data.append(data_mutations, ignore_index=True)

    for key, value in control_dict.items():
        print("loading " + value + " as RNA control RND")
        data_control = pd.read_table(value)
        data_control["label"] = key
        data_control["passage"] = 0
        if key == "RNA Control_RND":
            data_control["replica"] = 3
        else:
            data_control["replica"] = 1
        data = data.append(data_control, ignore_index=True)
        # data_control["passage"] = 0
        # data_control["replica"] = 3


    # data = pd.concat([data_mutations_lst, data_controls_lst], sort=False)

    data['Base'].replace('T', 'U', inplace=True)
    data['Ref'].replace('T', 'U', inplace=True)
    # generate consensus of both sample and control
    consensus_data = data[['Pos', 'Base', 'label', "Read_count", "Freq"]][data['Rank'] == 0]
    consensus_data = consensus_data[consensus_data['Pos'] == np.round(consensus_data['Pos'])]  # removes insertions
    consensus_data['Next'] = consensus_data['Base'] + consensus_data.shift(-1)['Base']
    consensus_data['Prev'] = consensus_data.shift(1)['Base'] + consensus_data['Base']
    consensus_data['Pos'] = consensus_data[['Pos']].apply(pd.to_numeric)
    consensus_data["Context"] = consensus_data.shift(1)['Base'] + consensus_data['Base'] + consensus_data.shift(-1)[
        'Base']
    consensus_data["major_read_count"] = consensus_data['Read_count'] * consensus_data['Freq']
    consensus_data = consensus_data.rename(columns={"Base": "Consensus"})
    data = pd.merge(data, consensus_data[['Pos', 'Prev', 'Next', 'Context', 'Consensus', "label", "major_read_count"]],
                    on=["Pos", "label"])
    data["mutation_type"] = data["Consensus_y"] + data["Base"]
    data["Mutation"] = data["Consensus_y"] + ">" + data["Base"]
    data = data[(data["Read_count"] > min_coverage)]
    data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
    data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
    start_pos, end_pos = sequnce_utilities.find_coding_region(ncbi_id)

    # No internet connection - for RV
    # start_pos, end_pos = 629, 7212

    start_pos_mod3 = start_pos % 3

    if (start_pos_mod3 == 0) | (start_pos_mod3 == 1):
        data["Codon_Pos"] = np.where(data["Pos"] % 3 == start_pos_mod3, 1,
                                           np.where(data["Pos"] % 3 == start_pos_mod3 + 1, 2,
                                                   3))
    elif start_pos_mod3 == 2:
        data["Codon_Pos"] = np.where(data["Pos"] % 3 == start_pos_mod3, 1,
                                           np.where(data["Pos"] % 3 == start_pos_mod3 + 1,
                                                   2,
                                                    np.where(data["Pos"] % 3 == 0, 2,
                                                             np.where(data["Pos"] % 3 == 1,
                                                                      3,
                                                                      0))))
    data["Consensus>Mutated_codon"] = data["Consensus_codon"] + ">" + data["Mutated_codon"]
    data["Type"] = data["Type"].fillna(value="NonCodingRegion")
    region_lst = [629, 835, 1621, 2329, 3196, 3634, 3925, 4915, 5170, 5239, 5785, 7165]
    data = add_Protein_to_pd_df_func(data, region_lst)
    data.to_csv(output_dir + "/q38_data_mutation.csv", sep=',', encoding='utf-8')
    data.to_pickle(output_dir + "/q38_data_mutation.pkl")


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
    min_coverage = 5000
    virus = "RVB14"
    date = "20201012"
    q = "q38"
    control_file_rnd = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/RVB14_RNA-Control/q38_3UTR/" \
                       "RVB14-RNA-Control.merged.with.mutation.type.freqs"
    label_control1 = "RNA Control_RND"
    control_file_spe = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/controls/" \
                       "IVT_3_Control/20201012_q38/IVT-3-Control.merged.with.mutation.type.freqs"
    label_control2 = "RNA Control\nPrimer ID"

    control_dict = {label_control1: control_file_rnd, label_control2: control_file_spe}

    creting_data_mutation_df(input_dir, min_coverage, virus, date, q, control_dict)

    """CV"""

    """RV-Capsid_Free"""

    """RV-Patients"""

if __name__ == "__main__":
    main()
