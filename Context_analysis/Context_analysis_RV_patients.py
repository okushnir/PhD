#! /usr/local/python/anaconda_python-3.6.1

"""
@Author: odedkushnir

"""


import os.path
from Utilities import sequnce_utilities
import glob
import pandas as pd
import numpy as np
import urllib



def checkKey(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        raise Exception()

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


    #for Local

    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/patients"
    output_dir = input_dir + "/Rank0_data_mutation"
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
    min_coverage = 5000
    dirs = glob.glob(input_dir + "/*")
    lst_srr = []
    for passage in dirs:
        # Checks if the file .merged.with.mutation.type.freqs file exists
        file_path = glob.glob(passage + "/20201017_q30_consensusX5/*.merged*freqs") #
        if len(file_path) >= 1:
            if "merged.with.mutation.type.freqs" in str(file_path):
                for file in file_path:
                    if (file.split(".")[-2] == "type") & (len(file.split(".")) > 3):
                        lst_srr.append(file)
                        ncbi_id = "NC_001490"
            else:
                print(file_path[0].split('/')[-1].split(".")[0].split("-")[0])
                # virus = os.path.basename(file_path[0].split('/')[-1].split(".")[0].split("-")[0])
                virus = "RVA"
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

    append_mutation2 = sequnce_utilities.find_mutation_type(
        "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/controls/IVT_5_"
        "Control/20201012_q38/IVT-5-Control.merged.freqs", "NC_001490", min_coverage)
    append_mutation3 = sequnce_utilities.find_mutation_type(
        "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/controls/p3_Control/"
        "20201012_q38/p3-Control.merged.freqs", "NC_001490", min_coverage)
# creating data_mutation.csv file

    sample_file0 = lst_srr[0]
    sample_file1 = lst_srr[1]
    sample_file2 = lst_srr[2]
    sample_file3 = lst_srr[3]
    sample_file4 = lst_srr[4]
    sample_file5 = lst_srr[5]
    sample_file6 = lst_srr[6]

    #Controls

    # control_file_rnd = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/RVB14_RNA-Control_L001-ds.850f2223182a41fcaec41c4f3735e428/q38_3UTR/RVB14-RNA-Control.merged.with.mutation.type.freqs"
    # label_control1 = "RNA Control_RND"

    control_file_id = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/controls/IVT_5_Control/20201012_q38/IVT-5-Control.merged.with.mutation.type.freqs"
    label_control2 = "RNA Control\nPrimer ID"

    control_file_cell = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/controls/p3_Control/20201012_q38/p3-Control.merged.with.mutation.type.freqs"
    label_control3 = "Cell Cultureֿ\nControl"

    # next_file = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/RVB14_Next_RNA_Control/q35/RVB14-Next_Control.merged.with.mutation.type.freqs"
    # label_next = "RVB14-Next-RNA Control"

    label_sample0 = sample_file0.split("/")[-1].split(".")[0]

    label_sample1 = sample_file1.split("/")[-1].split(".")[0]

    label_sample2 = sample_file2.split("/")[-1].split(".")[0]

    label_sample3 = sample_file3.split("/")[-1].split(".")[0]

    label_sample4 = sample_file4.split("/")[-1].split(".")[0]

    label_sample5 = sample_file5.split("/")[-1].split(".")[0]

    label_sample6 = sample_file6.split("/")[-1].split(".")[0]

    print ("loading " + sample_file0 + " as sample")
    data_mutations0 = pd.read_table(sample_file0)
    data_mutations0["label"] = label_sample0


    print ("loading " + sample_file1 + " as sample")
    data_mutations1 = pd.read_table(sample_file1)
    data_mutations1["label"] = label_sample1


    print("loading " + sample_file2 + " as sample")
    data_mutations2 = pd.read_table(sample_file2)
    data_mutations2["label"] = label_sample2

    print("loading " + sample_file3 + " as sample")
    data_mutations3 = pd.read_table(sample_file3)
    data_mutations3["label"] = label_sample3

    print("loading " + sample_file4 + " as sample")
    data_mutations4 = pd.read_table(sample_file4)
    data_mutations4["label"] = label_sample4

    print("loading " + sample_file5 + " as sample")
    data_mutations5 = pd.read_table(sample_file5)
    data_mutations5["label"] = label_sample5

    print("loading " + sample_file6 + " as sample")
    data_mutations6 = pd.read_table(sample_file6)
    data_mutations6["label"] = label_sample6

    print("loading " + control_file_id + " as RNA control")
    data_control2 = pd.read_table(control_file_id)
    data_control2["label"] = label_control2

    print("loading " + control_file_cell + " as RNA control")
    data_control3 = pd.read_table(control_file_cell)
    data_control3["label"] = label_control3

    data = pd.concat([data_mutations0, data_mutations1, data_mutations2, data_mutations3, data_mutations4,
                      data_mutations5, data_mutations6, data_control2, data_control3], sort=False)

    toPlot = [label_sample0, label_sample1, label_sample2, label_sample3, label_sample4, label_sample5, label_sample6,
              label_control2, label_control3]
    # filter_by_coverage_mutation(data, input_dir, min_read_count=min_coverage)



# Context data

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
    consensus_data["major_read_count"] = np.round(consensus_data['Read_count'] * consensus_data['Freq'])
    consensus_data = consensus_data.rename(columns={"Base": "Consensus"})
    data = pd.merge(data, consensus_data[['Pos', 'Prev', 'Next', 'Context', 'Consensus', "label", "major_read_count"]],
                    on=["Pos", "label"])
    data["mutation_type"] = data["Consensus_y"] + data["Base"]
    data["Mutation"] = data["Consensus_y"] + ">" + data["Base"]
    """Major Variant"""
    data = data[(data["label"].isin(toPlot)) & (data["Read_count"] > min_coverage)]
    """Minor variant"""
    # data = data[(data["label"].isin(toPlot)) & (data["Rank"] != 0) & (data["Read_count"] > min_coverage)] #

    data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
    data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
    # start_pos, end_pos = sequnce_utilities.find_coding_region(ncbi_id)
    # No internet connection - for RVA
    start_pos, end_pos = 503, 6728
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
    data["Protein"] = np.where(data["Pos"] <= 503, "5'UTR",
                                np.where(data["Pos"] <= 709, "1A",
                                    np.where(data["Pos"] <= 1495, "1B",
                                        np.where(data["Pos"] <= 2197, "1C",
                                            np.where(data["Pos"] <= 3094, "1D",
                                                np.where(data["Pos"] <= 3523, "2A",
                                                    np.where(data["Pos"] <= 3808, "2B",
                                                        np.where(data["Pos"] <= 4771, "2C",
                                                            np.where(data["Pos"] <= 5005, "3A",
                                                                np.where(data["Pos"] <= 5068, "3B",
                                                                   np.where(data["Pos"] <= 5617, "3C",
                                                                    np.where(data["Pos"] <= 6728, "3D", "3'UTR"))))))))))))
    data["Type"] = data["Type"].fillna(value="NonCodingRegion")
    data.to_csv(output_dir + "/q30_data_mutation.csv", sep=',', encoding='utf-8')
    data.to_pickle(output_dir + "/q30_data_mutation.pkl")

    mutation_for_rna = ["AG"]
    dataForPlotting_AG = data[(data["mutation_type"].isin(mutation_for_rna))]
    dataForPlotting_AG.to_csv(input_dir + "/Rank0_data_mutation/q30_data_XpA_by_mutation.csv", sep=',', encoding='utf-8')
    dataForPlotting_AG.to_pickle(input_dir + "/Rank0_data_mutation/q30_data_XpA_by_mutation.pkl")

    mutation_for_rna = ["UC"]
    dataForPlotting_UC = data[(data["mutation_type"].isin(mutation_for_rna))]
    dataForPlotting_UC.to_csv(input_dir + "/q30_data_UpX_by_mutation.csv", sep=',', encoding='utf-8')
    dataForPlotting_UC.to_pickle(input_dir + "/q30_data_UpX_by_mutation.pkl")

if __name__ == "__main__":
    main()