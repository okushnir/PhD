#! /usr/local/python/anaconda_python-3.6.1

"""
@Author: odedkushnir

"""


import os.path
from Utilities.sequnce_utilities import *
import glob
# import SRR_analysis_from_Cluster


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
    # SRR_dir = "/Users/odedkushnir/Projects/fitness/SRP/PV"
    # SRR_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/SRA/SRP006391_RV" #9
    # SRR_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/SRA1/SRP064468_PV" #10
    # SRR_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/SRA2/SRP137824_Entero_C" #10
    # SRR_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/SRA3/SRP155896_Entro_D" #6
    # SRR_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/SRA4/ERP014415_Entero_A" #10
    SRR_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/ADAR_Influenza/SRP091383" #6

    # SRR_dir = "/sternadi/home/volume3/okushnir/SRA4//ERP014415_Entero_A"
    org_df = pd.read_csv(SRR_dir + "/Top10SRR1.csv", sep=',', encoding='utf-8')
    org_df = org_df[["Organism", "Run"]]
    org_df = org_df.rename(columns={"Run": "label"})
    # print(org_df)
    min_coverage = 1
    #
    #for ERP014415_Entero_A
    # dirs = glob.glob(SRR_dir + "/RR*")

    dirs = glob.glob(SRR_dir + "/SRR*")
    lst_srr = []
    for SRR in dirs:
        file_path = glob.glob(SRR + "/q30_consensus_1e-03/*" + ".freqs")
        if len(file_path) > 1:
            for file in file_path:
                if file.split(".")[-2] == "type":
                    lst_srr.append(file)
        else:
            if file_path == []:
                break
            else:
                srr_dic = {"SRR2558452": "KJ170553", "SRR2558453": "KJ170542", "SRR2558454": "KJ170565",
                           "SRR2558455": "KJ170571", "SRR2558456": "KJ170550", "SRR2558457": "KJ170560",
                           "SRR2558518": "KJ170670", "SRR2558529": "KJ170579", "SRR2558449": "KJ170541",
                           "SRR2558492": "KJ170587", "ERR2352296": "MH732737", "ERR2352295": "JX976771",
                           "ERR2352351": "NC_001612", "ERR2352354": "NC_001612", "ERR2352352": "NC_001612",
                           "ERR2352292": "KX808644", "ERR2352353": "NC_001612", "ERR2352267": "AF081485",
                           "ERR2352317": "MH732737", "ERR2352290": "KJ957190", "SRR352020": "LC428177",
                           "SRR525047": "NC_001490", "SRR524812": "LC428177", "SRR363436": "LC428177",
                           "SRR351937": "LC428177", "SRR350571": "LC428177", "SRR359927": "LC428177",
                           "SRR359925": "LC428177", "SRR389042": "LC428177", "SRR350572": "LC428177",
                           "SRR6955546": "V01149", "SRR6955559": "V01149", "SRR6955555": "V01149", "SRR6955554": "V01149",
                           "SRR6955547": "V01149", "SRR6955549": "V01149", "SRR6955548": "V01149", "SRR6955539": "V01149",
                           "SRR6955537": "V01149", "SRR6955563": "V01149", "SRR7630054": "NC_001430",
                           "SRR8196193": "NC_001430", "SRR8196586": "NC_001490", "SRR7630002": "NC_001430",
                           "SRR7630178": "NC_001430", "SRR8196335": "NC_001430", "SRR7630107": "NC_001430",
                           "SRR8196251": "NC_001430", "SRR8196446": "NC_001430", "SRR4414064": "CY164045",
                           "SRR4414065": "CY164045", "SRR4414066": "CY164045", "SRR4414067": "CY164045",
                           "SRR4414068": "CY164045", "SRR4414069": "CY164045"}

                org_dic = {"CVB3": "M16572", "RVB14": "NC_001490", "Echovirus E7": "MH732737","Echovirus E6": "JX976771",
                           "Coxsackievirus A16": "NC_001612", "Enterovirus A": "NC_001612", "Echovirus E3": "KX808644",
                           "Coxsackievirus B2": "AF081485", "Echovirus E25": "KJ957190", "Human poliovirus 1": "V01149",
                           "Human poliovirus 2": "V01149", "Human poliovirus 3": "V01149", "Enterovirus C": "V01149",
                           "Enterovirus D68": "NC_001430", "Enterovirus D": "NC_001430", "Rhinovirus B": "NC_001490",
                           "Coxsackievirus B3 (strain Nancy)": "JN048468", "Rhinovirus C": "LC428177"}
                srr = os.path.basename(file_path[0]).split('.')[0]
                try:
                    ncbi_id = checkKey(srr_dic, srr)
                    print("Adding Mutation type to:%s" % (srr))
                except Exception as e:
                    print("type error: " + str(e))
                    print("The SRR is not the SRR dictionary")
                    ncbi_id = checkKey(org_dic, srr)
                    print("Adding Mutation type to:%s" % (srr))
                try:
                    append_mutation = find_mutation_type(file_path[0], ncbi_id)
                    lst_srr.append(file_path[0].split('.')[0] + ".with.mutation.type.freqs")
                except Exception as e:
                    print(e)
                continue
    print(lst_srr)


    # # creating data_mutation.csv file

    sample_file0 = lst_srr[0]
    sample_file1 = lst_srr[1]
    sample_file2 = lst_srr[2]
    sample_file3 = lst_srr[3]
    sample_file4 = lst_srr[4]
    sample_file5 = lst_srr[5]
    # sample_file6 = lst_srr[6]
    # sample_file7 = lst_srr[7]
    # sample_file8 = lst_srr[8]
    # sample_file9 = lst_srr[9]


    # # control_file = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/RV-IVT/q30_3UTR_new/RV-IVT.with.mutation.type.freqs"
    # #
    # # label_control = "RNA Control"
    #
    label_sample0 = sample_file0.split("/")[-1].split(".")[0]

    label_sample1 = sample_file1.split("/")[-1].split(".")[0]

    label_sample2 = sample_file2.split("/")[-1].split(".")[0]

    label_sample3 = sample_file3.split("/")[-1].split(".")[0]

    label_sample4 = sample_file4.split("/")[-1].split(".")[0]

    label_sample5 = sample_file5.split("/")[-1].split(".")[0]

    # label_sample6 = sample_file6.split("/")[-1].split(".")[0]
    #
    # label_sample7 = sample_file7.split("/")[-1].split(".")[0]
    #
    # label_sample8 = sample_file8.split("/")[-1].split(".")[0]
    #
    # label_sample9 = sample_file9.split("/")[-1].split(".")[0]


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

    # print("loading " + sample_file6 + " as sample")
    # data_mutations6 = pd.read_table(sample_file6)
    # data_mutations6["label"] = label_sample6
    #
    # print("loading " + sample_file7 + " as sample")
    # data_mutations7 = pd.read_table(sample_file7)
    # data_mutations7["label"] = label_sample7
    #
    # print("loading " + sample_file8 + " as sample")
    # data_mutations8 = pd.read_table(sample_file8)
    # data_mutations8["label"] = label_sample8
    #
    # print("loading " + sample_file9 + " as sample")
    # data_mutations9 = pd.read_table(sample_file9)
    # data_mutations9["label"] = label_sample9


    # print("loading " + control_file + " as homogeneous control")
    # data_control = pd.read_table(control_file)
    # data_control["source"] = label_control

    data = pd.concat([data_mutations0, data_mutations1, data_mutations2, data_mutations3, data_mutations4,
                      data_mutations5], sort=False) #, data_mutations6, data_mutations7, data_mutations8, data_mutations9
    data["passage"] = 0
    data["replica"] = 0

    data = pd.merge(left=data, right=org_df, on="label")

    toPlot = [label_sample0, label_sample1, label_sample2, label_sample3, label_sample4, label_sample5] #, label_sample6,
              # label_sample7, label_sample8, label_sample9

    filter_by_coverage_mutation(data, SRR_dir, min_read_count=min_coverage)

#Context data

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


    mutation_for_rna = ["AG"]
    dataForPlotting = data[(data["label"].isin(toPlot)) & (data["mutation_type"].isin(mutation_for_rna)) &
                           (data["Rank"] != 0) & (data["Read_count"] > min_coverage)]#& (data["Type"] == "Synonymous")]
    dataForPlotting.to_csv(SRR_dir + "/data_XpA_by_mutation.csv", sep=',', encoding='utf-8')

    mutation_for_rna = ["UC"]
    dataForPlotting_UC = data[(data["label"].isin(toPlot)) & (data["mutation_type"].isin(mutation_for_rna)) &
                           (data["Rank"] != 0) & (
                                       data["Read_count"] > min_coverage)]  # & (data["Type"] == "Synonymous")]
    dataForPlotting_UC.to_csv(SRR_dir + "/data_UpX_by_mutation.csv", sep=',', encoding='utf-8')

#     fig1 = plt.figure(figsize=(16, 9))
#     ax = plt.subplot()
#     make_boxplot_mutation(freqs_file_mutations, ax, out_plots_dir)
#     plt.savefig(out_plots_dir + suffix.split(sep='.')[0] + '_All_Mutations.png', dpi=300)
#     plt.close("all")
#     print("The All Mutations Plot is ready in the folder")
#
#     fig2 = plt.figure(figsize=(16, 9))
#     ax = plt.subplot()
#     make_boxplot_transition_mutation(freqs_file_mutations, ax, out_dir)
#     plt.savefig(out_plots_dir + suffix.split(sep='.')[0] + '_Transitions_Report.png', dpi=300)
#     plt.close("all")
#     print("The Transition Plot is ready in the folder")
#
#
# """Graphs"""
#     plots_dir = "/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/"
#
#
# #6. Mutation Rates
# def make_boxplot_mutation(freqs_file, ax, output_dir):
#     """
#     Plots the mutation frequencies boxplot
#     :param freqs_file: pandas DataFrame after find_mutation_type function
#     :param ax: ax location
#     :return:
#     """
#     data = pd.read_table(freqs_file)
#     data.reset_index(drop=True, inplace=True)
#     flag = '-' in data.Base.values
#     if flag is True:
#         data = data[data.Ref != '-']
#         data = data[data.Base != '-']
#         data.reset_index(drop=True, inplace=True)
#         # data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
#         # data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
#         # raise Exception("This script does not support freqs file with deletions, for support please contact Maoz ;)"
#     data['Base'].replace('T', 'U', inplace=True)
#     data['Ref'].replace('T', 'U', inplace=True)
#     min_read_count = 1000
#     data = data[data['Pos'] == np.round(data['Pos'])]  # remove insertions
#     data['Pos'] = data[['Pos']].apply(pd.to_numeric)
#     data = data[data['Read_count'] > min_read_count]
#     data['mutation_type'] = data['Ref'] + data['Base']
#     data = data[data['Ref'] != data['Base']]
#     data = data[data["Base"] != "-"]
#     data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
#     data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
#     data["Mutation"] = data["Ref"] + "->" + data["Base"]
#     data.to_csv(output_dir + "data_all_mutation.csv", sep=',', encoding='utf-8')
#     sns.set_palette(sns.color_palette("Paired", 12))
#     g1 = sns.boxplot(x="Mutation Type", y="Frequency", hue="Mutation", data=data,
#                      hue_order=["C->U", "U->C", "G->A", "A->G", "C->A", "G->U", "U->G", "U->A", "G->C", "A->C",
#                                 "A->U", "C->G"], order=["Synonymous", "Non-Synonymous", "Premature Stop Codon"], ax=ax)
#     g1.set(yscale="log")
#     plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0., fontsize=6)
#     g1.set_ylim(10 ** -6, 1)
#     g1.tick_params(labelsize=7)
#
# def make_boxplot_transition_mutation(freqs_file,ax, output_dir):
#     """
#     Plots the mutation frequencies boxplot
#     :param freqs_file: pandas DataFrame after find_mutation_type function
#     :param ax: ax location
#     :return:
#     """
#     data = pd.read_table(freqs_file)
#     data.reset_index(drop=True, inplace=True)
#     flag = '-' in data.Base.values
#     if flag is True:
#         data = data[data.Ref != '-']
#         data = data[data.Base != '-']
#         data.reset_index(drop=True, inplace=True)
#         # data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
#         # data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
#         # raise Exception("This script does not support freqs file with deletions, for support please contact Maoz ;)"
#     data['Base'].replace('T', 'U', inplace=True)
#     data['Ref'].replace('T', 'U', inplace=True)
#     min_read_count = 1000
#     data = data[data['Pos'] == np.round(data['Pos'])]  # remove insertions
#     data['Pos'] = data[['Pos']].apply(pd.to_numeric)
#     data = data[data['Read_count'] > min_read_count]
#     data['mutation_type'] = data['Ref'] + data['Base']
#     data = data[data['Ref'] != data['Base']]
#     data = data[data["Base"] != "-"]
#     data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
#     data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
#     data["Mutation"] = data["Ref"] + "->" + data["Base"]
#
#     data.to_csv(output_dir + "data_transitions_mutation.csv", sep=',', encoding='utf-8')
#
#     sns.set_palette(sns.color_palette("Paired", 12))
#
#     g1 = sns.factorplot(x="Mutation Type", y="Frequency", data=data, col="Mutation",
#                      col_order=["C->U", "U->C", "G->A", "A->G"], kind="box")
#     g1.set_xticklabels(["Synonymous", "Non\nSynonymous", "Premature\nStop Codon"], fontsize=10)
#     g1.set_xlabels('')
#     g1.set(yscale="log", ylim=(0.000001, 0.01))


if __name__ == "__main__":
    main()
