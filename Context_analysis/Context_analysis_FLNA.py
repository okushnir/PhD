#! /usr/local/python/anaconda_python-3.6.1

"""
@Author: odedkushnir

"""


import os.path
from sequnce_utilities import *
import glob
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

    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/FLNA"

    min_coverage = 10000
    dirs = glob.glob(input_dir + "/FLNA_RV*")
    lst_srr=[]
    for passage in dirs:
        file_path = glob.glob(passage + "/q30_3UTR/*" + ".freqs")
        if len(file_path) > 1:
            for file in file_path:
                if file.split(".")[-2] == "type":
                    lst_srr.append(file)
        else:
            if file_path == []:
                break
            else:

                org_dic = {"CVB3": "M16572", "RVB14": "NC_001490", "RV": "NC_001490", "Echovirus E7": "MH732737",
                           "Coxsackievirus A16": "NC_001612", "Enterovirus A": "NC_001612", "Echovirus E3": "KX808644",
                           "Coxsackievirus B2": "AF081485", "Echovirus E25": "KJ957190", "Human poliovirus 1": "V01149",
                           "Human poliovirus 2": "V01149", "Human poliovirus 3": "V01149", "Enterovirus C": "V01149",
                           "Enterovirus D68": "NC_001430", "Enterovirus D": "NC_001430", "Rhinovirus B": "NC_001490",
                           "Coxsackievirus B3 (strain Nancy)": "JN048468", "Rhinovirus C": "LC428177",
                           "Echovirus E6": "JX976771", "FLNA": "NM_001456", "ATCG1": "NM_001614"}
                print(file_path[0].split('/')[-1].split(".")[0].split("-")[0])
                virus = os.path.basename(file_path[0].split('/')[-1].split(".")[0].split("-")[0])
                try:
                    ncbi_id = checkKey(org_dic, virus)
                    print("Adding Mutation type to:%s" % (file_path[0].split('/')[-1].split(".")[0]))
                except Exception as e:
                    print("type error: " + str(e))
                try:
                    append_mutation = find_mutation_type(file_path[0], ncbi_id, min_read_count=10000)
                    lst_srr.append(file_path[0].split('.')[0] + ".with.mutation.type.freqs")
                except Exception as e:
                    print(e)
                continue
    print(lst_srr)


    # creating data_mutation.csv file

    sample_file0 = lst_srr[0]
    # sample_file1 = lst_srr[1]



    label_sample0 = sample_file0.split("/")[-1].split(".")[0]

    # label_sample1 = sample_file1.split("/")[-1].split(".")[0]



    print("loading " + sample_file0 + " as sample")
    data_mutations0 = pd.read_table(sample_file0)
    data_mutations0["label"] = label_sample0

    # print("loading " + sample_file1 + " as sample")
    # data_mutations1 = pd.read_table(sample_file1)
    # data_mutations1["label"] = label_sample1



    data = pd.concat([data_mutations0], sort=False)

    data["passage"] = 0
    data["replica"] = 0

    # data["passage"] = np.where(data["label"] == "RV-p71", 7, (np.where(data["label"] == "RV-p11", 1, (np.where(data["label"] == "RV-p12",
    #                                                            1, 0)))))
    # data["replica"] = np.where(data["label"] == "RV-p71", 1, (np.where(data["label"] == "RV-p11", 1, (np.where(data["label"] == "RV-p12",
    #                                                            2, 1)))))
    print(data.to_string())

    toPlot = [label_sample0]

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
    data = data[(data["label"].isin(toPlot)) & (data["Rank"] != 0) & (data["Read_count"] > min_coverage)]
    data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
    data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
    data.to_csv(input_dir + "/data_mutation.csv", sep=',', encoding='utf-8')

    mutation_for_rna = ["AG"]
    dataForPlotting_AG = data[(data["mutation_type"].isin(mutation_for_rna))]#& (data["Type"] == "Synonymous")]
    dataForPlotting_AG.to_csv(input_dir + "/data_XpA_by_mutation.csv", sep=',', encoding='utf-8')

    mutation_for_rna = ["UC"]
    dataForPlotting_UC = data[(data["mutation_type"].isin(mutation_for_rna))]  # & (data["Type"] == "Synonymous")]
    dataForPlotting_UC.to_csv(input_dir + "/data_UpX_by_mutation.csv", sep=',', encoding='utf-8')

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
#
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

    # rv_df = pd.DataFrame(columns=("Pos", "Base", "Freq", "Ref", "Read_count", "Rank", "Prob"))
    # for passage in dirs:
    #     file_path = glob.glob(passage + "/q30_3UTR_new/*" + ".freqs")
    #     passage_df = freqs_to_dataframe(file_path[0])
    #     passage_df["label"] = file_path[0].split("/")[-1].split(".")[0]
    #     rv_df = rv_df.append(passage_df, ignore_index=True, sort=False)
    #
    #
    # rv_df = add_mutation_type_to_df(rv_df, "NC_001490", 10000)
    #
    # print(rv_df)