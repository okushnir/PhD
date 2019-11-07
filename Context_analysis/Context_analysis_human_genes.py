#! /usr/local/python/anaconda_python-3.6.1

"""
@Author: odedkushnir

"""


import os.path
from sequnce_utilities import *
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



    gene_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/Cirseq/RV/Human/"
    lst_gene = []
    min_coverage = 1000
    dirs = glob.glob(gene_dir + "/20190529_q23_*")

    for dir in dirs:
        file_path = glob.glob(dir + "/*.freqs")
        if len(file_path) > 1:
            for file in file_path:
                if file.split(".")[-2] == "type":
                    lst_gene.append(file)
        else:
            if file_path == []:
                break
            else:
                gene_dic = {"IFNG": "NM_000619", "AZIN1": "NM_148174", "FLNA": "NM_001456"}
                gene = os.path.basename(file_path[0]).split('.')[0]
                try:
                    ncbi_id = checkKey(gene_dic, gene)
                    print("Adding Mutation type to:%s" % (gene))
                except Exception as e:
                    print("type error: " + str(e))
                    print("The gene is not the gene dictionary")
                    # ncbi_id = checkKey(gene_dic, gene)
                    # print("Adding Mutation type to:%s" % (srr))
                try:
                    append_mutation = find_mutation_type(file_path[0], ncbi_id, min_coverage)
                    lst_gene.append(file_path[0].split('.')[0] + ".with.mutation.type.freqs")
                except Exception as e:
                    print(e)
                continue
    print(lst_gene)


    # # creating data_mutation.csv file

    sample_file0 = lst_gene[0]
    sample_file1 = lst_gene[1]
    # sample_file2 = lst_gene[2]
    # sample_file3 = lst_srr[3]
    # sample_file4 = lst_srr[4]
    # sample_file5 = lst_srr[5]
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
    #
    # label_sample2 = sample_file2.split("/")[-1].split(".")[0]

    # label_sample3 = sample_file3.split("/")[-1].split(".")[0]
    #
    # label_sample4 = sample_file4.split("/")[-1].split(".")[0]
    #
    # label_sample5 = sample_file5.split("/")[-1].split(".")[0]
    #
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
    #
    # print("loading " + sample_file2 + " as sample")
    # data_mutations2 = pd.read_table(sample_file2)
    # data_mutations2["label"] = label_sample2
    #
    # print("loading " + sample_file3 + " as sample")
    # data_mutations3 = pd.read_table(sample_file3)
    # data_mutations3["label"] = label_sample3
    #
    # print("loading " + sample_file4 + " as sample")
    # data_mutations4 = pd.read_table(sample_file4)
    # data_mutations4["label"] = label_sample4
    #
    # print("loading " + sample_file5 + " as sample")
    # data_mutations5 = pd.read_table(sample_file5)
    # data_mutations5["label"] = label_sample5
    #
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

    data = pd.concat([data_mutations0, data_mutations1], sort=False) #, data_mutations2, data_mutations3, data_mutations4,
                      # data_mutations5, data_mutations6, data_mutations7, data_mutations8, data_mutations9] #
    data["passage"] = 0
    data["replica"] = 0

    # data = pd.merge(left=data, right=org_df, on="label")

    toPlot = [label_sample0, label_sample1]# , label_sample2, label_sample3, label_sample4, label_sample5, label_sample6,
              # label_sample7, label_sample8, label_sample9]

    filter_by_coverage_mutation(data, gene_dir, min_read_count=min_coverage)

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
    data = pd.merge(data, consensus_data[['Pos', 'Prev', 'Next', 'Context', 'Consensus', 'label', "major_read_count"]],
                    on=["Pos", "label"])
    data["mutation_type"] = data["Consensus_y"] + data["Base"]
    data["Mutation"] = data["Consensus_y"] + ">" + data["Base"]


    mutation_for_rna = ["AG"]
    dataForPlotting = data[(data["label"].isin(toPlot)) & (data["mutation_type"].isin(mutation_for_rna)) &
                           (data["Rank"] != 0) & (data["Read_count"] > min_coverage)]#& (data["Type"] == "Synonymous")]
    dataForPlotting.to_csv(gene_dir + "/data_XpA_by_mutation.csv", sep=',', encoding='utf-8')

    mutation_for_rna = ["UC"]
    dataForPlotting_UC = data[(data["label"].isin(toPlot)) & (data["mutation_type"].isin(mutation_for_rna)) &
                           (data["Rank"] != 0) & (
                                       data["Read_count"] > min_coverage)]  # & (data["Type"] == "Synonymous")]
    dataForPlotting_UC.to_csv(gene_dir + "/data_UpX_by_mutation.csv", sep=',', encoding='utf-8')

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
