
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
import os.path
import pathlib
import numpy as np
from file_utilities import *
from optparse import OptionParser
import pandas as pd
from scipy.stats import ttest_ind
sns.set_context("talk")




def main():
    # for Cluster
    # parser = OptionParser("usage: %prog [options]")
    # parser.add_option("-f", "--freqs_file_path", dest="freqs_file_path", help="path of the freqs file")
    # parser.add_option("-v", "--virus", dest="virus", help="Virus name: CVB3 for CV; RVB14 for RV; PV for PV")
    # (options, args) = parser.parse_args()
    #
    # freqs_file = options.freqs_file_path
    # virus = options.virus
    #
    # freqs_file = check_filename(freqs_file)

    #for Local
    sample1 = "RV-p11"
    suffix1 = "%s.with.mutation.type.freqs" % sample1
    freqs_file1 = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/%s/q30_3UTR_new/%s" % \
                 (sample1, suffix1)

    sample2 = "RV-p12"
    suffix2 = "%s.with.mutation.type.freqs" % sample2
    freqs_file2 = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/%s/q30_3UTR_new/%s" % \
                 (sample2, suffix2)
    sample3 = "RV-IVT"
    suffix3 = "%s.with.mutation.type.freqs" % sample3
    freqs_file3 = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/%s/q30_3UTR_new/%s" % \
                 (sample3, suffix3)
    virus = "RVB14"
    seq_meth = "AccuNGS"
    freqs_list = [freqs_file2, freqs_file1, freqs_file3]


    # 1. Get freqs file and CirSeq running directory.
    path = freqs_file1.split('/')[0:-1]
    out_dir = '' #the freqs directory
    for i in path:
        out_dir += str(i + '/')
    tmp_cirseq_dir = out_dir + 'tmp/'
    pathlib.Path(out_dir + 'plots/').mkdir(parents=True, exist_ok=True)
    out_plots_dir = out_dir + 'plots/'

    if virus == "CVB3":
        ncbi_id ="M16572"
    if virus == "RVB14":
        ncbi_id = "NC_001490"
    if virus == "PV":
        ncbi_id ="V01149"


    """ 3. Adding mutation types to the freqs file"""
    # for i in freqs_list:
    #     if not os.path.isfile(i[0:-5] + "with.mutation.type.freqs"):
    #          append_mutation = find_mutation_type(i, ncbi_id)
    #     freqs_file_mutations = i[0:-5] + "with.mutation.type.freqs"

    fig2 = plt.figure(figsize=(16, 9))
    gs = gridspec.GridSpec(4, 1)
    ax0 = plt.subplot(gs[0, 0])

    gs.tight_layout(fig2)
    # ax = plt.subplot()

    # make_boxplot_sample_control(freqs_list, ax0)
    #
    # plt.savefig(out_plots_dir + suffix1.split(sep='.')[0] + '_Transitions_Report_Control_sample_strech_new_1.png',
    #             dpi=600)
    # plt.close("all")
    make_boxplot_mutation_type(freqs_list, ax0)
    plt.savefig(out_plots_dir + suffix1.split(sep='.')[0] + '_Transitions_Report_Control_sample_silent_mis_non.png',
                dpi=600)
    plt.close("all")
    print("The Transition Plots is ready in the folder")

"""Graphs"""

def arrange_freqs_file(freqs_file):
    """
    Plots the mutation frequencies boxplot
    :param freqs_file: pandas DataFrame after find_mutation_type function
    :param ax: ax location
    :return:
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
    min_read_count = 100000
    data = data[data['Pos'] == np.round(data['Pos'])]  # remove insertions
    data['Pos'] = data[['Pos']].apply(pd.to_numeric)
    data = data[data['Read_count'] > min_read_count]
    data['mutation_type'] = data['Ref'] + data['Base']
    data = data[data['Ref'] != data['Base']]
    data = data[data["Base"] != "-"]
    data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
    data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
    data["Mutation"] = data["Ref"] + "->" + data["Base"]
    data["Source"] = freqs_file.split('/')[-1].split('.')[0]
    return data

def make_boxplot_sample_control(freqs_list, ax):
    data = pd.DataFrame()
    for k, i in enumerate(freqs_list):
        data_freqs = arrange_freqs_file(i)
        data = pd.concat([data_freqs, data], axis=0)
        data["Source"].replace("RV-p11", "Replica #1", inplace=True)
        data["Source"].replace("RV-p12", "Replica #2", inplace=True)
        data["Source"].replace("RV-IVT", "RNA Control", inplace=True)
    data.to_csv("/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/RV-p11/"
                "q30_3UTR_new/data_Mutation_type.csv", sep=',', encoding='utf-8')

    # By Mutation
    g1 = sns.factorplot(x="Mutation Type", y="Frequency", data=data, col=" ", hue="Source", order=["Synonymous",
                        "Non-Synonymous", "Premature Stop Codon"], col_order=["A->G", "U->C", "G->A", "C->U"],
                        kind="box", color=sns.color_palette("Paired", 12)[7], legend_out=True)#, col_wrap=2)
    new_title = ''
    g1._legend.set_title(new_title)

    titles = ["A->G", "U->C", "G->A", "C->U"]
    for ax, title in zip(g1.axes.flat, titles):
        ax.set_title(title)
    g1.set_xticklabels(["Silent", "Missense", "Nonsense"], fontsize=10)
    g1.set(yscale="log", ylim=(0.00001, 0.1))
    g1.set_xlabels('')
    handles = g1._legend_data.values()
    labels = g1._legend_data.keys()
    g1.fig.legend(handles=handles, labels=labels, loc='center right', ncol=1)


def make_boxplot_mutation_type(freqs_list, ax):
    data = pd.DataFrame()
    for k, i in enumerate(freqs_list):
        data_freqs = arrange_freqs_file(i)
        data = pd.concat([data_freqs, data], axis=0)
        data["Source"].replace("RV-p11", "Replica #1", inplace=True)
        data["Source"].replace("RV-p12", "Replica #2", inplace=True)
        data["Source"].replace("RV-IVT", "RNA Control", inplace=True)
    data.to_csv("/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/RV-p11/"
                "q30_3UTR_new/data_Mutation_type.csv", sep=',', encoding='utf-8')


    g1 = sns.factorplot(x="Source", y="Frequency", data=data, col="Mutation Type", hue="Mutation",
                        order=["Replica #1", "Replica #2", "RNA Control"],
                        col_order=["Synonymous", "Non-Synonymous", "Premature Stop Codon"],
                        hue_order=["A->G", "U->C", "C->U", "G->A"],
                        kind="box") #color=sns.color_palette("Paired", 12)[8])
    titles = ["Silent", "Missense", "Nonsense"]
    for ax, title in zip(g1.axes.flat, titles):
        ax.set_title(title)
    g1.set_xticklabels(["Replica #1", "Replica #2", "RNA Control"], fontsize=10)
    g1.set(yscale="log", ylim=(0.00001, 0.1))
    g1.set_xlabels('')
    new_title = ''
    g1._legend.set_title(new_title)

if __name__ == "__main__":
    main()
