#! /usr/local/python/anaconda_python-3.6.1


import pandas as pd
import numpy as np
import scipy.stats as stats
import time
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import Entrez
from Bio import SeqIO
from Utilities import cirseq_utilities
import re
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns; sns.set_style("ticks")
# sns.set_style("ticks")
from statsmodels.graphics.mosaicplot import mosaic

def make_consensus(data):
    """

    :param dataframe:
    :return:
    """

    consensus_data = data[['Pos', 'Base', "Read_count", "Freq"]][data['Rank'] == 0]
    consensus_data["Base"].replace('T', 'U', inplace=True)

    consensus_data = consensus_data[consensus_data['Pos'] == np.round(consensus_data['Pos'])]  # removes insertions
    consensus_data['Prev'] = consensus_data.shift(1)['Base'] + consensus_data['Base']
    consensus_data['Next'] = consensus_data['Base'] + consensus_data.shift(-1)['Base']
    consensus_data['Pos'] = consensus_data[['Pos']].apply(pd.to_numeric)
    consensus_data["Context"] = consensus_data.shift(1)['Base'] + consensus_data['Base'] + consensus_data.shift(-1)[
        'Base']
    consensus_data["major_read_count"] = np.round(consensus_data['Read_count'] * consensus_data['Freq'])
    consensus_data = consensus_data.rename(columns={"Base": "Consensus"})
    return consensus_data


def find_upstream_base(record_dict, node_name, pos):
    """

    :param record_dict:
    :param node_name:
    :param pos:
    :return:
    """
    seq_record = record_dict.get(node_name)
    try:
        base = seq_record.seq[pos-2]
    except Exception:
        return None
    return base



def parasel_context(parasel_output_path, fasta_file):
    """
    :param parasel_output_path:
    :return:
    """
    record_dict = {}
    with open(fasta_file, 'r') as fasta_obj:
        record_dict = SeqIO.to_dict(SeqIO.parse(fasta_obj, "fasta"))
    data = pd.read_table(parasel_output_path)
    data["letter1"].replace(0, 'A', inplace=True)
    data["letter1"].replace(1, 'C', inplace=True)
    data["letter1"].replace(2, 'G', inplace=True)
    data["letter1"].replace(3, 'U', inplace=True)

    data["letter2"].replace(0, 'A', inplace=True)
    data["letter2"].replace(1, 'C', inplace=True)
    data["letter2"].replace(2, 'G', inplace=True)
    data["letter2"].replace(3, 'U', inplace=True)

    data["BaseToBase"] = data["letter1"] + data["letter2"]
    data["Mutation"] = data["letter1"] + "->" + data["letter2"]

    data["Filter"] = data.node_name.str.contains('^B')
    data = data[data["Filter"] == True]

    data["Prev"] = data.apply(lambda x: find_upstream_base(record_dict, x["node_name"], x["Pos(query)"]), axis=1)
    data.dropna(subset=["Prev"], inplace=True)
    data["Prev"].replace('T', 'U', inplace=True)
    data["Filter_bases"] = data.Prev.str.contains("A") | data.Prev.str.contains("C") | data.Prev.str.contains("G") | \
                         data.Prev.str.contains("U")
    data = data[data["Filter_bases"] == True]

    data["Context"] = data["Prev"] + data["letter1"]
    data["Context"].replace('T', 'U', inplace=True)
    return data


def checkKey(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        raise Exception()

def main():
    parasel_file = "/Users/odedkushnir/Projects/signatures/ADAR/Phylogenetic/RVB/output/parasel_output.self.txt"
    fasta_file = "/Users/odedkushnir/Projects/signatures/ADAR/Phylogenetic/RVB/Ann.full13.B.NAMES.txt"

    # find_upstream_base(record_dict, "B_79_FJ445155", 4291)


    parasel_df = parasel_context(parasel_file, fasta_file)

    # # For A->G Mutations
    # parasel_df_ag = parasel_df[parasel_df["letter1"] == "A"]
    #
    # parasel_df_ag["Filter_AG"] = parasel_df_ag.BaseToBase.str.contains("AG") | parasel_df_ag.BaseToBase.str.contains("AA")
    # parasel_df_ag = parasel_df_ag[parasel_df_ag["Filter_AG"] == True]
    #
    # df_list = [parasel_df, parasel_df_ag]
    # for i, df in enumerate(df_list):


        # df_list[i] = df
    # print(df_list)

    parasel_df_ag = parasel_df[parasel_df["letter1"] == "A"]

    parasel_df_ag["Filter_AG"] = parasel_df_ag.BaseToBase.str.contains("AG") | parasel_df_ag.BaseToBase.str.contains("AA")
    parasel_df_ag = parasel_df_ag[parasel_df_ag["Filter_AG"] == True]

    #save to dataframes to csv
    # df_list[0].to_csv("/Users/odedkushnir/Projects/signatures/ADAR/Phylogeny/parasel_all_mutations_self.csv", sep=',',
    #                   encoding='utf-8')
    # df_list[1].to_csv("/Users/odedkushnir/Projects/signatures/ADAR/Phylogeny/parasel_ag_self.csv", sep=',',
    #                   encoding='utf-8')



    # parasel_df_ag = df_list[1]
    # parasel_df = df_list[0]



    parasel_df_ag["ADAR_Context"] = parasel_df_ag.Context.str.contains('UA') | parasel_df_ag.Context.str.contains('AA')
    parasel_df_ag["A>G"] = parasel_df_ag.Mutation.str.contains("A->G")
    # parasel_df_ag = parasel_df_ag[parasel_df_ag["Filter_adar"] == True]
    # print(parasel_df_ag.to_string)

    # extreme1 = pd.read_csv("/Users/odedkushnir/Projects/signatures/ADAR/Phylogenetic/RVB/output/extrem1.csv")
    # extreme2 = pd.read_csv("/Users/odedkushnir/Projects/signatures/ADAR/Phylogenetic/RVB/output/extrem2.csv")
    # extreme1.set_index("Mutation", inplace=True)
    # extreme2.set_index("Mutation", inplace=True)
    # print(extreme1.to_string())
    # oddsratio1, pval1 = stats.fisher_exact(extreme1, alternative="greater")
    # print(oddsratio1, pval1)
    #
    # print("extreme2")
    # print(extreme2.to_string())
    # oddsratio2, pval2 = stats.fisher_exact(extreme2, alternative="greater")
    # print(oddsratio2, pval2)
    #
    crossta_df = pd.crosstab(index=parasel_df_ag["ADAR_Context"], columns=parasel_df_ag["A>G"])
    # crossta_df.to_csv("/Users/odedkushnir/Projects/signatures/ADAR/Phylogeny/Crosstab_df.csv", sep=',',
    #                   encoding='utf-8')

    print(crossta_df.to_string())
    # # parasel_df_ag_asen = parasel_df_ag.sort_values(by=["Mutation", "ADAR_Context"], ascending=True)

    oddsratio, pval = stats.fisher_exact(crossta_df, alternative="less")
    print(oddsratio, pval)

    # oddsratio, pvalue = stats.fisher_exact([[8, 2], [1, 5]], alternative="less")
    # print(oddsratio, pvalue)
    #

    # data = parasel_df_ag[parasel_df_ag["ADAR_Context"] == True]
    # data = pd.DataFrame({"ADAR_Context": parasel_df_ag["ADAR_Context"], "Mutation": parasel_df_ag["Mutation"]})
    # print(data.to_string())
    # mosaic(data, ["ADAR_Context", "Mutation"])

#y-axis from 60000
    # crossta_df = crossta_df[crossta_df.index == "A->G"]
    # sns.set_context("poster")

    g = crossta_df.plot.bar(stacked=False, color=("#998ec3", "#f1a340"))
    g.set(yscale="log")
    # g.get_legand().set_title()
    # g.set_ylim(60000, 1.1*10**5)
    plt.yticks(fontsize=18)
    plt.xticks(rotation=0, fontsize=18)
    plt.ylabel("Counts", fontsize=20)
    plt.xlabel("ADAR-like\nContext", fontsize=20)
    plt.legend(title="A>G", fontsize=12)
    sns.despine()
    # plt.legend(loc='center left', bbox_to_anchor=(1.03, 1))
    # plt.show()
    plt.tight_layout()
    plt.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/AG_Context_parasel_result_fig4.png",
                dpi=300)
    # plt.savefig("/Users/odedkushnir/Projects/signatures/ADAR/Phylogenetic/RVB/plots/AG_Context_parasel_result_fig3.png", dpi=300)



#breaking y axis
    # fig, axes = plt.subplots(2, sharex=True)
    # p1 = crossta_df.plot.bar(stacked=True, color=("#998ec3", "#f1a340"), ax=axes[0]) #top
    # p2 = crossta_df.plot.bar(stacked=True,color=("#998ec3", "#f1a340"),  ax=axes[1]) #bottom
    # axes[0].set_ylim(ymin=60000)
    # axes[1].set_ylim(ymax=5000)
    #
    # axes[0].spines['bottom'].set_visible(False)
    # axes[1].spines['top'].set_visible(False)
    # # axes[0].set_xticks([])
    # # axes[0].xaxis.tick_top()
    # axes[0].tick_params(axis='x',          # changes apply to the x-axis
    # which='both',      # both major and minor ticks are affected
    # bottom=False,      # ticks along the bottom edge are off
    # top=False,         # ticks along the top edge are off
    # labelbottom=False)
    # axes[1].xaxis.tick_bottom()
    #
    # axes[1].get_legend().remove()
    # # axes[0].legend(title="ADAR Context")
    # # fig.suptitle("ParaSel Result", fontsize=14)
    # # plt.legend(title="ADAR Context")
    # axes[0].set_ylabel("Count")
    # # axes[1].set_ylabel("Count")
    #
    # plt.xticks(rotation=0)
    # # plt.tight_layout()
    # plt.show()
    # plt.savefig("/Users/odedkushnir/Projects/signatures/ADAR/Phylogeny/plots/AG_Context_parasel_result.png", dpi=300)

    #Upstream nuc counter
    # fisher_df = parasel_df_ag.groupby(["Context"]).count()
    # columns= ["Pos(query)"]
    # fisher_df = pd.DataFrame(fisher_df, columns=columns)
    # fisher_df.rename(columns={'Pos(query)': 'Count'}, inplace=True)
    # print(fisher_df.to_string())

    # merge_df = pd.merge(left=rv_p71_df, right=parasel_df, left_on="Pos", right_on="Pos", how="outer")
    #
    # merge_df.dropna(subset=["node_id"], inplace=True)
    # merge_df = merge_df[merge_df["Consensus"] == "A"]

    # grouped = parasel_df.groupby('Context', as_index=False).agg("count")
    # contex_df = pd.DataFrame(grouped)
    # contex_df = contex_df[['Context', 'Pos']]
    # contex_df = contex_df.rename(columns={"Pos": "Count"})
    # contex_df["Sum"] = contex_df["Count"].sum()
    # contex_df["Fraction"] = contex_df.apply(lambda x: (x["Count"] / x["Sum"])*100, axis=1)
    # print(contex_df)
    # contex_df.to_csv("/Users/odedkushnir/Projects/signatures/ADAR/Phylogeny/parasel_precentage.csv", sep=',',
    #                  encoding='utf-8')

    #bases percentage
    # raw_data = pd.read_table("/Users/odedkushnir/Projects/signatures/ADAR/Phylogeny/parasel_output.txt")
    # raw_data["Filter"] = raw_data.node_name.str.contains('^B')
    # raw_data = raw_data[raw_data["Filter"] == True]
    # grouped_raw = raw_data.groupby('letter1', as_index=False).agg("count")
    # grouped_df = pd.DataFrame(grouped_raw)
    # grouped_df = grouped_df[['letter1', 'Pos']]
    # grouped_df = grouped_df.rename(columns={"Pos": "Count"})
    # grouped_df["Sum"] = grouped_df["Count"].sum()
    # grouped_df["Fraction"] = grouped_df.apply(lambda x: (x["Count"] / x["Sum"])*100, axis=1)
    # grouped_df.to_csv("/Users/odedkushnir/Projects/signatures/ADAR/Phylogeny/bases_precentage.csv", sep=',', encoding='utf-8')


if __name__ == "__main__":
    main()