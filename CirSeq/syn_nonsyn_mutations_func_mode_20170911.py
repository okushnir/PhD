"""
@Author: odedkushnir

"""

import pandas as pd
import numpy as np
import time
from Bio.Seq import Seq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Entrez
from Bio import SeqIO
import os.path
import pathlib
from optparse import OptionParser
import file_utilities



start_time = time.time()


def main():
    # For Local Run
    path = "/Volumes/STERNADILABTEMP$/volume1/okushnir/Cirseq/PV/Mahoney/P3/20170907_q20r3_blastn/PV-p3.1036617.freqs"
    virus = "PV"

    # file_name = "PV-p3.1036617.freqs"
    # virus = file_name.split(sep='-')[0]
    # virus += file_name.split(sep='-')[1].split(sep='.')[0]

    # For Cluster Run
    # parser = OptionParser("usage: %prog [options]")
    # parser.add_option("-f", "--freqs_file_path", dest="freqs_file_path", help="path of the freqs file")
    # parser.add_option("-v", "--virus", dest="virus", help="Virus name: CVB3 for CV; RVB14 for RV; PV for PV")
    # (options, args) = parser.parse_args()
    # path = options.freqs_file_path
    # virus = options.virus

    path_list = path.split(sep="-")[0].split(sep="/")
    plots_path = ''
    for i in range(0, len(path_list) - 1, 1):
        plots_path += path_list[i] + "/"



    if virus == "CVB3":
        ncbi_id ="M16572"
    if virus == "RVB14":
        ncbi_id = "NC_001490"
    if virus == "PV":
        ncbi_id ="V01149"

    if not os.path.isfile(path + ".with.mutation.type.func2.freqs"):
        append_mutation = find_mutation_type(path, ncbi_id)


    mutation_file = path + ".with.mutation.type.func2.freqs"
    mutation_rates = freqs_to_dataframe(mutation_file)

    pathlib.Path(plots_path + 'plots/').mkdir(parents=True, exist_ok=True)

    df = make_boxplot_mutation(mutation_rates, plots_path, virus)
    make_boxplot_mutation_median(df, plots_path, virus)




def find_mutation_type(freqs_file, ncbi_id):
    """
    This function adds Mutation type to the freqs file
    :param freqs_file:  The path of the relevant freqs file
    :return:DataFrame of the frqs file with mutation type column, save it in txt file
    """
    file_name = freqs_file
    data = freqs_to_dataframe(freqs_file)
    start_pos, end_pos = find_coding_region(ncbi_id)
    data = data.loc[data['Pos'] >= start_pos]
    check_12mer(data)
    data['Codon'] = ""
    data['Ref_Protein'] = ""
    data['Potential_Protein'] = ""
    data['Mutation Type'] = ""
    data['Pos'] = data['Pos'].astype(int)
    data.reset_index(drop=True, inplace=True)
    for kmer in range(data.first_valid_index(), (len(data)), 12):
        print ("going over kmer: " + str(kmer) + "/" + str((len(data))))
        kmer_df = data[:].loc[kmer:kmer+11]
        print("finding codons....")
        find_codon(kmer_df)
        print("translating codons....")
        translate_codon(kmer_df)
        print("Sets the Mutation type....")
        kmer_df['Mutation Type'] = kmer_df[['Ref_Protein', 'Potential_Protein']].apply(
            lambda protein: check_mutation_type(protein[0], protein[1]), axis=1)
    print("After a long for loop")
    file_name += ".with.mutation.type.func2.freqs"
    data.to_csv(file_name, sep='\t', encoding='utf-8')
    print("The File is ready in the folder")
    print("--- %s sec ---" % (time.time() - start_time))
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
    :param ncbi_id:
    :return:
    """
    try:
        Entrez.email = "A.N.Other@example.com"
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=ncbi_id)
        ncbi_gb = SeqIO.read(handle, "gb")
        handle.close()
        location = list(ncbi_gb.features[1].location)
        start_pos = location[0]
        end_pos = location[-4]
        return start_pos, end_pos
    except:
        print('Failed to fetch record.')


def check_12mer(data):
    """
      :param data: pandas DataFrame of the freqs file
      :return: None - checks that DataFrame is in 12 kmer
      """
    if (len(data)) % 12 != 0:
        print('The data is not in 12 mer. arrange the data')
    else:
        print('The data is in 12 mer. All OK.')
        print('The length of data is:', (len(data)))


def find_codon(data):
    """
    :param data: pandas DataFrame of the freqs file
    :return: pandas DataFrame with codon column
    """
    data['Codon'] = ""
    x = 4
    first_index_kmer = data.first_valid_index()
    for i in range(first_index_kmer, first_index_kmer+x, 1):
        data['Codon'][i] = data["Base"].loc[i] + data["Ref"].loc[first_index_kmer+x] + data["Ref"].loc[(first_index_kmer+(x*2))]
    for i in range((first_index_kmer+x), (first_index_kmer+(x*2)), 1):
            data['Codon'][i] = data["Ref"].loc[first_index_kmer] + data["Base"].loc[i] + data["Ref"].loc[first_index_kmer+(x*2)]
    for i in range((first_index_kmer+(x*2)), (first_index_kmer+(x * 3)), 1):
            data['Codon'][i] = data["Ref"].loc[first_index_kmer] + data["Ref"].loc[first_index_kmer+x] + data["Base"].loc[i]
    return data


def translate(seq):
    """
    :param seq: string of ATGC
    :return: protein sequence of seq
    """
    seq = Seq(str(seq))
    protein = seq.translate()
    return protein


def translate_codon(data):
    """
    :param data: pandas DataFrame of the freqs file
    :return: pandas DataFrame with Reference protein column and Potential protein column
    """
    data['Ref_Protein'] = ""
    data['Potential_Protein'] = ""
    protein_ref = Seq(data["Codon"][data.first_valid_index()]).translate()
    data['Ref_Protein'].loc[data.first_valid_index():data.first_valid_index()+11] = protein_ref
    data['Ref_Protein'].replace('(', '', inplace=True)
    data['Ref_Protein'].replace(')', '', inplace=True)

    data["Potential_Protein"] = data["Codon"].apply(translate)
    data["Potential_Protein"].replace('(', '', inplace=True)
    data["Potential_Protein"].replace(')', '', inplace=True)
    return data


def check_mutation_type(protein1, protein2):
    """
    :param protein1: amino acid 1
    :param protein2: amino acid 2
    :return: The mutation type
    """
    Mutation_Type = ""
    if protein1 == protein2:
        Mutation_Type = "Synonymous"
    elif protein1 != protein2:
        Mutation_Type = "Non-Synonymous"
    if protein2 == '*':
        Mutation_Type = "Premature Stop Codon"
    return Mutation_Type

# from Maoz
def make_boxplot_mutation(data, out_dir, virus_name):
    """
        :param data: pandas DataFrame after find_mutation_type function
        :return: pandas DataFrame ready for plotting
        """
    data['Base'].replace('T', 'U', inplace=True)
    data['Ref'].replace('T', 'U', inplace=True)
    min_read_count = 100000
    data = data[data['Pos'] == np.round(data['Pos'])]  # remove insertions
    data['Pos'] = data[['Pos']].apply(pd.to_numeric)
    data = data[data['Read_count'] > min_read_count]
    data['mutation_type'] = data['Ref'] + data['Base']
    data = data[data['Ref'] != data['Base']]
    data = data[data["Base"] != "-"]
    # data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
    # data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
    data["Mutation"] = data["Ref"] + "->" + data["Base"]
    sns.set_palette(sns.color_palette("Paired", 12))
    g1 = sns.boxplot(x="Mutation Type", y="Freq", hue="Mutation", data=data,
                     hue_order=["C->U", "U->C", "G->A", "A->G", "C->A", "G->U", "U->G", "U->A", "G->C", "A->C",
                                "A->U", "C->G"], order=["Synonymous", "Non-Synonymous", "Premature Stop Codon"])
    g1.set(yscale="log")
    sns.set_style("darkgrid")
    plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0., fontsize=7)
    g1.set_ylim(10 ** -6, 1)
    g1.tick_params(labelsize=7)
    plt.title(virus_name + ' Mutations Frequencies', fontsize=22)
    plt.savefig(out_dir + "plots/freqs_type_%s.png" % str(min_read_count), dpi=300)
    plt.close("all")
    print("The Plot is ready in the folder")
    return data


def make_boxplot_mutation_median(data, out_dir, virus_name):
        """
            :param data: pandas DataFrame after find_mutation_type function
            :return: pandas DataFrame ready for plotting
            """
        min_read_count = 100000
        non_syn = data[data['Mutation Type'] == 'Non-Synonymous']
        # non_syn = non_syn[non_syn['Pos'] > 1000]
        # non_syn = non_syn[non_syn['Pos'] < 2800]
        g1 = sns.boxplot(x="Mutation", y="Freq", data=non_syn,
                         order=["A->C", "A->G", "A->U", "C->A", "C->G", "C->U", "G->A", "G->C", "G->U", "U->A",
                                "U->C", "U->G"], color="DarkOrchid")
        sns.set_style("darkgrid")
        grouped_df = non_syn.groupby(["Mutation"])["Freq"]
        for key, item in grouped_df:
            print(grouped_df.get_group(key), "\n\n")

        medians = grouped_df.median()
        print(medians)
        median_labels = [str(np.round(s, 6)) for s in medians]
        pos = range(len(medians))
        for tick, label in zip(pos, g1.get_xticklabels()):
            g1.text(pos[tick], medians[tick]*1.01, median_labels[tick],
                    horizontalalignment='center', size='x-small', color='w', weight='semibold', fontsize=5)

        len_col = grouped_df.size()
        for tick, label in zip(pos, g1.get_xticklabels()):
            g1.text(pos[tick], medians[tick]+0.5, len_col[tick],
                    horizontalalignment='center', size='x-small', color='navy', weight='semibold', fontsize=5)
        g1.set(yscale="log")
        plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0., fontsize=7)
        g1.set_ylim(10 ** -6, 1)
        g1.tick_params(labelsize=7)
        plt.title(virus_name + ' Non-Synonymous Mutations Frequencies', fontsize=19)
        plt.savefig(out_dir + "plots/non_syn_median_%s.png" % str(min_read_count), dpi=300)
        plt.close("all")
        print("The Plot is ready in the folder")


if __name__ == "__main__":
    main()










