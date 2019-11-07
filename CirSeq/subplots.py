#! /usr/local/python-anaconda-3.5//bin/python

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.gridspec as gridspec
import glob
import re


def repeat_len_graphs(results, ax):
    '''makes a repeat length graph'''
    x_values = results.keys()
    x_values.remove("bad")
    y_repeat_len = [results[x]["repeat_length"] for x in x_values]

    medianprops = {'color': 'black', 'linewidth': 2}
    whiskerprops = {'color': 'black', 'linestyle': '-'}

    repeat_plt = ax.boxplot(y_repeat_len, patch_artist=True, medianprops=medianprops, whiskerprops=whiskerprops)
    for box in repeat_plt['boxes']:
        # change outline color
        box.set(color='DarkOrchid', linewidth=2)
        # change fill color
        box.set(facecolor='DarkOrchid')
    ax.set_xlabel("Number of Repeats", fontsize=16)
    ax.set_ylabel("Fragment Size [bp]", fontsize=16)
    labelsy = np.arange(0, 350, 50)
    labelsx = np.arange(1, 11, 1)
    ax.set_yticklabels(labelsy, fontsize=14)
    ax.set_xticklabels(labelsx, fontsize=14)
    sns.set_style("darkgrid")


def clusterd_column(val, ax):
    ''' makes distribution graph'''
    column_number = 3

    ind = np.arange(column_number)  # the x locations for the groups
    width = 0.35  # the width of the bars
    rects1 = ax.bar(ind, val, width, color='DarkOrchid')
    # add some text for labels, title and axes ticks
    ax.set_ylabel('% of Reads', fontsize=16)
    # ax.set_yticklabels()
    # ax.set_title('Overall Alignment Rate')
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(('CV', 'mRNA', 'rRNA'), fontsize=16)
    labelsy = np.arange(0, 50, 10)
    ax.set_yticklabels(labelsy, fontsize=14)
    sns.set_style("darkgrid")
    ax.set_xlim(-0.5, 3)

    # for bar in rects1:
    #     ax.text(ha='center', va='bottom')


def coverage_graph(freqs, ax):
    # show a unified graph otherwise

    data = parse_reads(freqs)
    pos = data[0]
    reads = data[1]
    graph = plt.plot(pos, reads, color="DarkOrchid")

    plt.xlabel("Position In The Genome [bp]", fontsize=16)
    plt.ylabel("Number Of Reads", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    sns.set_style("darkgrid")
    plt.xlim(0, 6550)
    plt.ylim(1000, 1000000)
    plt.yscale("log")


def parse_reads(freqs):
    path = freqs
    df = pd.read_csv(path, sep='\t')

    # remove all duplicates from Pos except the first occurrence
    # remove all x.number duplicates
    df[["Pos"]] = df[["Pos"]].astype(int)
    df = df.drop_duplicates("Pos")

    pos = df["Pos"]  # a vector of all positions
    reads = df["Read_count"]
    return pos, reads


def get_repeats_num(in_dir):
    files = glob.glob(in_dir + "/*.fasta.blast.freqs.stats")
    repeat_summery = {}
    for file in files:
        pattern = re.compile("(\d+\t{1}\d+\n{1})", re.MULTILINE)
        text = open(file, "r").read()
        reads = pattern.findall(text)
        for r in reads:
            key = int(r.split("\t")[0])
            value = int(r.split("\t")[1].split("\n")[0])
            if key in repeat_summery:
                repeat_summery[key] += value
            else:
                repeat_summery[key] = value
    return repeat_summery


def read_repeat_graph(repeat_summery, ax):
    ''' makes read VS. repeat graph'''
    keys = []
    values = []
    for key, val in repeat_summery.items():
        keys.append(key)
        values.append(val)
    graph = plt.bar(keys, values, color='DarkOrchid')
    plt.xlabel("Number of Repeats", fontsize=16)
    plt.ylabel("Number of Reads", fontsize=16)
    plt.xticks(list(range(1, 11)))
    # plt.xlim(min(x)-0.5, max(x)+0.5)
    plt.title('Amount of Reads per Repeat', fontsize=22)
    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    sns.set_style("darkgrid")



#RV
values = (19.32, 40.52, 45.17)
# results = np.load('Z:/volume1/okushnir/Cirseq/RV/20170322_output_all_23_qscore/results.npy').item()
freqs_dir = "Z:/volume1/okushnir/Cirseq/RV/20170322_output_all_23_qscore/tmp" #graph 3
repeats_dict = get_repeats_num(freqs_dir)
freqs_file_path = "Z:/volume1/okushnir/Cirseq/20170322_output_all_23_qscore/RVB14p2.freqs"


gs = gridspec.GridSpec(2, 2)
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1])
ax2 = plt.subplot(gs[2])
ax3 = plt.subplot(gs[3])



clusterd_column(values, ax0)
ax0.set_title('Reads Distribution', fontsize=22)
repeat_len_graphs(results, ax1)
ax1.set_title('Multiple Tandem Repeat Degree', fontsize=22)
read_repeat_graph(repeats_dict, ax2)
ax2.set_title('Amount of reads per repeat', fontsize=22)
coverage_graph(freqs_file_path, ax3)
ax3.set_title('Coverage', fontsize=22)


#CV
# values = (96.09, 1.66, 1.18) #graph 1
# results = np.load(r'Z:\volume1\okushnir\Cirseq\CV\20170426_output_all_23_qscore_01\results.npy').items() #graph 2
# freqs_dir = "Z:/volume1/okushnir/Cirseq/CV/20170426_output_all_23_qscore_01/tmp" #graph 3
# repeats_dict = get_repeats_num(freqs_dir)
# freqs_file_path = "Z:/volume1/okushnir/Cirseq/CV/20170426_output_all_23_qscore_01/CVB3-p2.freqs" #graph 4
#
#
# gs = gridspec.GridSpec(2, 2)
# ax0 = plt.subplot(gs[0])
# ax1 = plt.subplot(gs[1])
# ax2 = plt.subplot(gs[2])
# ax3 = plt.subplot(gs[3])
#
# clusterd_column(values, ax0)
# ax0.set_title('Reads Distribution', fontsize=22)
# repeat_len_graphs(results, ax1)
# ax1.set_title('Multiple Tandem Repeat Degree', fontsize=22)
# read_repeat_graph(repeats_dict, ax2)
# ax2.set_title('Amount of reads per repeat', fontsize=22)
# coverage_graph(freqs_file_path, ax3)
# ax3.set_title('Coverage', fontsize=22)

# Four axes, returned as a 2-d array

# f, axarr = plt.subplots(2, 2)
# repeat_len_graphs(results, axarr[0, 0])
# axarr[0, 0].set_title('Multiple Tandem Repeat Degree')
# clusterd_column(values,  axarr[0, 1])
# axarr[0, 1].set_title('Reads Distribution')
# coverage_graph(freqs_file_path, axarr[1, :1])
# axarr[1, 0].set_title('Coverage')
# coverage_graph(freqs_file_path, axarr[1, 1])
# axarr[1, 1].set_title('Coverage')



plt.show()

