#! /usr/local/python-anaconda-2.7//bin/python

import sys
sys.path.insert(0, '/sternadi/home/volume1/okushnir/scripts/')
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def main():


#RV
    # values = (19.32, 40.52, 45.17)
    # lables = ('RV', 'mRNA', 'rRNA')
    # out_dir = "Z:/volume1/okushnir/Cirseq/RV/20170322_output_all_23_qscore/"
    # clusterd_column(values, lables, out_dir)

#CV
    values = (96.09, 1.66, 1.18)
    lables = ('CV', 'mRNA', 'rRNA')
    out_dir = "Z:/volume1/okushnir/Cirseq/CV/20170426_output_all_23_qscore_01/"
    clusterd_column(values, lables, out_dir)


def clusterd_column(values, labels, out_dir):
    column_number = 3

    ind = np.arange(column_number)  # the x locations for the groups
    width = 0.35  # the width of the bars

    ax = plt.subplot()
    rects1 = ax.bar(ind, values, width, color='DarkOrchid')


    # add some text for labels, title and axes ticks

    ax.set_ylabel('% of Reads', fontsize=16)
    # ax.set_title('Overall Alignment Rate')
    ax.set_xticks(ind)
    ax.set_xticklabels(labels, fontsize=16)
    sns.set_style("darkgrid")
    ax.set_xlim(-0.75, 3)
    ax.set_ylim(0, 100, 20)
    plt.title('Reads Distribution', fontsize=22)

    # Write the value above the bar
    # for rect in rects1:
    #
    #     height = rect.get_height()
    #     ax.text(rect.get_x() + rect.get_width() / 2., 1.05 * height, '%d' % int(height), ha='center', va='bottom')
    #
    path = out_dir + "plots/overall_alignment_rate_test.png"
    plt.savefig(path, dpi=680)
    plt.close('all')


if __name__ == "__main__":
        main()
