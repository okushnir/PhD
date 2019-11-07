#! /usr/local/python_anaconda/bin/python3.4

"""
@Author: odedkushnir

"""

import sys
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



def main():

    freqs_position_plot(r"C:\Users\Oded\Desktop\cooccuring_mutation.txtfreq.txt", r"C:\Users\Oded\Desktop\plot.png")
    freqs_position_plot(r"C:\Users\Oded\Desktop\CVB3-p2.edited.freqs", r"C:\Users\Oded\Desktop\plot_ori_freqs.png")

def freqs_position_plot(file_name, out_file):
    with open(file_name, "r") as freq:
        data = freqs_to_dataframe(freq)
        print("Data in pandas DataFrame")
        data["Mutation"] = data["Ref"] + "->" + data["Base"]
        print(data)
        fg = sns.FacetGrid(data=data, hue="Mutation", size=12, aspect=2, palette="husl")
        fg = (fg.map(plt.scatter, "Pos", "Freq").add_legend())
        plt.yscale('log')
        plt.ylim(10**-6, 1)
        # for i in data["Mutation"]:
        #     mutation_dict[i] = mutation_dict.get(i, 0) + 1
        print("plot is ready")
        plt.savefig(out_file, dpi=300)
        plt.close("all")


def freqs_to_dataframe(freqs_file):
    data = pd.read_table(freqs_file)
    data = data[data.Ref != '-']
    data = data[data.Base != '-']
    data.reset_index(drop=True, inplace=True)
    return data

if __name__ == "__main__":
    main()