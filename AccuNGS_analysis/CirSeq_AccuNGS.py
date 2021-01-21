#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import sys, argparse
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_theme("paper")
sns.set_style("ticks")
sns.despine()
# sns.set(font_scale=1.1)



def main():

    virus_table = pd.read_csv("/Users/odedkushnir/Google Drive/Studies/PhD/Prgress reports/20180730 Commite/Presentation/Workbook1.csv")
    virus_table = virus_table[["Method", "Starting Material", "Sequenced Viral particles"]]
    method_data = pd.melt(virus_table, id_vars="Method")
    method_data = method_data.rename(columns={"value": "Estimated Viral Particles"})
    print(method_data)
    with sns.axes_style("ticks"):
        g1 = sns.barplot(data=method_data, x="Method", y="Estimated Viral Particles", hue="variable")
        g1.set(yscale='log')
        g1.legend(title="")
        # plt.legend(bbox_to_anchor=(0.5, 1.2), loc='upper left')
        # plt.tight_layout()
        plt.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/Thesis/Figs/cirseq_accungs.png", dpi=300)

if __name__ == "__main__":
    main()