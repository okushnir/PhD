#!/powerapps/share/python-anaconda-3.2019.7/bin/python

"""
@Author: odedkushnir

"""

import sys, argparse
import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re
import glob
import Bio.Seq as Seq
import matplotlib.gridspec as gridspec
import time


def plot_generation_time(file_name1):#1, file_name2):
    """
    :param file_name:
    :return: plot and data frame
    """
    data = pd.read_csv(file_name1)
    # data2 = pd.read_csv(file_name2)
    # data = pd.concat([data1, data2])
    data = data.rename(columns={"Cq Mean": "Cq_Mean"})
    data = data.rename(columns={"Biological Set Name": "Treatment"})

    #Amicon
    # data = data1
    # data["Sample"] = data["Sample"].apply(lambda x: x.replace("Amicon ", "Amicon\n"))
    # data["Sample"] = np.where(data["Sample"] == "NUCase+", "RNase25\nDNase10U", data["Sample"])
    # sample_order = ["Amicon\nViralRNA", "Amicon\nQuickRNA", "NUCase-", "RNase25\nDNase10U"]

    #20191126
    # data = data.loc[(data.Sample == "RNase25")|(data.Sample == "NUCase-")|(data.Sample == "RNase100")|
    #                 (data.Sample == "RNase250")|(data.Sample == "RNase500")]
    # data = data.loc[(data.Cq_Mean > 0)]
    # data["Sample"] = np.where(data["Sample"] == "RNase25", "RNase25\nDNase10U", np.where(data["Sample"] == "RNase100",
    #                  "RNase100\nDNase40U", np.where(data["Sample"] == "RNase250", "RNase250\nDNase625U", np.where(
    #         data["Sample"] == "RNase500", "RNase500\nDNase300U", "NUCase-"))))
    # sample_order = ["NUCase-", "RNase25\nDNase10U", "RNase100\nDNase40U", "RNase250\nDNase625U", "RNase500\nDNase300U"]

    #20191205
    data = data.loc[(data.Sample == "RNase62.5") | (data.Sample == "RNase125") |
                    (data.Sample == "RNase250") | (data.Sample == "RNase500") | (data.Sample == "DNase6.25") |
                    (data.Sample == "DNase12.5") | (data.Sample == "DNase25") | (data.Sample == "DNase50") |
                    (data.Sample == "NUCase-")]

    data = data.loc[(data.Cq_Mean > 0)]
    unit_regex = re.compile(r'[\d][\d][.][\d]|[\d][\d][\d]|[\d][.][\d][\d]|[\d][\d]')
    data["Units"] = data["Sample"].apply(lambda x: unit_regex.search(x).group() if x != "NUCase-" else 0)
    data["Units"] = data["Units"].astype(float)

    sample_order = ["NUCase-", "DNase6.25", "DNase12.5", "DNase25", "DNase50", "RNase62.5", "RNase125", "RNase250", "RNase500"]

    target_order = ["RV", "GAPDH", "rRNA"]

    g0 = sns.pointplot(x="Sample", y="Cq", data=data, hue="Target", order=sample_order, hue_order=target_order)
    g0.legend(loc=2)
    g0.set_xticklabels(g0.get_xticklabels(), rotation=30)
    output_dir = (file_name1.split("/")[0:-1])
    output_dir = '/'.join(output_dir)
    plt.savefig(output_dir + "/Cq_one_plot.png", dpi=300)

    g1 = sns.relplot(x="Units", y="Cq", data=data, hue="Target", hue_order=target_order, kind="line", col="Treatment",
                     facet_kws={"sharex": False})
    g1.axes.flat[0].set_xlabel("mg")
    g1.axes.flat[1].set_xlabel("Units")
    # g1.legend(loc=2)
    # g1.set_xticklabels(g1.get_xticklabels(), rotation=30)

    # g1.set(xlabel="Time[hr]", ylabel="Copy number of RNA")
    # g1.set_yscale("log")
    output_dir = (file_name1.split("/")[0:-1])
    output_dir = '/'.join(output_dir)
    plt.savefig(output_dir + "/Cq.png", dpi=300)
    return data

# def main(args):
def main():
    # Amicon
    # file_name1 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/RV/RVB14/RealTime Results xlsx/20191124 capsid/Oded_2019-11-24 12-39-02_BR003827 -  Quantification Cq Results_0.csv"
    # plot_generation_time(file_name1)

    # file_name1 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/RV/RVB14/RealTime Results xlsx/20191107 capsid/Oded_2019-11-07 12-39-22_BR003827 -  Quantification Cq Results_0.csv"
    # file_name2 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/RV/RVB14/RealTime Results xlsx/20191126 capsid/Oded_2019-11-26 17-43-48_BR003827 -  Quantification Cq Results_0.csv"
    # plot_generation_time(file_name1, file_name2)

    # file_name1 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/RV/RVB14/RealTime Results xlsx/20191107 capsid/Oded_2019-11-07 12-39-22_BR003827 -  Quantification Cq Results_0.csv"
    file_name2 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/RV/RVB14/RealTime Results xlsx/20191205 capsid/Oded_2019-12-05 10-35-18_BR003827 -  Quantification Cq Results_0.csv"
    plot_generation_time(file_name2)
if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    # parser.add_argument("virus", type=str, help="name of the virus, RVB14")
    # parser.add_argument("passages", type=str, help="from which passages, p0-p12")
    # parser.add_argument("without", type=int, help="Exclude passage no.")
    # parser.add_argument("input_dir", type=str, help="the path to the directory that contains data_mutation.csv")
    # parser.add_argument("quality", type=str, help="what is the prefix for the data_mutation.csv file; quality of the pipline ; for example: q38")
    # parser.add_argument("mutation_type", type=str, help="all/syn")
    # args = parser.parse_args(sys.argv[1:])
    # main(args)
    main()