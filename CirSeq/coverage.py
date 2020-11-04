
"""
@Author: daniellem1

"""

'''preprocess .freqs files in order to get for each genome position it's num of reads'''

import os.path
import os
import pandas as pd
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("ticks")
sns.despine()


def coverage_graph(freqs, out_dir):
    sample = freqs.split("/")[-1].split(".")[0]
    df = pd.read_table(freqs, sep="\t")

    # remove all duplicates from Pos except the first occurrence
    # remove all x.number duplicates
    df[["Pos"]] = df[["Pos"]].astype(int)
    df = df.drop_duplicates("Pos")
    df = df[df["Base"] != "-"]
    df = df[df["Read_count"] > 1000]
    med_reads = df["Read_count"].median()
    print("%s: %s" % (sample, str(med_reads)))

    graph = sns.lineplot(x="Pos", y="Read_count", data=df)

    graph.set_xlabel("Position In The Genome[bp]")
    graph.set_ylabel("Number Of Reads")
    graph.set_title("Coverage")
    graph.set_xlim(0, (df["Pos"].values[-1]+100))
    graph.set_ylim(1, (df["Read_count"].max()+1000))
    plt.axhline(med_reads, color='red', ls='--')
    plt.text(0, med_reads+100, str(med_reads), color="red")
    plt.tight_layout()
    plt.savefig(out_dir + "/coverage_%s.png" % sample, dpi=300)
    plt.close('all')
    return df


def parse_reads(freqs):

    ''' this method returns a vector of reads corresponding to genome positions.
input:
        freqs file
output:
        an integer vector containing for each position in the genome it's num of reads.
'''

    path = freqs
    df = pd.read_table(path, sep="\t")

    # remove all duplicates from Pos except the first occurrence
    # remove all x.number duplicates
    df[["Pos"]] = df[["Pos"]].astype(int)
    df = df.drop_duplicates("Pos")
    df = df[df["Base"] != "-"]
    df = df[df["Read_count"] > 1000]

    pos = df["Pos"]  # a vector of all positions
    reads = df["Read_count"]
    med_reads = reads.median()
    print(med_reads)
    return pos, reads

def main():
#RV
    # sample = "Free_33_Amicon"
    # base_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/capsid"
    # freqs_file_path = base_dir + "/%s/20201012_q38/%s.freqs" % (sample, sample.replace("_", "-"))
    # out_dir = base_dir + "/Coverage_plots"
    # coverage_graph(freqs_file_path, out_dir)
#CV
    sample = "CVB3_RNA_Control"
    base_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/190627_RV_CV/merged/CVB3/"
    freqs_file_path = base_dir + "/%s/q38_3UTR/%s.freqs" % (sample, sample.replace("_", "-"))
    out_dir = base_dir + "/Coverage_plots"
    coverage_graph(freqs_file_path, out_dir)


if __name__ == "__main__":
    main()







