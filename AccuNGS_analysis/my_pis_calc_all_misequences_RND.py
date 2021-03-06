import pandas as pd
import os
import matplotlib.pyplot as plt
from FITS_analysis import fits_new_plotter
import numpy as np
from scipy import stats
import pis_calac
import seaborn as sns

sns.set_style("ticks")

def main():
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/"
    output_dir = input_dir + "/20201120_plots"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)




    data_passages = pd.read_pickle(input_dir + "passages/Rank0_data_mutation/q38_data_mutation.pkl")
    data_passages["replica"] = data_passages["replica"].astype(int)
    data_passages = data_passages[data_passages["replica"] == 3]
    data_passages["Source"] = "RV\nPassages"
    # data_passages["Source"] = np.where((data_passages["label"] == "RNA Control_RND"), "Control", data_passages["Source"])
    data_passages["Source"] = np.where((data_passages["passage"] == 0), "Control",
                                       data_passages["Source"])

    data_cvb3 = pd.read_pickle("/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/CVB3/Rank0_data_mutation/"
                               "q38_data_mutation.pkl")
    data_cvb3["Source"] = "CV\nPassages"
    data_cvb3["Source"] = np.where((data_cvb3["label"] == "CVB3-RNA Control"), "Control", data_cvb3["Source"])


    columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
               "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts",
              "Consensus>Mutated_codon", "Source", "method", "passage", "replica"]
    data = pd.concat([data_passages, data_cvb3], sort=False)
    data = pd.DataFrame(data, columns=columns)

    pis_data = pis_calac.pis_calc(data, pivot_cols=["label", "Source"])

    print(pis_data.to_string())

    source_order = ["Control",  "RV\nPassages", "CV\nPassages"]
    pi_palette = (sns.color_palette()[0], sns.color_palette()[3], sns.color_palette()[6])
    sns.set_context("paper", font_scale=1.2)
    g = sns.catplot(x="Source", y="Pi", data=pis_data, order=source_order, hue="Source", hue_order=source_order,
                    size=3.5, palette=pi_palette, height=4)
    g.set(yscale="log")
    g.set_ylabels("Nucleotide diversity $\pi$")
    g.set_xlabels("")
    g.set(ylim=(2*10**-4, 2*10**-3))
    plt.tight_layout()
    plt.savefig(output_dir + "/Pi_plot_rnd_primer", dpi=300)


if __name__ == "__main__":
    main()
