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
    output_dir = input_dir + "/20201112_plots"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)


    data_capsid = pd.read_pickle(input_dir + "capsid/Rank0_data_mutation/q38_data_mutation.pkl")
    data_capsid["RNA"] = np.where((data_capsid["label"] == "RNA Control\nPrimer ID"), "Control",
                                    data_capsid["RNA"])
    data_capsid["label"] = np.where((data_capsid["label"] == "RNA Control\nPrimer ID"), "3`RNA Control\nPrimer ID",
                                    data_capsid["label"])
    data_capsid = data_capsid[data_capsid["label"] != "p8 Mixed Population"]
    data_capsid["RNA"] = np.where((data_capsid["RNA"] == "Capsid"), "RV\nCapsid", data_capsid["RNA"])
    data_capsid["RNA"] = np.where((data_capsid["RNA"] == "Free"), "RV\nFree", data_capsid["RNA"])

    data_capsid = data_capsid.rename(columns={"RNA": "Source"})

    data_passages = pd.read_pickle(input_dir + "passages/Rank0_data_mutation/q38_data_mutation.pkl")
    data_passages["replica"] = data_passages["replica"].astype(int)
    data_passages = data_passages[data_passages["replica"] != 3]
    data_passages = data_passages[data_passages["label"] != "RNA Control\nPrimer ID"]
    data_passages["Source"] = "RV\nPassages"

    data_patients = pd.read_pickle(input_dir + "patients/Rank0_data_mutation/q30_data_mutation.pkl")
    data_patients["Source"] = "RV\nPatients"
    data_patients["Source"] = np.where((data_patients["label"] == "RNA Control\nPrimer ID"), "Control",
                                       data_patients["Source"])
    data_patients["label"] = np.where((data_patients["label"] == "RNA Control\nPrimer ID"), "5`RNA Control\nPrimer ID",
                                       data_patients["label"])
    data_patients["Source"] = np.where((data_patients["label"] == "p3 Cell Culture\nControl"), "RV\nPassages",
                                       data_patients["Source"])


    columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
               "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts",
              "Consensus>Mutated_codon", "Source", "method", "passage", "replica"]
    data = pd.concat([data_capsid, data_passages, data_patients], sort=False)
    data = pd.DataFrame(data, columns=columns)
    # print(data.to_string())

    pis_data = pis_calac.pis_calc(data, pivot_cols=["label", "Source"])

    print(pis_data.to_string())

    source_order = ["Control", "RV\nCapsid", "RV\nFree", "RV\nPassages", "RV\nPatients"]
    pi_palette = (sns.color_palette()[0], sns.color_palette()[1], sns.color_palette()[2], sns.color_palette()[3],
                  sns.color_palette()[4])
    sns.set_context("paper", font_scale=1)
    g = sns.catplot(x="Source", y="Pi", data=pis_data, order=source_order, hue="Source", hue_order=source_order,
                    size=3.5, palette=pi_palette)
    g.set(yscale="log")
    g.set_ylabels("Nucleotide diversity $\pi$")
    g.set(ylim=(2*10**-4, 2*10**-3))
    plt.tight_layout()
    plt.savefig(output_dir + "/Pi_plot_PrimerID", dpi=300)


if __name__ == "__main__":
    main()
