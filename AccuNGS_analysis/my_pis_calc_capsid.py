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
    output_dir = input_dir + "/20201028_plots"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)


    data_capsid = pd.read_pickle(input_dir + "capsid/Rank0_data_mutation/q38_data_mutation.pkl")
    data_capsid["RNA"] = np.where((data_capsid["label"] == "RNA Control\nPrimer ID"), "Control",
                                    data_capsid["RNA"])
    data_capsid = data_capsid[data_capsid["label"] != "Mix Populationֿ\nControl"]

    data_passages = pd.read_pickle(input_dir + "passages/Rank0_data_mutation/q38_data_mutation.pkl")
    data_passages = data_passages[data_passages["label"] != "RNA Control\nPrimer ID"]
    data_passages["RNA"] = "Passages"
    data_passages["RNA"] = np.where((data_passages["label"] == "RNA Control_RND"), "Control", data_passages["RNA"])

    data_patients = pd.read_pickle(input_dir + "patients/Rank0_data_mutation/q30_data_mutation.pkl")
    data_patients = data_patients[data_patients["label"] != "RNA Control\nPrimer ID"]
    data_patients = data_patients[data_patients["label"] != "Cell Cultureֿ\nControl"]
    data_patients["RNA"] = "Patients"

    columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
               "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts",
              "Consensus>Mutated_codon", "RNA", "method", "passage", "replica"]
    data = pd.concat([data_capsid, data_passages, data_patients], sort=False)
    data = pd.DataFrame(data, columns=columns)

    # data["RNA"] = np.where(data["label"] == "p2-1", "P",
    #                        data["RNA"])
    # data["passage"] = data["passage"].astype(int)
    data = data.rename(columns={"RNA": "Source"})

    pis_data = pis_calac.pis_calc(data, pivot_cols=["label", "Source"])

    print(pis_data.to_string())

    source_order = ["Control", "Capsid", "Free", "Passages", "Patients"]

    sns.set_context("paper", font_scale=1.2)
    g = sns.catplot(x="Source", y="Pi", data=pis_data, order=source_order, hue="Source", hue_order=source_order)
    g.set(yscale="log")
    # g.set_yscale('log')
    g.set(ylim=(2*10**-4, 2*10**-3))
    # plt.tight_layout()
    plt.savefig(output_dir + "/Pi_plot", dpi=300)


if __name__ == "__main__":
    main()
