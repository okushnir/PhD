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
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/capsid"
    output_dir = input_dir + "/20201025_plots"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)


    data_mutations = pd.read_pickle(input_dir + "/q38_data_mutation.pkl")

    columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
               "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts",
              "Consensus>Mutated_codon", "RNA", "method", "replica"]
    data = pd.DataFrame(data_mutations, columns=columns)
    data["RNA"] = np.where(data["label"] == "RNA Control\nPrimer ID", "RNA Control\nPrimer ID",
                                    data["RNA"])
    # data["passage"] = data["passage"].astype(int)

    pis_data = pis_calac.pis_calc(data, pivot_cols=["RNA"])

    print(pis_data.to_string())
    replica_order = ("1", "2", "3")
    rna_order = ["RNA Control\nPrimer ID", "Mix Populationֿ\nControl", "Capsid", "Free"]

    sns.set_context("paper", font_scale=0.8)
    g = sns.pointplot(x="RNA", y="Pi", data=pis_data, join=False, order=rna_order)
    g.set_yscale('log')
    g.set(ylim=(2*10**-5, 4*10**-3))
    # plt.xticks(fontsize=10)
    # plt.yticks(fontsize=10)
    # g.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(output_dir + "/Pi_plot", dpi=300)


if __name__ == "__main__":
    main()
