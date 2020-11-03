
import pandas as pd
import os
import matplotlib.pyplot as plt
from FITS_analysis import fits_new_plotter
import numpy as np
from scipy import stats
from AccuNGS_analysis import add_Protein_to_pd_df
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker
import seaborn as sns
from statannot import add_stat_annotation
sns.set_style("ticks")

# print(plt.style.available)

def weighted_varaint(x, **kws):
    var, count = map(np.asarray, zip(*x))
    return var.sum() / count.sum()


def main():
    # input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/"
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/capsid/"
    prefix = "inosine_predict_context"
    output_dir = input_dir + "20201104_%s" % prefix
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    data_filter = pd.read_pickle(input_dir + prefix + "/data_filter.pkl")
    data_filter_ag = pd.read_pickle(input_dir + prefix + "/data_filter_ag.pkl")
    data_filter_uc = pd.read_pickle(input_dir + prefix +"/data_filter_uc.pkl")

    data_filter_ag_replica1 = data_filter_ag[(data_filter_ag["label"] == "Capsid-31-Amicon") |
                                             (data_filter_ag["label"] == "Free-31-Amicon") |
                                             (data_filter_ag["label"] == "RNA Control\nPrimer ID") |
                                             (data_filter_ag["label"] == "Mix Populationֿ\nControl")]
    data_filter_ag_replica1["RNA"] = np.where((data_filter_ag_replica1["RNA"] == "Capsid"), "Capsid #1",
                                              data_filter_ag_replica1["RNA"])
    data_filter_ag_replica1["RNA"] = np.where((data_filter_ag_replica1["RNA"] == "Free"), "Free #1",
                                              data_filter_ag_replica1["RNA"])

    #Plots
    capsid_order = ["RNA Control\nPrimer ID", "Mix Populationֿ\nControl", "Capsid-31-Amicon", "Capsid-32-Ultra",
                    "Capsid-33-Ultra",
                    "Free-31-Amicon", "Free-33-Amicon", "Free-32-Ultra", "Free-33-Ultra"]  #
    rna_order = ["RNA Control\nPrimer ID", "Mix Populationֿ\nControl", "Capsid", "Free"]
    rna_order_replica1 = ["RNA Control\nPrimer ID", "Mix Populationֿ\nControl", "Capsid #1", "Free #1"]
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U"]
    transition_order = ["A>G", "U>C", "G>A", "C>U"]
    type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]
    type_order_ag = ["Synonymous", "Non-Synonymous"]
    context_order_uc = ["UpU", "UpA", "UpC", "UpG"]
    adar_preference = ["High", "Intermediate", "Low"]

    # g1 = sns.catplot(x="label", y="frac_and_weight", data=data_filter, hue="Mutation", order=capsid_order, palette="tab20",
    #                     kind="point", dodge=False, hue_order=mutation_order, join=True, estimator=weighted_varaint,
    #                  orient="v")
    # g1.set_axis_labels("", "Variant Frequency")
    # g1.set_xticklabels(fontsize=9, rotation=45)
    # g1.set(yscale='log')
    # g1.set(ylim=(10**-7, 10**-3))
    # g1.savefig(output_dir + "/All_Mutations_point_plot", dpi=300)
    # plt.close()

    g2 = sns.catplot(x="label", y="frac_and_weight", data=data_filter, hue="Mutation", order=capsid_order, palette="tab20"
                        ,kind="point", dodge=True, hue_order=transition_order, join=False, estimator=weighted_varaint,
                     orient="v", legend=True)
    g2.set_axis_labels("", "Variant Frequency")
    g2.set(yscale='log')
    g2.set(ylim=(10 ** -6, 10 ** -2))
    # g2.set_yticklabels(fontsize=12)
    g2.set_xticklabels(fontsize=10, rotation=90)
    # plt.show()
    # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/Transition_Mutations_point_plot_RV", dpi=300)
    g2.savefig(output_dir + "/Transition_Mutations_point_plot", dpi=300)
    plt.close()

    g_rna = sns.catplot(x="RNA", y="frac_and_weight", data=data_filter, hue="Mutation", order=rna_order,
                     palette="tab20", kind="point", dodge=True, hue_order=transition_order, join=False, estimator=weighted_varaint,
                     orient="v", legend=True)
    g_rna.set_axis_labels("", "Variant Frequency")
    g_rna.set(yscale='log')
    g_rna.set(ylim=(10 ** -6, 10 ** -2))
    g_rna.savefig(output_dir + "/Transition_Mutations_point_RNA_plot", dpi=300)
    plt.close()

    # A>G Prev Context
    g4 = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="5`_ADAR_Preference", order=capsid_order,
                     palette="rocket", kind="point", dodge=True, hue_order=adar_preference, estimator=weighted_varaint,
                     orient="v", col="Type", join=False, col_order=type_order_ag)
    g4.set_axis_labels("", "Variant Frequency")
    g4.set(yscale='log')
    g4.set(ylim=(7 * 10 ** -7, 4 * 10 ** -3))
    g4.set_xticklabels(rotation=90)
    # plt.show()
    g4.savefig(output_dir + "/Context_label_point_plot", dpi=300)
    plt.close()

    mutation_g8 = sns.catplot("RNA", "frac_and_weight", data=data_filter_ag_replica1, hue="5`_ADAR_Preference",
                     palette="rocket", kind="point", dodge=True, estimator=weighted_varaint, order=rna_order_replica1,
                     orient="v", col="Type", join=False, col_order=type_order_ag, hue_order=adar_preference)
    mutation_g8.set(yscale="log")
    # mutation_g8.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    mutation_g8.set_axis_labels("", "Variant Frequency")
    mutation_g8.set(ylim=(7 * 10 ** -7, 4 * 10 ** -3))
    # plt.show()
    mutation_g8.savefig(output_dir + "/ADAR_like_AG_Mutation_col.png", dpi=300)
    plt.close()

    data_filter_ag_grouped = data_filter_ag.groupby(["5`_ADAR_Preference", "label", "Type", "Pos", "Protein", "RNA"])["frac_and_weight"].agg(lambda x: weighted_varaint(x))
    data_filter_ag_grouped = data_filter_ag_grouped.reset_index()
    data_filter_ag_grouped = data_filter_ag_grouped.rename(columns={"frac_and_weight": "Frequency"})
    data_filter_ag_grouped["Frequency"] = data_filter_ag_grouped["Frequency"].astype(float)
    # data_filter_ag_grouped = data_filter_ag_grouped[data_filter_ag_grouped["RNA"] = "Capsid"]
    print(data_filter_ag_grouped.to_string())

    data_filter_ag_grouped_silent = data_filter_ag_grouped[data_filter_ag_grouped["Type"] == "Synonymous"]
    data_filter_ag_grouped_silent = data_filter_ag_grouped_silent[data_filter_ag_grouped_silent["Protein"] != "2A"]
    data_filter_ag_grouped_silent = data_filter_ag_grouped_silent[data_filter_ag_grouped_silent["Protein"] != "3'UTR"]

    position_mutation = sns.relplot(x="Pos", y="Frequency", data=data_filter_ag_grouped_silent, hue="5`_ADAR_Preference",
                                    col="RNA", col_wrap=2, style="5`_ADAR_Preference", palette="rocket",
                                    hue_order=adar_preference, style_order=["High", "Low", "Intermediate"], height=4)

    position_mutation.set_axis_labels("", "Variant Frequency")
    position_mutation.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -4)
    position_mutation.axes.flat[0].set_ylim(10**-5, 10 ** -2)
    plt.savefig(output_dir + "/position_mutation.png", dpi=300)
    plt.close()

    # g5 = sns.catplot("RNA", "frac_and_weight", data=data_filter_ag_replica1, hue="ADAR_like", order=rna_order_replica1,
    #                  palette=flatui, kind="point", dodge=True, hue_order=[True, False], estimator=weighted_varaint,
    #                  orient="v", col="Type", join=False, col_order=type_order_ag)
    # g5.set_axis_labels("", "Variant Frequency")
    # g5.set(yscale='log')
    # g5.set(ylim=(7 * 10 ** -7, 4 * 10 ** -3))
    # # g5.set_xticklabels("RNA Control\nPrimer ID", "Mix Populationֿ\nControl", "Capsid #1", "Free #1")
    # # plt.show()
    # g5.savefig(output_dir + "/Context_RNA_point_plot", dpi=300)
    # plt.close()
    #
    # g6 = sns.catplot("RNA", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=rna_order, palette=flatui,
    #                  kind="point", dodge=True, hue_order=[True, False], estimator=weighted_varaint, orient="v",
    #                  join=False)
    # g6.set_axis_labels("", "Variant Frequency")
    # g6.set(yscale='log')
    # g6.set(ylim=(7 * 10 ** -7, 4 * 10 ** -3))
    # # g6.set_xticklabels(rotation=45)
    # # plt.show()
    # g6.savefig(output_dir + "/Context_point_all_mutations_type_plot", dpi=300)
    # plt.close()


# g9 = sns.catplot("RNA", "frac_and_weight", data=data_filter_uc, hue="ADAR_like", order=rna_order, palette=flatui,
    #                     estimator=weighted_varaint, orient="v", dodge=True, kind="point",
    #                  col="Type", join=False, col_order=type_order_ag)
    # g9.set_axis_labels("", "Variant Frequency")
    # g9.set(yscale='log')
    # g9.set(ylim=(10 ** -7, 10 ** -3))
    # g9.set_xticklabels(rotation=45)
    # # plt.show()
    # g9.savefig(output_dir + "/UC_Context_point_plot", dpi=300)
    # plt.close()

if __name__ == "__main__":
    main()

