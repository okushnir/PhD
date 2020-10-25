
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
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/capsid"
    output_dir = input_dir + "/20201025_plots"
    data_filter = pd.read_pickle(output_dir + "/data_filter.pkl")
    data_filter_ag = pd.read_pickle(output_dir + "/data_filter_ag.pkl")
    data_filter_uc = pd.read_pickle(output_dir + "/data_filter_uc.pkl")


    capsid_order = ["RNA Control\nPrimer ID", "Mix Populationֿ\nControl","Capsid-31-Amicon", "Capsid-33-Ultra",
                    "Free-31-Amicon", "Free-33-Amicon", "Free-32-Ultra", "Free-33-Ultra"]#, "Capsid-32-Ultra"
    rna_order = ["RNA Control\nPrimer ID", "Mix Populationֿ\nControl", "Capsid", "Free"]
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U"]
    transition_order = ["A>G", "U>C", "G>A", "C>U"]
    type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]
    context_order = ["UpA", "ApA", "CpA", "GpA"]
    type_order = ["Synonymous", "Non-Synonymous"]
    # data_filter_ag["weighted_varaint"] = data_filter_ag["frac_and_weight"].apply(lambda x: weighted_varaint(x))
    #

    #Plots
    g1 = sns.catplot(x="label", y="frac_and_weight", data=data_filter, hue="Mutation", order=capsid_order, palette="tab20",
                        kind="point", dodge=False, hue_order=mutation_order, join=True, estimator=weighted_varaint,
                     orient="v")
    g1.set_axis_labels("", "Variant Frequency")
    g1.set_xticklabels(fontsize=9, rotation=45)
    g1.set(yscale='log')
    g1.set(ylim=(10**-7, 10**-3))

    # plt.show()
    g1.savefig(output_dir + "/All_Mutations_point_plot", dpi=300)
    plt.close()
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
    # g2.set_yticklabels(fontsize=12)
    g_rna.set_xticklabels(fontsize=10, rotation=45)
    g_rna.savefig(output_dir + "/Transition_Mutations_point_RNA_plot", dpi=300)
    plt.close()

    # A>G Prev Context
    flatui = ["#3498db", "#9b59b6"]
    g5 = sns.catplot("RNA", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=rna_order, palette=flatui,
                kind="point", dodge=True, hue_order=[True, False], estimator=weighted_varaint, orient="v",
                     col="Type", join=False, col_order=type_order)
    g5.set_axis_labels("", "Variant Frequency")
    g5.set(yscale='log')
    g5.set(ylim=(7*10**-7, 4*10**-3))
    g5.set_xticklabels(rotation=45)
    g5.savefig(output_dir + "/Context_RNA_point_plot", dpi=300)
    plt.close()

    # data_filter_ag_synon = data_filter_ag[data_filter_ag["Type"] == "Synonymous"]
    # g5_synon = sns.pointplot("RNA", "frac_and_weight", data=data_filter_ag_synon, hue="ADAR_like", order=rna_order,
    #                        palette=flatui, kind="point", dodge=True, hue_order=[True, False],
    #                        estimator=weighted_varaint, orient="v", col="Type", join=False, col_order=type_order)
    # add_stat_annotation(g5_synon, data=data_filter_ag_synon, x="RNA", y="Frequency", hue="ADAR_like", order=rna_order,
    #                     box_pairs=[(("Capsid", True), ("Free", True))],
    #                     test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
    # g5_synon.set_yscale('log')
    # g5_synon.set_ylabel("Variant Frequency")
    # g5_synon.set(ylim=(7*10**-7, 4*10**-3))
    # sns.despine()
    # plt.tight_layout()
    # plt.savefig(output_dir + "/Context_RNA_synon_point_plot", dpi=300)
    # plt.close()

    g6 = sns.catplot("RNA", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=rna_order, palette=flatui,
                kind="point", dodge=True, hue_order=[True, False], estimator=weighted_varaint, orient="v",
                     join=False)
    g6.set_axis_labels("", "Variant Frequency")
    g6.set(yscale='log')
    g6.set(ylim=(7*10**-7, 4*10**-3))
    g6.set_xticklabels(rotation=45)
    # plt.show()
    g6.savefig(output_dir + "/Context_point_all_mutations_type_plot", dpi=300)
    plt.close()

    # g9 = sns.catplot("RNA", "frac_and_weight", data=data_filter_uc, hue="Next", order=rna_order, palette="tab20",
    #                  hue_order=context_order_uc, estimator=weighted_varaint, orient="v", dodge=True, kind="point",
    #                  col="Type", join=False, col_order=type_order)
    # g9.set_axis_labels("", "Variant Frequency")
    # g9.set(yscale='log')
    # g9.set(ylim=(10 ** -7, 10 ** -3))
    # g9.set_xticklabels(rotation=45)
    # # plt.show()
    # g9.savefig(output_dir + "/UC_Context_point_plot", dpi=300)
    # plt.close()

    data_filter_ag_grouped = data_filter_ag.groupby(["ADAR_like", "label", "Type", "Pos", "Protein", "RNA"])["frac_and_weight"].agg(lambda x: weighted_varaint(x))
    data_filter_ag_grouped = data_filter_ag_grouped.reset_index()
    data_filter_ag_grouped = data_filter_ag_grouped.rename(columns={"frac_and_weight": "Frequency"})
    data_filter_ag_grouped["Frequency"] = data_filter_ag_grouped["Frequency"].astype(float)
    data_filter_ag_grouped = data_filter_ag_grouped[data_filter_ag_grouped["RNA"] == "Capsid"]
    print(data_filter_ag_grouped.to_string())

    data_filter_ag_grouped_silent = data_filter_ag_grouped[data_filter_ag_grouped["Type"] == "Synonymous"]
    data_filter_ag_grouped_silent = data_filter_ag_grouped_silent[data_filter_ag_grouped_silent["Protein"] != "2A"]
    position_g = sns.scatterplot("Pos", "Frequency", data=data_filter_ag_grouped_silent, hue="Protein",
                                 palette="tab10", style="ADAR_like", style_order=[True, False], legend="full")

    # position_g.set_axis_labels("", "Variant Frequency")
    position_g.set_yscale('log')
    position_g.set_ylim(10 ** -4, 10 ** -1)
    position_g.set(xlim=(3500, 7500))
    position_g.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(output_dir + "/position.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    main()

