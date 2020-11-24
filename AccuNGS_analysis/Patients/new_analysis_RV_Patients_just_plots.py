
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
from AccuNGS_analysis.adar_mutation_palette import mutation_palette

sns.set(font_scale=1.2)
sns.set_style("ticks")
sns.despine()

# print(plt.style.available)

def weighted_varaint(x, **kws):
    var, count = map(np.asarray, zip(*x))
    return var.sum() / count.sum()


def main():
    # input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/"
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/patients/"
    prefix = "inosine_predict_context_freq0.01"
    output_dir = input_dir + "20201124_%s" % prefix
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    data_filter = pd.read_pickle(input_dir + prefix + "/data_filter.pkl")
    data_filter_ag = pd.read_pickle(input_dir + prefix + "/data_filter_ag.pkl")
    data_filter_uc = pd.read_pickle(input_dir + prefix +"/data_filter_uc.pkl")

    #Plots
    label_order = ["RNA Control\nPrimer ID", "p3 Cell Culture\nControl", "Patient-1", "Patient-4", "Patient-5",
                   "Patient-9", "Patient-16", "Patient-17", "Patient-20"]
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U"]
    transition_order = ["A>G", "U>C", "G>A", "C>U"]
    type_order1 = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]
    context_order = ["UpA", "ApA", "CpA", "GpA"]
    type_order2 = ["Synonymous", "Non-Synonymous"]
    context_order_uc = ["UpA", "UpU", "UpG",  "UpC"]
    type_order_ag = ["Synonymous", "Non-Synonymous", "NonCodingRegion"]
    adar_preference = ["High", "Intermediate", "Low"]

    g1 = sns.catplot(x="label", y="frac_and_weight", data=data_filter, hue="Mutation", order=label_order, palette="tab20",
                        kind="point", dodge=True, hue_order=mutation_order, join=False, estimator=weighted_varaint,
                     orient="v")
    g1.set_axis_labels("", "Variant Frequency")
    g1.set_xticklabels(fontsize=9, rotation=90)
    g1.set(yscale='log')
    # g1.set(ylim=(10**-7, 10**-3))

    # plt.show()
    g1.savefig(output_dir + "/All_Mutations_point_plot", dpi=300)
    plt.close()
    g2 = sns.catplot(x="label", y="frac_and_weight", data=data_filter, hue="Mutation", order=label_order, palette=mutation_palette(4)
                        ,kind="point", dodge=True, hue_order=transition_order, join=False, estimator=weighted_varaint,
                     orient="v", legend=True)
    g2.set_axis_labels("", "Variant Frequency")
    g2.set(yscale='log')
    g2.set(ylim=(10 ** -5, 10 ** -3))
    # g2.set_yticklabels(fontsize=12)
    g2.set_xticklabels(fontsize=10, rotation=90)
    # plt.show()
    # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/Transition_Mutations_point_plot_RV", dpi=300)
    g2.savefig(output_dir + "/Transition_Mutations_point_plot", dpi=300)
    # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/Prgress reports/20200913 Final report/plots" +
    #                   "/Fig9a_Transition_Mutations_point_plot_Patients", dpi=300)
    plt.close()

    # g_rna = sns.catplot(x="RNA", y="frac_and_weight", data=data_filter, hue="Mutation", order=rna_order,
    #                  palette="tab20", kind="point", dodge=True, hue_order=transition_order, join=False, estimator=weighted_varaint,
    #                  orient="v", legend=True)
    # g_rna.set_axis_labels("", "Variant Frequency")
    # g_rna.set(yscale='log')
    # g_rna.set(ylim=(10 ** -6, 10 ** -2))
    # # g2.set_yticklabels(fontsize=12)
    # g_rna.set_xticklabels(fontsize=10, rotation=45)
    # plt.show()
    # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/Transition_Mutations_point_plot_RV", dpi=300)
    # g_rna.savefig(output_dir + "/Transition_Mutations_point_RNA_plot", dpi=300)
    # plt.close()

    # A>G Prev Context
    flatui = ["#3498db", "#9b59b6"]
    g5 = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=label_order,
                     palette=mutation_palette(2), kind="point", dodge=True, hue_order=[True, False],
                     estimator=weighted_varaint, orient="v", col="Type", join=False, col_order=type_order2)
    g5.set_axis_labels("", "Variant Frequency")
    g5.set(yscale='log')
    g5.set(ylim=(7*10**-7, 4*10**-3))
    g5.set_xticklabels(rotation=90)
    # plt.show()
    g5.savefig(output_dir + "/Context_point_plot", dpi=300)
    # g5.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/Prgress reports/20200913 Final report/plots" +
    #            "/Fig9b_Context_point_plot_Patients", dpi=300)
    plt.close()

    mutation_ag = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="5`_ADAR_Preference",
                              palette=mutation_palette(3, adar=True, ag=True), kind="point", dodge=True,
                              estimator=weighted_varaint, order=label_order,
                              orient="v", col="Type", join=False, col_order=type_order_ag, hue_order=adar_preference)
    mutation_ag.set(yscale="log")
    mutation_ag.set(ylim=(1 * 10 ** -5, 1 * 10 ** -2))
    mutation_ag.set_xticklabels(rotation=90)
    mutation_ag.fig.suptitle("A>G ADAR_like Mutation in RV patients", y=0.99)
    plt.subplots_adjust(top=0.85)
    mutation_ag.set_axis_labels("", "Variant Frequency")
    mutation_ag.savefig(output_dir + "/ag_ADAR_like_Mutation_col_patients.png", dpi=300)
    plt.close()

    g6 = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=label_order,
                     palette=mutation_palette(2), kind="point", dodge=True, hue_order=[True, False],
                     estimator=weighted_varaint, orient="v", join=False)
    g6.set_axis_labels("", "Variant Frequency")
    g6.set(yscale='log')
    g6.set(ylim=(7*10**-7, 4*10**-3))
    g6.set_xticklabels(rotation=90)
    # plt.show()
    g6.savefig(output_dir + "/Context_point_all_mutations_type_plot", dpi=300)
    plt.close()

    g9 = sns.catplot("label", "frac_and_weight", data=data_filter_uc, hue="Next", order=label_order, palette="tab20",
                     hue_order=context_order_uc, estimator=weighted_varaint, orient="v", dodge=True, kind="point",
                     col="Type", join=False, col_order=type_order2)
    g9.set_axis_labels("", "Variant Frequency")
    g9.set(yscale='log')
    g9.set(ylim=(10 ** -5, 10 ** -2))
    g9.set_xticklabels(rotation=90)
    # plt.show()
    g9.savefig(output_dir + "/UC_Context_point_plot", dpi=300)
    plt.close()

    data_filter_ag_grouped = data_filter_ag.groupby(["ADAR_like", "label", "Type"])["frac_and_weight"].agg(lambda x: weighted_varaint(x))
    data_filter_ag_grouped = data_filter_ag_grouped.reset_index()
    data_filter_ag_grouped = data_filter_ag_grouped.rename(columns={"frac_and_weight": "Frequency"})
    data_filter_ag_grouped["Frequency"] = data_filter_ag_grouped["Frequency"].astype(float)
    print(data_filter_ag_grouped.to_string())

    data_filter_ag_grouped_silent = data_filter_ag_grouped[data_filter_ag_grouped["Type"] == "Synonymous"]
    data_filter_ag_grouped_silent = data_filter_ag_grouped_silent[data_filter_ag_grouped_silent["label"] == "Cell Cultureֿ\nControl"]

    # position_mutation = sns.relplot(x="Pos", y="Frequency", data=data_filter_ag_grouped_silent, hue="5`_ADAR_Preference",
    #                                 col="RNA", col_wrap=2, style="5`_ADAR_Preference", palette=mutation_palette(3, adar=True, ag=True),
    #                                 hue_order=adar_preference, style_order=["High", "Low", "Intermediate"], height=4)
    #
    # position_mutation.set_axis_labels("", "Variant Frequency")
    # position_mutation.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -4)
    # position_mutation.axes.flat[0].set_ylim(10**-5, 10 ** -2)
    # plt.savefig(output_dir + "/position_mutation.png", dpi=300)
    # plt.close()

if __name__ == "__main__":
    main()

