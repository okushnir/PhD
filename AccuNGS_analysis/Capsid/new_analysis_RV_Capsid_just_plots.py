
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
from datetime import datetime
from statannotations.Annotator import Annotator
import contextlib

sns.set_style("ticks")

# print(plt.style.available)

def weighted_varaint(x, **kws):
    var, count = map(np.asarray, zip(*x))
    return var.sum() / count.sum()


def main():
    # input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/"
    # input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/capsid/"
    # input_dir = "/Users/odedkushnir/PhD_Projects/After_review/AccuNGS/RV/capsid/"
    input_dir = "C:/Users/odedku/PhD_Projects/After_review/AccuNGS/RV/capsid/"
    prefix = "inosine_predict_context"
    date = datetime.today().strftime("%Y%m%d")
    output_dir = input_dir + "{0}_{1}".format(date, prefix)
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    data_filter = pd.read_pickle(input_dir + prefix + "/data_filter.pkl")
    data_filter_ag = pd.read_pickle(input_dir + prefix + "/data_filter_ag.pkl")
    data_filter_uc = pd.read_pickle(input_dir + prefix + "/data_filter_uc.pkl")

    data_filter_replica1 = data_filter[(data_filter["label"] == "Capsid-31-Amicon") |
                                             (data_filter["label"] == "Free-31-Amicon") |
                                             (data_filter["label"] == "RNA Control\nPrimer ID") |
                                             (data_filter["label"] == "p8 Mixed Population")]
    data_filter_replica1["RNA"] = np.where((data_filter_replica1["RNA"] == "Capsid"), "p9 Capsid #1",
                                              data_filter_replica1["RNA"])
    data_filter_replica1["RNA"] = np.where((data_filter_replica1["RNA"] == "Free"), "p9 Free #1",
                                              data_filter_replica1["RNA"])
    data_filter_replica1["RNA"] = np.where((data_filter_replica1["RNA"] == "p8 Mixed Population"),
                                           "p8 Mixed\nPopulation", data_filter_replica1["RNA"])


    data_filter_replica1.to_csv(input_dir + prefix + "/data_filter_rep1.csv", sep=",", encoding='utf-8')

    data_filter_ag_replica1 = data_filter_ag[(data_filter_ag["label"] == "Capsid-31-Amicon") |
                                             (data_filter_ag["label"] == "Free-31-Amicon") |
                                             (data_filter_ag["label"] == "RNA Control\nPrimer ID") |
                                             (data_filter_ag["label"] == "p8 Mixed Population")]
    data_filter_ag_replica1["RNA"] = np.where((data_filter_ag_replica1["RNA"] == "Capsid"), "p9 Capsid #1",
                                              data_filter_ag_replica1["RNA"])
    data_filter_ag_replica1["RNA"] = np.where((data_filter_ag_replica1["RNA"] == "Free"), "p9 Free #1",
                                              data_filter_ag_replica1["RNA"])
    data_filter_ag_replica1["RNA"] = np.where((data_filter_ag_replica1["RNA"] == "RNA Control\nPrimer ID"),
                                           "RNA\nControl", data_filter_ag_replica1["RNA"])

    data_filter_ag_replica1["RNA"] = np.where((data_filter_ag_replica1["RNA"] == "p8 Mixed Population"),
                                           "p8 Mixed\nPopulation", data_filter_ag_replica1["RNA"])
    data_filter_ag_replica1.to_csv(input_dir + prefix +"/data_filter_ag_rep1.csv", sep=",", encoding='utf-8')
    #Plots
    label_order = ["RNA\n Control", "p8 Mixed Population", "Capsid-31-Amicon", "Capsid-32-Ultra",
                   "Capsid-33-Ultra", "Free-31-Amicon", "Free-32-Ultra", "Free-33-Ultra"]  #

    rna_order_replica1 = ["RNA\nControl", "p8 Mixed\nPopulation", "p9 Capsid #1", "p9 Free #1"]
    transition_order = ["A>G", "U>C", "G>A", "C>U"]
    type_order_ag = ["Synonymous", "Non-Synonymous"]
    adar_preference = ["High", "Intermediate", "Low"]
    rna_order = ["RNA\n Control", "p8 Mixed Population", "Capsid", "Free"]
    type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U"]
    context_order_uc = ["UpU", "UpA", "UpC", "UpG"]
    context_order = ["UpA", "ApA", "CpA", "GpA"]
    type_order = ["Synonymous", "Non-Synonymous"]
    plus_minus = u"\u00B1"
    pairs = [(("RNA\nControl", "High"), ("RNA\nControl", "Intermediate")),
             (("RNA\nControl", "High"), ("RNA\nControl", "Low")),
             (("p8 Mixed\nPopulation", "High"), ("p8 Mixed\nPopulation", "Intermediate")),
             (("p8 Mixed\nPopulation", "High"), ("p8 Mixed\nPopulation", "Low")),
             (("p9 Capsid #1", "High"), ("p9 Capsid #1", "Intermediate")),
             (("p9 Capsid #1", "High"), ("p9 Capsid #1", "Low")),
             (("p9 Free #1", "High"), ("p9 Free #1", "Intermediate")),
             (("p9 Free #1", "High"), ("p9 Free #1", "Low"))]

    # g1 = sns.catplot(x="label", y="frac_and_weight", data=data_filter, hue="Mutation", order=label_order, palette="tab20",
    #                     kind="point", dodge=False, hue_order=mutation_order, join=True, estimator=weighted_varaint,
    #                  orient="v")
    # g1.set_axis_labels("", "Variant Frequency")
    # g1.set_xticklabels(fontsize=9, rotation=45)
    # g1.set(yscale='log')
    # g1.set(ylim=(10**-7, 10**-3))
    # g1.savefig(output_dir + "/All_Mutations_point_plot", dpi=300)
    # plt.close()

    g2 = sns.catplot(x="label", y="frac_and_weight", data=data_filter, hue="Mutation", order=label_order,
                     palette=mutation_palette(4), kind="point", dodge=True, hue_order=transition_order, join=False,
                     estimator=weighted_varaint, orient="v", legend=True)
    g2.set_axis_labels("", "Variant Frequency {} CI=95%".format(plus_minus))
    g2.set(yscale='log')
    g2.set(ylim=(10 ** -6, 10 ** -2))
    # g2.set_yticklabels(fontsize=12)
    g2.set_xticklabels(fontsize=10, rotation=90)
    # plt.show()
    # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/Transition_Mutations_point_plot_RV", dpi=300)
    g2.savefig(output_dir + "/Transition_Mutations_label_point_plot", dpi=300)
    plt.close()

    g_rna = sns.catplot(x="RNA", y="frac_and_weight", data=data_filter_replica1, hue="Mutation", order=rna_order_replica1,
                     palette=mutation_palette(4), kind="point", dodge=True, hue_order=transition_order, join=False,
                        estimator=weighted_varaint, orient="v", legend=True)
    g_rna.set_axis_labels("", "Variant Frequency {} CI=95%".format(plus_minus))
    g_rna.set(yscale='log')
    g_rna.set(ylim=(10 ** -6, 10 ** -2))
    g_rna.savefig(output_dir + "/Transition_Mutations_RNA_point_plot", dpi=300)
    plt.close()

    # A>G Prev Context
    g4 = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="5`_ADAR_Preference", order=label_order,
                     palette=mutation_palette(3, adar=True, ag=True), kind="point", dodge=True,
                     hue_order=adar_preference, estimator=weighted_varaint,
                     orient="v", col="Type", join=False, col_order=type_order_ag)
    g4.set_axis_labels("", "Variant Frequency {} CI=95%".format(plus_minus))
    g4.set(yscale='log')
    g4.set(ylim=(7 * 10 ** -7, 4 * 10 ** -3))
    g4.set_xticklabels(rotation=90)
    # plt.show()
    g4.savefig(output_dir + "/Context_label_point_plot", dpi=300)
    plt.close()

    mutation_g8 = sns.catplot("RNA", "frac_and_weight", data=data_filter_ag_replica1, hue="5`_ADAR_Preference",
                     palette=mutation_palette(3, adar=True, ag=True), kind="point", dodge=True,
                              estimator=weighted_varaint, order=rna_order_replica1,
                              orient="v", col="Type", join=False, col_order=type_order_ag, hue_order=adar_preference)
    mutation_g8.set(yscale="log")
    # mutation_g8.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    mutation_g8.set_axis_labels("", "Variant Frequency {} CI=95%".format(plus_minus))
    mutation_g8.set(ylim=(1 * 10 ** -5, 1 * 10 ** -2))
    # plt.show()
    mutation_g8.savefig(output_dir + "/ag_ADAR_like_Mutation_col.png", dpi=300)
    plt.close()

    mutation_type_g1 = sns.boxplot(x="RNA", y="Frequency", data=data_filter_ag_replica1, hue="5`_ADAR_Preference",
                                   order=rna_order_replica1, palette=mutation_palette(3, adar=True, ag=True), dodge=True,
                                   hue_order=adar_preference)
    mutation_type_g1.set_yscale('log')
    mutation_type_g1.set_ylim(10 ** -5, 10 ** -2)
    mutation_type_g1.set(xlabel="", ylabel="Variant Frequency")
    annot = Annotator(mutation_type_g1, pairs, x="RNA", y="Frequency", hue="5`_ADAR_Preference",
                      data=data_filter_ag_replica1, order=rna_order_replica1, hue_order=adar_preference)
    annot.configure(test='t-test_welch', text_format='star', loc='outside', verbose=2,
                    comparisons_correction="Bonferroni")
    annot.apply_test()
    file_path = output_dir + "/sts_adar.csv"
    with open(file_path, "w") as o:
        with contextlib.redirect_stdout(o):
            passage_g1, test_results = annot.annotate()
    plt.legend(bbox_to_anchor=(1.05, 0.5), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(output_dir + "/ag_ADAR_like_Mutation_box_sts.png", dpi=300)
    plt.close()

    mutation_g8_box = sns.catplot("RNA", "Frequency", data=data_filter_ag_replica1, hue="5`_ADAR_Preference",
                                  palette=mutation_palette(3, adar=True, ag=True), order=rna_order_replica1,
                     col="Type",col_order=type_order_ag, hue_order=adar_preference, kind="box")
    mutation_g8_box.set(yscale="log")
    mutation_g8_box.savefig(output_dir + "/ag_ADAR_like_Mutation_col_box.png", dpi=300)

    data_filter_ag_grouped = data_filter_ag.groupby(["5`_ADAR_Preference", "label", "Type", "RNA", "Pos"])["frac_and_weight"].agg(lambda x: weighted_varaint(x))
    data_filter_ag_grouped = data_filter_ag_grouped.reset_index()
    data_filter_ag_grouped = data_filter_ag_grouped.rename(columns={"frac_and_weight": "Frequency"})
    data_filter_ag_grouped["Frequency"] = data_filter_ag_grouped["Frequency"].astype(float)
    # data_filter_ag_grouped = data_filter_ag_grouped[data_filter_ag_grouped["RNA"] = "Capsid"]
    print(data_filter_ag_grouped.to_string())

    data_filter_ag_grouped_silent = data_filter_ag_grouped[data_filter_ag_grouped["Type"] == "Synonymous"]
    # data_filter_ag_grouped_silent = data_filter_ag_grouped_silent[data_filter_ag_grouped_silent["Protein"] != "2A"]
    # data_filter_ag_grouped_silent = data_filter_ag_grouped_silent[data_filter_ag_grouped_silent["Protein"] != "3'UTR"]

    position_mutation = sns.relplot(x="Pos", y="Frequency", data=data_filter_ag_grouped_silent, hue="5`_ADAR_Preference",
                                    col="RNA", col_wrap=2, style="5`_ADAR_Preference", palette=mutation_palette(3, adar=True, ag=True),
                                    hue_order=adar_preference, style_order=["High", "Low", "Intermediate"], height=4)

    position_mutation.set_axis_labels("", "Variant Frequency")
    position_mutation.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -4)
    position_mutation.axes.flat[0].set_ylim(10**-4, 10 ** -1)
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

