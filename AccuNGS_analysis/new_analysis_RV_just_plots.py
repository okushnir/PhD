
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

# tips = sns.load_dataset("tips")
# print(tips.to_string())
# tips["weight"] = 10 * np.random.rand(len(tips))
#
# tips["tip_and_weight"] = zip(tips.tip, tips.weight)

def weighted_varaint(x, **kws):
    var, count = map(np.asarray, zip(*x))
    return var.sum() / count.sum()


def main():
    flatui = ["#3498db", "#9b59b6"]
    # input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/"
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/new"
    output_dir = input_dir + "/20201019_new_plots"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    data_filter = pd.read_pickle(output_dir + "/data_filter.pkl")
    data_filter_ag = pd.read_pickle(output_dir + "/data_filter_ag.pkl")
    data_filter_uc = pd.read_pickle(output_dir + "/data_filter_uc.pkl")
    label_order = ["RNA Control\nRND", "RNA Control\nPrimer ID","p2-1", "p2-2", "p2-3", "p5-1", "p5-2", "p5-3", "p8-1",
                   "p8-2", "p8-3", "p10-2", "p10-3", "p12-1", "p12-2", "p12-3"]
    passage_order = ["0", "2", "5", "8", "10", "12"]
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U"]
    transition_order = ["A>G", "U>C", "G>A", "C>U"]
    # type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]

    #Plots
    # g1 = sns.catplot(x="label", y="frac_and_weight", data=data_filter, hue="Mutation", order=label_order,
    #                  palette="tab20",
    #                  kind="point", dodge=False, hue_order=mutation_order, join=True, estimator=weighted_varaint,
    #                  orient="v")
    # g1.set_axis_labels("", "Variant Frequency")
    # g1.set_xticklabels(fontsize=9, rotation=45)
    # g1.set(yscale='log')
    # g1.set(ylim=(10 ** -7, 10 ** -3))
    #
    # # plt.show()
    # g1.savefig(output_dir + "/All_Mutations_point_plot", dpi=300)
    # plt.close()
    #
    # g2 = sns.catplot(x="label", y="frac_and_weight", data=data_filter, hue="Mutation", order=label_order,
    #                  palette="tab20"
    #                  , kind="point", dodge=True, hue_order=transition_order, join=False, estimator=weighted_varaint,
    #                  orient="v", legend=True)
    # g2.set_axis_labels("", "Variant Frequency")
    # g2.set(yscale='log')
    # g2.set(ylim=(10 ** -6, 10 ** -2))
    # # g2.set_yticklabels(fontsize=12)
    # g2.set_xticklabels(fontsize=9, rotation=90)
    # # plt.show()
    # # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/Transition_Mutations_point_plot_RV", dpi=300)
    # g2.savefig(output_dir + "/Transition_Mutations_point_plot", dpi=300)
    # plt.close()
    #
    # data_filter["passage"] = data_filter["passage"].astype(str)
    # passage_g = sns.catplot(x="passage", y="frac_and_weight", data=data_filter, hue="Mutation", order=passage_order,
    #                         palette="tab20"
    #                         , kind="point", dodge=True, hue_order=transition_order, join=False,
    #                         estimator=weighted_varaint,
    #                         orient="v", legend=True)
    # passage_g.set_axis_labels("", "Variant Frequency")
    # passage_g.set(yscale='log')
    # passage_g.set(ylim=(10 ** -6, 10 ** -2))
    # # g2.set_yticklabels(fontsize=12)
    # passage_g.set_xticklabels(fontsize=10, rotation=45)
    # # plt.show()
    # # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/Transition_Mutations_point_plot_RV", dpi=300)
    # passage_g.savefig(output_dir + "/Transition_Mutations_point_plot_all_together", dpi=300)
    # plt.close()
    # data_filter["passage"] = data_filter["passage"].astype(int)
    #
    # data_pmsc = data_filter[data_filter["Type"] == "Premature Stop Codon"]
    # data_pmsc["mutation_type"] = data_pmsc.Mutation.str.contains("A>G") | data_pmsc.Mutation.str.contains("U>C") | \
    #                              data_pmsc.Mutation.str.contains("C>U") | data_pmsc.Mutation.str.contains("G>A")
    # data_pmsc_transition = data_pmsc[data_pmsc["mutation_type"] == True]
    # g3 = sns.catplot("label", "frac_and_weight", data=data_pmsc_transition, hue="Mutation", order=label_order,
    #                  palette="tab20", estimator=weighted_varaint, orient="v", dodge=True, kind="point", col="Type",
    #                  join=False)
    # g3.set_axis_labels("Sample", "Variant Frequency")
    # g3.set(yscale='log')
    # g3.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -5)
    # g3.set_xticklabels(rotation=45)
    # # plt.show()
    # g3.savefig(output_dir + "/Transitions_PMSC_Mutations_point_plot", dpi=300)
    # plt.close()
    #
    # g4 = sns.relplot("passage", "frac_and_weight", data=data_filter, hue="Mutation", palette="tab20",
    #                  hue_order=transition_order, estimator=weighted_varaint, col="Type", kind="line",
    #                  col_order=type_order)
    #
    # g4.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -5)
    # g4.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g4.savefig(output_dir + "/Time_Transition_Mutations_line_plot", dpi=300)
    # plt.close()
    #
    # # A>G Prev Context
    # g5 = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=label_order,
    #                  palette=flatui,
    #                  kind="point", dodge=True, hue_order=[True, False], estimator=weighted_varaint, orient="v",
    #                  col="Type", join=False, col_order=type_order)
    # g5.set_axis_labels("", "Variant Frequency")
    # g5.set(yscale='log')
    # g5.set(ylim=(7 * 10 ** -7, 4 * 10 ** -3))
    # g5.set_xticklabels(fontsize=9, rotation=90)
    # # plt.show()
    # g5.savefig(output_dir + "/Context_point_plot", dpi=300)
    # plt.close()
    #
    # adar_all_g = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=label_order,
    #                          palette=flatui, kind="point", dodge=True, hue_order=[True, False],
    #                          estimator=weighted_varaint, orient="v", join=False)
    # adar_all_g.set_axis_labels("", "Variant Frequency")
    # adar_all_g.set(yscale='log')
    # adar_all_g.set(ylim=(7 * 10 ** -7, 4 * 10 ** -3))
    # adar_all_g.set_xticklabels(rotation=45)
    # # plt.show()
    # adar_all_g.savefig(output_dir + "/Context_point_all_mutations_plot", dpi=300)
    # plt.close()
    # data_filter_ag["passage"] = data_filter_ag["passage"].astype(int)
    #
    # g_context = sns.catplot("ADAR_like", "frac_and_weight", data=data_filter_ag_pass5,
    #                         order=[True, False], palette="Set2", kind="point",
    #                         join=False, estimator=weighted_varaint, orient="v", dodge=True)
    # # add_stat_annotation(g_context, data=data_filter_ag_pass5, x="ADAR_like", y="Frequency", order=[True, False],
    # #                     boxPairList=[(True, False)], test='Mann-Whitney', textFormat='star', loc='inside', verbose=2)
    # g_context.set_axis_labels("ADAR-like\nContext", "Variant Frequency")
    # g_context.set_xticklabels(rotation=45)
    # g_context.set(yscale='log')
    # g_context.set(ylim=(5 * 10 ** -5, 5 * 10 ** -3))
    # plt.tight_layout()
    # # plt.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/context_p5_point_plot.png", dpi=300)
    # plt.savefig(output_dir + "/context_p5_point_plot.png", dpi=300)
    # plt.close()
    #
    # ax = sns.boxplot("ADAR_like", "Frequency", data=data_filter_ag_pass5, palette=flatui, order=(True, False))
    # ax = sns.stripplot("ADAR_like", "Frequency", data=data_filter_ag_pass5, color=".2", order=(True, False))
    # add_stat_annotation(ax, data=data_filter_ag_pass5, x="ADAR_like", y="Frequency",
    #                     boxPairList=[(True, False)], test='Mann-Whitney', textFormat='star', loc='outside', verbose=2)
    # ax.set_yscale('log')
    # ax.set_xlabel("ADAR-like\nContext")
    # ax.set_ylabel("Variant Frequency")
    # ax.set(ylim=(10 ** -4, 10 ** -2))
    # sns.despine()
    # plt.tight_layout()
    # # plt.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/context_p5_point_plot_v2.png", dpi=300)
    # plt.savefig(output_dir + "/context_p5_point_plot_v2.png", dpi=300)
    # plt.close()
    #
    # data_codons = data_filter_ag[data_filter_ag["pval"] < 0.01]
    # data_codons = data_codons[data_codons["ADAR_like"] == True]
    # data_codons = data_codons[data_codons["Type"] == "Synonymous"]
    # g_codons = sns.catplot("label", "frac_and_weight", data=data_codons, hue="Consensus>Mutated_codon",
    #                        order=label_order, palette="tab20", kind="point",
    #                        join=False, estimator=weighted_varaint, orient="v", dodge=True)
    # g_codons.set_axis_labels("", "Variant Frequency")
    # g_codons.set_xticklabels(fontsize=5, rotation=45)
    # g_codons.set(yscale='log')
    # # g_codons.set(ylim=(10**-7, 10**-3))
    #
    # # plt.show()
    # g_codons.savefig(output_dir + "/codons_point_plot", dpi=300)
    # plt.close()
    #
    # g_codons_3c = sns.catplot("label", "frac_and_weight", data=data_codons[data_codons["Protein"] == "3C"],
    #                           hue="Consensus>Mutated_codon",
    #                           order=label_order, palette="tab20", kind="point",
    #                           join=False, estimator=weighted_varaint, orient="v", dodge=True)
    # g_codons_3c.set_axis_labels("", "Variant Frequency")
    # g_codons_3c.set_xticklabels(fontsize=5, rotation=45)
    # g_codons_3c.set(yscale='log')
    # # g_codons.set(ylim=(10**-7, 10**-3))
    #
    # # plt.show()
    # g_codons_3c.savefig(output_dir + "/codons_point_plot_3c", dpi=300)
    # plt.close()
    #
    # g6 = sns.relplot("passage", "frac_and_weight", data=data_filter_ag, hue="Context", palette="tab20",
    #                  hue_order=context_order, estimator=weighted_varaint, col="Type", kind="line",
    #                  col_order=type_order)
    #
    # g6.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -5)
    # g6.set(ylim=(0, 10 ** -2))
    # yaxis = plt.gca().yaxis
    # yaxis.set_minor_locator(fits_new_plotter.MinorSymLogLocator(1e-1))
    # g6.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g6.savefig(output_dir + "/Context_miseq_sample_time_line_plot", dpi=300)
    # plt.close()
    #
    # g7 = sns.relplot(x="passage", y="frac_and_weight", hue="Context", data=data_filter_ag, palette="Paired",
    #                  kind="line",
    #                  style="Type", style_order=type_order, hue_order=context_order, estimator=weighted_varaint)
    # g7.set(yscale="log")
    # g7.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    # g7.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g7.savefig(output_dir + "/ADAR_like_AG_Mutation_Context_trajectories_Context_.png", dpi=300)
    # plt.close()
    #
    # g8 = sns.relplot(x="passage", y="frac_and_weight", hue="ADAR_like", data=data_filter_ag, palette="Paired",
    #                  kind="line", style="Type", style_order=type_order,
    #                  estimator=weighted_varaint)  # , hue_order=context_order)
    # g8.set(yscale="log")
    # g8.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    # g8.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g8.savefig(output_dir + "/ADAR_like_AG_Mutation_ADAR_trajectories_Context.png", dpi=300)
    # plt.close()
    #
    # mutation_g8 = sns.relplot(x="passage", y="frac_and_weight", hue="ADAR_like", data=data_filter_ag, palette=flatui,
    #                           kind="line", col="Type", col_order=type_order,
    #                           estimator=weighted_varaint)  # , hue_order=context_order)
    # mutation_g8.set(yscale="log")
    # mutation_g8.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    # mutation_g8.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # mutation_g8.savefig(output_dir + "/ADAR_like_AG_Mutation_ADAR_trajectories_Context_mutation_col.png", dpi=300)
    # plt.close()
    #
    # mutation_g9 = sns.relplot(x="passage", y="frac_and_weight", hue="ADAR_like", data=data_filter_ag, palette=flatui,
    #                           kind="line", estimator=weighted_varaint)  # , hue_order=context_order)
    # mutation_g9.set(yscale="log")
    # mutation_g9.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    # mutation_g9.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # mutation_g9.savefig(output_dir + "/ADAR_like_AG_Mutation_ADAR_trajectories_Context_all_mutation.png", dpi=300)
    # plt.close()
    #
    # g9 = sns.catplot("label", "frac_and_weight", data=data_filter_uc, hue="Next", order=label_order, palette="tab20",
    #                  hue_order=context_order_uc, estimator=weighted_varaint, orient="v", dodge=True, kind="point",
    #                  col="Type", join=False, col_order=type_order)
    # g9.set_axis_labels("", "Variant Frequency")
    # g9.set(yscale='log')
    # g9.set(ylim=(10 ** -7, 10 ** -3))
    # g9.set_xticklabels(rotation=45)
    # # plt.show()
    # g9.savefig(output_dir + "/UC_Context_point_plot", dpi=300)
    # plt.close()
    #
    # g10 = sns.relplot("passage", "frac_and_weight", data=data_filter_uc, hue="Next", palette="tab20",
    #                   hue_order=context_order_uc, estimator=weighted_varaint, col="Type", kind="line",
    #                   col_order=type_order)
    #
    # # g8.set(yscale="log")
    # g10.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -4)
    # # g8.set(ylim=(0, 10 ** -6))
    # # yaxis = plt.gca().yaxis
    # yaxis.set_minor_locator(fits_new_plotter.MinorSymLogLocator(1e-1))
    # g10.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g10.savefig(output_dir + "/UC_Context_miseq_sample_time_line_plot", dpi=300)
    # plt.close()
    #
    # data_filter_ag_grouped = data_filter_ag.groupby(["ADAR_like", "label", "Type", "Pos", "Protein", "passage"])[
    #     "frac_and_weight"].agg(lambda x: weighted_varaint(x))
    # data_filter_ag_grouped = data_filter_ag_grouped.reset_index()
    # # data_filter_ag_grouped["passage"] = data_filter_ag_grouped["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    # print(data_filter_ag_grouped.to_string())
    #
    # # data_filter_ag_grouped["frac_and_weight"] = np.log10(data_filter_ag_grouped["frac_and_weight"])
    #
    # data_filter_ag_grouped = data_filter_ag_grouped.rename(columns={"frac_and_weight": "Frequency"})
    # data_filter_ag_grouped["Frequency"] = data_filter_ag_grouped["Frequency"].astype(float)
    #
    # data_filter_ag_grouped = data_filter_ag_grouped[data_filter_ag_grouped["label"] != "RNA Control\nPrimer ID"]
    # data_filter_ag_grouped = data_filter_ag_grouped[data_filter_ag_grouped["label"] != "RNA Control\nRND"]
    # # data_filter_ag_grouped["passage"] = np.where(data_filter_ag_grouped["passage"] == "RNA Control", 0, data_filter_ag_grouped["passage"])
    # # data_filter_ag_grouped["passage"] = data_filter_ag_grouped["passage"].astype(int)
    #
    # data_reg_adar = data_filter_ag_grouped[data_filter_ag_grouped["ADAR_like"] == True]
    # data_reg_nonadar = data_filter_ag_grouped[data_filter_ag_grouped["ADAR_like"] == False]
    #
    # data_reg_adar_syn = data_reg_adar[data_reg_adar["Type"] == "Synonymous"]
    # data_reg_nonadar_syn = data_reg_nonadar[data_reg_nonadar["Type"] == "Synonymous"]
    #
    # data_reg_adar_nonsyn = data_reg_adar[data_reg_adar["Type"] == "Non-Synonymous"]
    # data_reg_nonadar_non_syn = data_reg_nonadar[data_reg_nonadar["Type"] == "Non-Synonymous"]
    #
    # fig, axes = plt.subplots(2, 2, sharey=True, sharex=True)
    # stat_slope1, stat_intercept1, r_value1, p_value1, std_err1 = stats.linregress(data_reg_adar_syn['passage'],
    #                                                                               data_reg_adar_syn[
    #                                                                                   'Frequency'])  # passage
    # print("stats - ADAR SYN slope: %s" % stat_slope1)
    # print("stats - ADAR SYN intercept: %s" % stat_intercept1)
    #
    # # data_reg_adar_nonsyn
    # stat_slope2, stat_intercept2, r_value1, p_value1, std_err1 = stats.linregress(data_reg_adar_nonsyn['passage'],
    #                                                                               data_reg_adar_nonsyn['Frequency'])
    #
    # print("stats - ADAR NON-SYN slope: %s" % stat_slope2)
    # print("stats - ADAR NON-SYN intercept: %s" % stat_intercept2)
    #
    # from sklearn import linear_model
    # regr = linear_model.LinearRegression()
    # X = data_reg_adar_syn.passage.values.reshape(-1, 1)
    # y = data_reg_adar_syn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope1 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept1 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("ADAR SYN slope: %s" % str(regr.coef_[0]).split("[")[1].split("]")[0])
    # print("ADAR SYN intercept: %s" % str(regr.intercept_).split("[")[1].split("]")[0])
    #
    # g1 = sns.regplot(x="passage", y="Frequency", data=data_reg_adar_syn, ax=axes[0, 0],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope1, intercept1)})
    # g1.set(title="ADAR-like Synonymous")
    #
    # X = data_reg_nonadar_syn.passage.values.reshape(-1, 1)
    # y = data_reg_nonadar_syn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope2 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept2 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("NON-ADAR SYN slope: %s" % str(regr.coef_[0]).split("[")[1].split("]")[0])
    # print("NON-ADAR SYN intercept: %s" % str(regr.intercept_).split("[")[1].split("]")[0])
    #
    # g2 = sns.regplot(x="passage", y="Frequency", data=data_reg_nonadar_syn, ax=axes[0, 1],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope2, intercept2)})
    # g2.set(title="Non ADAR-like Synonymous")
    #
    # X = data_reg_adar_nonsyn.passage.values.reshape(-1, 1)
    # y = data_reg_adar_nonsyn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope3 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept3 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("ADAR NON-SYN slope: %s" % str(regr.coef_[0]).split("[")[1].split("]")[0])
    # print("ADAR NON-SYN intercept: %s" % str(regr.intercept_).split("[")[1].split("]")[0])
    #
    # g3 = sns.regplot(x="passage", y="Frequency", data=data_reg_adar_nonsyn, ax=axes[1, 0],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope3, intercept3)})
    # g3.set(title="ADAR-like Non Synonymous")
    #
    # X = data_reg_nonadar_non_syn.passage.values.reshape(-1, 1)
    # y = data_reg_nonadar_non_syn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope4 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept4 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("NON-ADAR NON-SYN slope: %s" % str(regr.coef_[0]).split("[")[1].split("]")[0])
    # print("NON-ADAR NON-SYN intercept: %s" % str(regr.intercept_).split("[")[1].split("]")[0])
    #
    # g4 = sns.regplot(x="passage", y="Frequency", data=data_reg_nonadar_non_syn, ax=axes[1, 1],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope4, intercept4)})
    # g4.set(title="Non ADAR-like Non Synonymous")
    #
    # axes[0, 0].legend()
    # axes[0, 0].legend(loc=2)
    # axes[0, 0].set_yscale("log")
    # axes[0, 0].set_xlabel('')
    # axes[0, 0].set_ylabel('')
    # axes[0, 0].set_xlim(0, 14, 2)
    # axes[0, 0].set_ylim(10 ** -6, 10 ** -2)
    # axes[0, 1].legend()
    # axes[0, 1].legend(loc=2)
    # axes[0, 1].set_yscale("log")
    # axes[0, 1].set_xlabel('')
    # axes[0, 1].set_ylabel('')
    # axes[0, 1].set_xlim(0, 14, 2)
    # axes[0, 1].set_ylim(10 ** -6, 10 ** -2)
    # axes[1, 0].legend()
    # axes[1, 0].legend(loc=2)
    # axes[1, 0].set_yscale("log")
    # axes[1, 0].set_xlabel('')
    # axes[1, 0].set_ylabel('')
    # axes[1, 0].set_xlim(0, 14, 2)
    # axes[1, 0].set_ylim(10 ** -6, 10 ** -2)
    # axes[1, 1].legend()
    # axes[1, 1].legend(loc=2)
    # axes[1, 1].set_yscale("log")
    # axes[1, 1].set_xlabel('')
    # axes[1, 1].set_ylabel('')
    # axes[1, 1].set_xlim(0, 14, 2)
    # axes[1, 1].set_ylim(10 ** -6, 10 ** -2)
    # #
    # # for ax in axes:
    # #     ax.set_yscale("log")
    # #     ax.set_ylim(10**-4,10**-2)
    # #     # # l = ax.get_ylabel()
    # #     # # ax.set_ylabel(l, fontsize=8)
    # #     ax.set_xlabel('')
    # #     ax.set_ylabel('')
    # fig.text(0.5, 0.01, 'Population', ha='center')
    # fig.text(0.001, 0.5, 'Frequency', va='center', rotation='vertical')
    # # fig.suptitle("Mutation trajectories", fontsize=14)
    # fig.tight_layout()
    # # plt.show()
    # plt.savefig(output_dir + "/regplot_AG_Mutation_Context_trajectories.png", dpi=300)
    # plt.close()
    #
    # X = data_reg_adar_nonsyn.passage.values.reshape(-1, 1)
    # y = data_reg_adar_nonsyn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # print("ADAR NON-SYN slope: %s" % str(regr.coef_[0]))
    # print("ADAR NON-SYN intercept: %s" % str(regr.intercept_))
    #
    # g11 = sns.lmplot(x="passage", y="Frequency", data=data_filter_ag_grouped, hue="ADAR_like", markers=["o", "x"],
    #                  hue_order=[True, False], fit_reg=True, col="Type", col_order=type_order, palette=flatui)
    # g11.set(xlim=(0, 13))
    # # g11.set(yscale="log")
    # # g11.set(ylim=(10**-6, 10**-2))
    # # props = dict(boxstyle='round', alpha=0.5, color=sns.color_palette()[0])
    # # textstr = "ADAR-like SYN: y={0:.3g}x+{1:.3g}\nnon ADAR-like SYN: y={2:.3g}x+{3:.3g}".format(slope1, intercept1, slope2, intercept2)
    # # g11.ax.text(0.7, 0.9, textstr, transform=g11.axes.transAxes, fontsize=14, bbox=props)
    # # plt.show()
    # g11.savefig(output_dir + "/lmplot_ADAR_Context", dpi=300)
    # plt.close()
    #
    # data_filter_ag_grouped_silent = data_filter_ag_grouped[data_filter_ag_grouped["Type"] == "Synonymous"]
    # data_filter_ag_grouped_silent = data_filter_ag_grouped_silent[data_filter_ag_grouped_silent["Protein"] != "2A"]
    # position_g = sns.scatterplot("Pos", "Frequency", data=data_filter_ag_grouped_silent, hue="Protein",
    #                              palette="tab10", style="ADAR_like", style_order=[True, False], legend="full")
    #
    # # position_g.set_axis_labels("", "Variant Frequency")
    # position_g.set_yscale('log')
    # position_g.set(xlim=(3500, 7500))
    # position_g.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # plt.tight_layout()
    # # position_g.set(ylim=(10 ** -6, 10 ** -2))
    # # g2.set_yticklabels(fontsize=12)
    # # position_g.set_xticklabels(fontsize=10, rotation=45)
    # # plt.show()
    # # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/Transition_Mutations_point_plot_RV", dpi=300)
    # plt.savefig(output_dir + "/position.png", dpi=300)
    # plt.close()
    #
    data_filter["passage"] = data_filter["passage"].astype(str)
    position_mutation = sns.catplot(x="Pos", y="Frequency", data=data_filter, hue="Mutation", palette="tab20"
                    , kind="point", dodge=True, hue_order=transition_order, join=False,
                                    legend=True, col="passage", col_order=passage_order, col_wrap=3)

    position_mutation.set_axis_labels("", "Variant Frequency")
    position_mutation.set_yscale('log')
    position_mutation.set(xlim=(3500, 7500))
    position_mutation.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # plt.tight_layout()
    # position_g.set(ylim=(10 ** -6, 10 ** -2))
    # g2.set_yticklabels(fontsize=12)
    # position_g.set_xticklabels(fontsize=10, rotation=45)
    # plt.show()

    plt.savefig(output_dir + "/position_mutation.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    main()

