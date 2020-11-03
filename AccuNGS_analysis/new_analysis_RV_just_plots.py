
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
from AccuNGS_analysis import old_statannot

sns.set(font_scale=1.2)
sns.set_style("ticks")
sns.despine()


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
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/"
    prefix = "inosine_predict_context"
    output_dir = input_dir + "20201103_%s" % prefix
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    data_filter = pd.read_pickle(input_dir + prefix + "/data_filter.pkl")
    data_filter_ag = pd.read_pickle(input_dir + prefix + "/data_filter_ag.pkl")
    data_filter_uc = pd.read_pickle(input_dir + prefix +"/data_filter_uc.pkl")
    data_filter["passage"] = data_filter["passage"].astype(int)

    #Plots
    label_order = ["RNA Control\nRND", "RNA Control\nPrimer ID","p2-1", "p2-2", "p2-3", "p5-1", "p5-2", "p5-3", "p8-1",
                   "p8-2", "p8-3", "p10-2", "p10-3", "p12-1", "p12-2", "p12-3"]
    passage_order = ["0", "2", "5", "8", "10", "12"]
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U"]
    transition_order = ["A>G", "U>C", "G>A", "C>U"]
    type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]
    type_order_ag = ["Synonymous", "Non-Synonymous"]
    context_order = ["UpA", "ApA", "CpA", "GpA"]
    context_order_uc = ["UpU", "UpA", "UpC", "UpG"]
    adar_preference = ["High", "Intermediate", "Low"]

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

    # data_filter["passage"] = data_filter["passage"].astype(str)
    # passage_g = sns.catplot(x="passage", y="frac_and_weight", data=data_filter, hue="Mutation", order=passage_order,
    #                         palette="tab20"
    #                         , kind="point", dodge=True, hue_order=transition_order, join=False,
    #                         estimator=weighted_varaint,
    #                         orient="v", legend=True)
    # passage_g.set_axis_labels("Passage", "Variant Frequency")
    # passage_g.set(yscale='log')
    # passage_g.set(ylim=(10 ** -6, 10 ** -2))
    # # g2.set_yticklabels(fontsize=12)
    # # passage_g.set_xticklabels(fontsize=10, rotation=45)
    # # plt.show()
    # # passage_g.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/Prgress reports/20200913 Final report/plots" +
    # #                   "/Transition_Mutations_point_plot_RVB14", dpi=300)
    # passage_g.savefig(output_dir + "/Transition_Mutations_point_plot_RVB14", dpi=300)
    # plt.close()
    # data_filter["passage"] = data_filter["passage"].astype(int)

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
    #                  palette=flatui, kind="point", dodge=True, hue_order=[True, False], estimator=weighted_varaint,
    #                  orient="v", col="Type", join=False, col_order=type_order_ag)
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
    # data_filter_ag_pass5 = data_filter_ag.loc[data_filter_ag.passage == 5]
    # data_filter_ag_pass5 = data_filter_ag_pass5[data_filter_ag_pass5["pval"] < 0.01]
    # data_filter_ag_pass5 = data_filter_ag_pass5.loc[data_filter_ag_pass5.Type == "Synonymous"]
    #
    # g_context = sns.catplot("ADAR_like", "frac_and_weight", data=data_filter_ag_pass5,
    #                         order=[True, False], palette="rocket", kind="point",
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


    mutation_g8 = sns.catplot("passage", "frac_and_weight", data=data_filter_ag, hue="5`_ADAR_Preference",
                     palette="rocket", kind="point", dodge=True, estimator=weighted_varaint,
                     orient="v", col="Type", join=False, col_order=type_order_ag, hue_order=adar_preference)
    mutation_g8.set(yscale="log")
    # mutation_g8.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    mutation_g8.set_axis_labels("Passage", "Variant Frequency")
    # plt.show()
    mutation_g8.savefig(output_dir + "/ADAR_like_AG_Mutation_col.png", dpi=300)
    plt.close()


    data_filter_ag_grouped = data_filter_ag.groupby(["label", "Type", "Pos", "Protein", "passage",
                                                     "ADAR_grade_five", "5`_ADAR_Preference"])["frac_and_weight"].agg(
        lambda x: weighted_varaint(x))
    data_filter_ag_grouped = data_filter_ag_grouped.reset_index()

    print(data_filter_ag_grouped.to_string())

    data_filter_ag_grouped = data_filter_ag_grouped.rename(columns={"frac_and_weight": "Frequency"})
    data_filter_ag_grouped["Frequency"] = data_filter_ag_grouped["Frequency"].astype(float)
    data_filter_ag_grouped = data_filter_ag_grouped[data_filter_ag_grouped["label"] != "RNA Control\nPrimer ID"]
    data_filter_ag_grouped = data_filter_ag_grouped[data_filter_ag_grouped["label"] != "RNA Control\nRND"]

    data_reg_full_adar = data_filter_ag_grouped[data_filter_ag_grouped["ADAR_grade_five"] == 1]
    data_reg_semi_adar = data_filter_ag_grouped[data_filter_ag_grouped["ADAR_grade_five"] == 0.5]
    data_reg_nonadar = data_filter_ag_grouped[data_filter_ag_grouped["ADAR_grade_five"] == 0]

    data_reg_full_adar_syn = data_reg_full_adar[data_reg_full_adar["Type"] == "Synonymous"]
    data_reg_semi_adar_syn = data_reg_semi_adar[data_reg_semi_adar["Type"] == "Synonymous"]
    data_reg_nonadar_syn = data_reg_nonadar[data_reg_nonadar["Type"] == "Synonymous"]

    data_reg_full_adar_non_syn = data_reg_full_adar[data_reg_full_adar["Type"] == "Non-Synonymous"]
    data_reg_semi_adar_non_syn = data_reg_semi_adar[data_reg_semi_adar["Type"] == "Non-Synonymous"]
    data_reg_nonadar_non_syn = data_reg_nonadar[data_reg_nonadar["Type"] == "Non-Synonymous"]

    # fig, axes = plt.subplots(2, 3, sharey="all", sharex="all")
    #
    # from sklearn import linear_model
    # regr = linear_model.LinearRegression()
    # X = data_reg_full_adar_syn.passage.values.reshape(-1, 1)
    # y = data_reg_full_adar_syn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope1 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept1 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("100-ADAR SYN: y=%sx + %s" % (str(regr.coef_[0]).split("[")[1].split("]")[0], str(regr.intercept_).split("[")[1].split("]")[0]))
    #
    # g1 = sns.regplot(x="passage", y="Frequency", data=data_reg_full_adar_syn, ax=axes[0, 0],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope1, intercept1)})
    # g1.set(title="100-ADAR-like Synonymous")
    #
    # X = data_reg_semi_adar_syn.passage.values.reshape(-1, 1)
    # y = data_reg_semi_adar_syn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope2 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept2 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("50-ADAR SYN: y=%sx + %s" % (str(regr.coef_[0]).split("[")[1].split("]")[0], str(regr.intercept_).split("[")[1].split("]")[0]))
    #
    # g2 = sns.regplot(x="passage", y="Frequency", data=data_reg_semi_adar_syn, ax=axes[0, 1],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope2, intercept2)})
    # g2.set(title="50-ADAR-like Synonymous")
    #
    # X = data_reg_nonadar_syn.passage.values.reshape(-1, 1)
    # y = data_reg_nonadar_syn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope3 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept3 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("NON-ADAR SYN: y=%sx + %s" % (str(regr.coef_[0]).split("[")[1].split("]")[0], str(regr.intercept_).split("[")[1].split("]")[0]))
    #
    # g3 = sns.regplot(x="passage", y="Frequency", data=data_reg_nonadar_syn, ax=axes[0, 2],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope3, intercept3)})
    # g3.set(title="Non ADAR-like Synonymous")
    #
    # X = data_reg_full_adar_non_syn.passage.values.reshape(-1, 1)
    # y = data_reg_full_adar_non_syn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope4 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept4 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("100-ADAR NON-SYN: y=%sx + %s" % (str(regr.coef_[0]).split("[")[1].split("]")[0], str(regr.intercept_).split("[")[1].split("]")[0]))
    #
    #
    # g4 = sns.regplot(x="passage", y="Frequency", data=data_reg_full_adar_non_syn, ax=axes[1, 0],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope4, intercept4)})
    # g4.set(title="100-ADAR-like Non Synonymous")
    #
    # X = data_reg_semi_adar_non_syn.passage.values.reshape(-1, 1)
    # y = data_reg_semi_adar_non_syn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope5 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept5 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("50-ADAR NON-SYN: y=%sx + %s" % (str(regr.coef_[0]).split("[")[1].split("]")[0], str(regr.intercept_).split("[")[1].split("]")[0]))
    #
    # g5 = sns.regplot(x="passage", y="Frequency", data=data_reg_semi_adar_non_syn, ax=axes[1, 1],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope5, intercept5)})
    # g5.set(title="50-ADAR-like Non Synonymous")
    #
    # X = data_reg_nonadar_non_syn.passage.values.reshape(-1, 1)
    # y = data_reg_nonadar_non_syn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope6 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept6 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("NON-ADAR NON-SYN: y=%sx + %s" % (str(regr.coef_[0]).split("[")[1].split("]")[0], str(regr.intercept_).split("[")[1].split("]")[0]))
    #
    # g6 = sns.regplot(x="passage", y="Frequency", data=data_reg_nonadar_non_syn, ax=axes[1, 2],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope6, intercept6)})
    # g6.set(title="Non ADAR-like Non Synonymous")
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
    # axes[0, 2].legend()
    # axes[0, 2].legend(loc=2)
    # axes[0, 2].set_yscale("log")
    # axes[0, 2].set_xlabel('')
    # axes[0, 2].set_ylabel('')
    # axes[0, 2].set_xlim(0, 14, 2)
    # axes[0, 2].set_ylim(10 ** -6, 10 ** -2)
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
    # axes[1, 2].legend()
    # axes[1, 2].legend(loc=2)
    # axes[1, 2].set_yscale("log")
    # axes[1, 2].set_xlabel('')
    # axes[1, 2].set_ylabel('')
    # axes[1, 2].set_xlim(0, 14, 2)
    # axes[1, 2].set_ylim(10 ** -6, 10 ** -2)
    # fig.text(0.5, 0.01, 'Passage', ha='center')
    # fig.text(0.001, 0.5, 'Frequency', va='center', rotation='vertical')
    # plt.savefig(output_dir + "/regplot_AG_Mutation_Context_trajectories.png", dpi=300)
    # plt.close()

    stat_slope1, stat_intercept1, r_value1, p_value1, std_err1 = stats.linregress(data_reg_full_adar_syn['passage'],
                                                                                  data_reg_full_adar_syn[
                                                                                      'Frequency'])
    stat_slope2, stat_intercept2, r_value2, p_value2, std_err2 = stats.linregress(data_reg_semi_adar_syn['passage'],
                                                                                  data_reg_semi_adar_syn[
                                                                                      'Frequency'])
    stat_slope3, stat_intercept3, r_value3, p_value3, std_err3 = stats.linregress(data_reg_nonadar_syn['passage'],
                                                                                  data_reg_nonadar_syn[
                                                                                      'Frequency'])
    # data_reg_adar_nonsyn
    stat_slope4, stat_intercept4, r_value4, p_value4, std_err4 = stats.linregress(data_reg_full_adar_non_syn['passage'],
                                                                                  data_reg_full_adar_non_syn[
                                                                                      'Frequency'])
    stat_slope5, stat_intercept5, r_value5, p_value5, std_err5 = stats.linregress(data_reg_semi_adar_non_syn['passage'],
                                                                                  data_reg_semi_adar_non_syn[
                                                                                      'Frequency'])
    stat_slope6, stat_intercept6, r_value6, p_value6, std_err6 = stats.linregress(data_reg_nonadar_non_syn['passage'],
                                                                                  data_reg_nonadar_non_syn[
                                                                                      'Frequency'])
    data_filter_ag_grouped = data_filter_ag_grouped.rename(columns={"passage": "Passage"})
    g11 = sns.lmplot(x="Passage", y="Frequency", data=data_filter_ag_grouped, hue="5`_ADAR_Preference",
                     hue_order=adar_preference, markers=["o", "v", "x"], fit_reg=True, col="5`_ADAR_Preference",
                     col_order=adar_preference, row="Type", row_order=type_order_ag, palette="rocket",
                     line_kws={'label': "Linear Reg"}, legend=True, height=6)
    g11.fig.subplots_adjust(wspace=.02)
    ax = g11.axes[0, 0]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    label_line_1 = "y={0:.3g}x+{1:.3g}".format(stat_slope1, stat_intercept1)
    label_line_2 = "y={0:.3g}x+{1:.3g}".format(stat_slope2, stat_intercept2)
    label_line_3 = "y={0:.3g}x+{1:.3g}".format(stat_slope3, stat_intercept3)
    label_line_4 = "y={0:.3g}x+{1:.3g}".format(stat_slope4, stat_intercept4)
    label_line_5 = "y={0:.3g}x+{1:.3g}".format(stat_slope5, stat_intercept5)
    label_line_6 = "y={0:.3g}x+{1:.3g}".format(stat_slope6, stat_intercept6)
    L_labels[0].set_text(label_line_1)
    ax = g11.axes[0, 1]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    L_labels[0].set_text(label_line_2)
    ax = g11.axes[0, 2]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    L_labels[0].set_text(label_line_3)
    ax = g11.axes[1, 0]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    L_labels[0].set_text(label_line_4)
    ax = g11.axes[1, 1]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    L_labels[0].set_text(label_line_5)
    ax = g11.axes[1, 2]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    L_labels[0].set_text(label_line_6)

    g11.set(xlim=(0, 13))
    g11.set(ylim=(0.000, 0.010))
    plt.tight_layout()
    g11.savefig(output_dir + "/lmplot_ADAR_Context", dpi=300)
    plt.close()

    data_filter_ag = data_filter_ag[data_filter_ag["Protein"] != "2A"]
    data_filter_ag = data_filter_ag[data_filter_ag["Protein"] != "3'UTR"]
    data_filter_ag = data_filter_ag[data_filter_ag["Type"] == "Synonymous"]

    position_mutation = sns.relplot(x="Pos", y="Frequency", data=data_filter_ag, hue="5`_ADAR_Preference",
                                    col="passage", col_wrap=3, palette="rocket",
                                    hue_order=adar_preference, height=4,
                                    style="5`_ADAR_Preference", style_order=["High", "Low", "Intermediate"])

    position_mutation.set_axis_labels("", "Variant Frequency")
    position_mutation.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -4)
    position_mutation.axes.flat[0].set_ylim(10 ** -5, 10 ** -2)
    plt.savefig(output_dir + "/position_mutation.png", dpi=300)
    plt.close()


    """U>C Context"""
    mutation_uc = sns.catplot("passage", "frac_and_weight", data=data_filter_uc, hue="3`_ADAR_Preference",
                              palette="rocket", kind="point", dodge=True, estimator=weighted_varaint,
                              orient="v", col="Type", join=False, hue_order=adar_preference,
                              col_order=type_order_ag)
    mutation_uc.set(yscale="log")
    mutation_uc.set_axis_labels("Passage", "Variant Frequency")
    # plt.show()
    mutation_uc.savefig(output_dir + "/uc_ADAR_like_Mutation_col.png", dpi=300)
    plt.close()

    data_filter_uc_grouped = data_filter_uc.groupby(["label", "Type", "Pos", "Protein", "passage",
                                                     "ADAR_grade_three", "3`_ADAR_Preference"])["frac_and_weight"].agg(
        lambda x: weighted_varaint(x))
    data_filter_uc_grouped = data_filter_uc_grouped.reset_index()

    print(data_filter_uc_grouped.to_string())

    data_filter_uc_grouped = data_filter_uc_grouped.rename(columns={"frac_and_weight": "Frequency"})
    data_filter_uc_grouped["Frequency"] = data_filter_uc_grouped["Frequency"].astype(float)
    data_filter_uc_grouped = data_filter_uc_grouped[data_filter_uc_grouped["label"] != "RNA Control\nPrimer ID"]
    data_filter_uc_grouped = data_filter_uc_grouped[data_filter_uc_grouped["label"] != "RNA Control\nRND"]

    data_reg_full_adar_uc = data_filter_uc_grouped[data_filter_uc_grouped["ADAR_grade_three"] == 1]
    data_reg_semi_adar_uc = data_filter_uc_grouped[data_filter_uc_grouped["ADAR_grade_three"] == 0.5]
    data_reg_nonadar_uc = data_filter_uc_grouped[data_filter_uc_grouped["ADAR_grade_three"] == 0]

    data_reg_full_adar_syn_uc = data_reg_full_adar_uc[data_reg_full_adar_uc["Type"] == "Synonymous"]
    data_reg_semi_adar_syn_uc = data_reg_semi_adar_uc[data_reg_semi_adar_uc["Type"] == "Synonymous"]
    data_reg_nonadar_syn_uc = data_reg_nonadar_uc[data_reg_nonadar_uc["Type"] == "Synonymous"]

    data_reg_full_adar_non_syn_uc = data_reg_full_adar_uc[data_reg_full_adar_uc["Type"] == "Non-Synonymous"]
    data_reg_semi_adar_non_syn_uc = data_reg_semi_adar_uc[data_reg_semi_adar_uc["Type"] == "Non-Synonymous"]
    data_reg_nonadar_non_syn_uc = data_reg_nonadar_uc[data_reg_nonadar_uc["Type"] == "Non-Synonymous"]

    stat_slope7, stat_intercept7, r_value7, p_value7, std_err7 = stats.linregress(data_reg_full_adar_syn_uc['passage'],
                                                                                  data_reg_full_adar_syn_uc[
                                                                                      'Frequency'])
    stat_slope8, stat_intercept8, r_value8, p_value8, std_err8 = stats.linregress(data_reg_semi_adar_syn_uc['passage'],
                                                                                  data_reg_semi_adar_syn_uc[
                                                                                      'Frequency'])
    stat_slope9, stat_intercept9, r_value9, p_value9, std_err9 = stats.linregress(data_reg_nonadar_syn_uc['passage'],
                                                                                  data_reg_nonadar_syn_uc[
                                                                                      'Frequency'])
    # data_reg_adar_nonsyn
    stat_slope10, stat_intercept10, r_value10, p_value10, std_err10 = stats.linregress(data_reg_full_adar_non_syn_uc['passage'],
                                                                                  data_reg_full_adar_non_syn_uc[
                                                                                      'Frequency'])
    stat_slope11, stat_intercept11, r_value11, p_value11, std_err11 = stats.linregress(data_reg_semi_adar_non_syn_uc['passage'],
                                                                                  data_reg_semi_adar_non_syn_uc[
                                                                                      'Frequency'])
    stat_slope12, stat_intercept12, r_value12, p_value12, std_err12 = stats.linregress(data_reg_nonadar_non_syn_uc['passage'],
                                                                                  data_reg_nonadar_non_syn_uc[
                                                                                      'Frequency'])
    data_filter_uc_grouped = data_filter_uc_grouped.rename(columns={"passage": "Passage"})
    uc_lmplot = sns.lmplot(x="Passage", y="Frequency", data=data_filter_uc_grouped, hue="3`_ADAR_Preference",
                           markers=["o", "v", "x"], hue_order=adar_preference, fit_reg=True, col="3`_ADAR_Preference",
                           col_order=adar_preference, row="Type", row_order=type_order_ag, palette="rocket",
                           line_kws={'label': "Linear Reg"}, legend=True, height=6)
    uc_lmplot.fig.subplots_adjust(wspace=.02)
    ax = uc_lmplot.axes[0, 0]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    label_line_7 = "y={0:.3g}x+{1:.3g}".format(stat_slope7, stat_intercept7)
    label_line_8 = "y={0:.3g}x+{1:.3g}".format(stat_slope8, stat_intercept8)
    label_line_9 = "y={0:.3g}x+{1:.3g}".format(stat_slope9, stat_intercept9)
    label_line_10 = "y={0:.3g}x+{1:.3g}".format(stat_slope10, stat_intercept10)
    label_line_11 = "y={0:.3g}x+{1:.3g}".format(stat_slope11, stat_intercept11)
    label_line_12 = "y={0:.3g}x+{1:.3g}".format(stat_slope12, stat_intercept12)
    L_labels[0].set_text(label_line_7)
    ax = uc_lmplot.axes[0, 1]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    L_labels[0].set_text(label_line_8)
    ax = uc_lmplot.axes[0, 2]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    L_labels[0].set_text(label_line_9)
    ax = uc_lmplot.axes[1, 0]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    L_labels[0].set_text(label_line_10)
    ax = uc_lmplot.axes[1, 1]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    L_labels[0].set_text(label_line_11)
    ax = uc_lmplot.axes[1, 2]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    L_labels[0].set_text(label_line_12)

    uc_lmplot.set(xlim=(0, 13))
    uc_lmplot.set(ylim=(0.000, 0.010))
    plt.tight_layout()
    uc_lmplot.savefig(output_dir + "/uc_lmplot_ADAR_Context", dpi=300)
    plt.close()

    data_filter_uc = data_filter_uc[data_filter_uc["Protein"] != "2A"]
    data_filter_uc = data_filter_uc[data_filter_uc["Protein"] != "3'UTR"]
    data_filter_uc = data_filter_uc[data_filter_uc["Type"] == "Synonymous"]

    position_mutation_uc = sns.relplot(x="Pos", y="Frequency", data=data_filter_uc, hue="3`_ADAR_Preference",
                                    col="passage", col_wrap=3, palette="rocket",
                                    hue_order=adar_preference, height=4,
                                       style="3`_ADAR_Preference", style_order=["High", "Low", "Intermediate"])

    position_mutation_uc.set_axis_labels("", "Variant Frequency")
    position_mutation_uc.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -4)
    position_mutation_uc.axes.flat[0].set_ylim(10**-5, 10 ** -2)
    plt.savefig(output_dir + "/uc_position_mutation.png", dpi=300)
    plt.close()


    # ax = sns.boxplot("ADAR_like", "Frequency", data=data_filter_ag_pass5, palette="rocket", order=(True, False))
    # ax = sns.stripplot("ADAR_like", "Frequency", data=data_filter_ag_pass5, color=".2", order=(True, False))
    # old_statannot.add_stat_annotation(ax, data=data_filter_ag_pass5, x="ADAR_like", y="Frequency",
    #                     boxPairList=[(True, False)], test='Mann-Whitney', textFormat='star', loc='outside', verbose=2)
    # ax.set_yscale('log')
    # ax.set_xlabel("ADAR-like\nContext")
    # ax.set_ylabel("Variant Frequency")
    # ax.set(ylim=(10 ** -5, 10 ** -2))
    # sns.despine()
    # plt.tight_layout()
    # # plt.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/context_p5_point_plot_v2.png", dpi=300)
    # plt.savefig(output_dir + "/context_p5_box_plot_RVB14.png", dpi=300)
    # plt.close()

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
    # g6 = sns.relplot("passage", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", palette="tab20",
    #                  hue_order=[True, False], estimator=weighted_varaint, col="Type", kind="line",
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

    # g7 = sns.relplot(x="passage", y="frac_and_weight", hue="Prev", data=data_filter_ag, palette="Paired",
    #                  kind="line",style="Type", style_order=type_order, hue_order=context_order,
    #                  estimator=weighted_varaint)
    # g7.set(yscale="log")
    # g7.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    # g7.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g7.savefig(output_dir + "/ADAR_like_AG_Mutation_Context_trajectories_Context_.png", dpi=300)
    # plt.close()

    # g8 = sns.relplot(x="passage", y="frac_and_weight", hue="ADAR_like", data=data_filter_ag, palette="Paired",
    #                  kind="line", style="Type", style_order=type_order,
    #                  estimator=weighted_varaint)  # , hue_order=context_order)
    # g8.set(yscale="log")
    # g8.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    # g8.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g8.savefig(output_dir + "/ADAR_like_AG_Mutation_ADAR_trajectories_Context.png", dpi=300)
    # plt.close()


    # mutation_g9 = sns.relplot(x="passage", y="frac_and_weight", hue="ADAR_like", data=data_filter_ag, palette=flatui,
    #                           kind="line", estimator=weighted_varaint, hue_order=[True, False])
    # mutation_g9.set(yscale="log")
    # mutation_g9.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    # mutation_g9.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # mutation_g9.savefig(output_dir + "/ADAR_like_AG_Mutation_ADAR_trajectories_Context_all_mutation.png", dpi=300)
    # plt.close()

    # g9 = sns.catplot("label", "frac_and_weight", data=data_filter_uc, hue="Next", order=label_order, palette="tab20",
    #                  hue_order=context_order_uc, estimator=weighted_varaint, orient="v", dodge=True, kind="point",
    #                  col="Type", join=False, col_order=type_order_ag)
    # g9.set_axis_labels("", "Variant Frequency")
    # g9.set(yscale='log')
    # g9.set(ylim=(10 ** -7, 10 ** -3))
    # g9.set_xticklabels(rotation=45)
    # # plt.show()
    # g9.savefig(output_dir + "/UC_Context_point_plot", dpi=300)
    # plt.close()

    # g10 = sns.relplot("passage", "frac_and_weight", data=data_filter_uc, hue="Next", palette="tab20",
    #                   hue_order=context_order_uc, estimator=weighted_varaint, col="Type", kind="line",
    #                   col_order=type_order_ag)
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


if __name__ == "__main__":
    main()

