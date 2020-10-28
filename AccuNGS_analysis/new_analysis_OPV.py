
import pandas as pd
import os
import matplotlib.pyplot as plt
from FITS_analysis import fits_new_plotter
import numpy as np
from scipy import stats
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker
import seaborn as sns
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
    input_dir = "/Users/odedkushnir/Projects/fitness/CirSeq/PV/OPV"
    output_dir = input_dir + "/20201027_plots"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)


    data_mutations = pd.read_csv(input_dir + "/q23_data_mutation.csv")

    columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
               "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts", "Consensus>Mutated_codon"]
    data_filter = pd.DataFrame(data_mutations, columns=columns)
    data_filter["pval"] = data_filter["pval"].fillna(1)
    data_filter["no_variants"] = data_filter["Frequency"] * data_filter["Read_count"]
    # filter based on pval<0.01 and Prob>0.95
    # data_filter["no_variants"] = np.where(data_filter["pval"] > 0.01, 0, data_filter["no_variants"])
    data_filter["no_variants"] = np.where(data_filter["Prob"] < 0.95, 0, data_filter["no_variants"])


    data_filter["frac_and_weight"] = list(zip(data_filter.no_variants, data_filter.Read_count))
    # data_filter = data_filter[data_filter["label"] != "RVB14-RNA Control"]
    # data_filter = data_filter[data_filter["label"] != "RVB14-p2"]
    #
    # data_miseq = data_filter[data_filter["label"] != "RVB14-Next-RNA Control"]
    # data_miseq = data_miseq[data_filter["label"] != "RVB14-p1"]

    data_filter["passage"] = data_filter["label"].apply(lambda x: x.split("-")[-1][1])
    # data_filter["passage"] = np.where(data_filter["passage"] == "RNA Control", 0, data_filter["passage"])
    data_filter["passage"] = data_filter["passage"].astype(int)
    data_filter["Type"] = data_filter["Type"].fillna("NonCodingRegion")
    data_filter.to_csv(output_dir + "/data_filter.csv", sep=',', encoding='utf-8')

    label_order = ["OPV-p1", "OPV-p3", "OPV-p5", "OPV-p6", "OPV-p7"]
    passage_order = ["1", "3", "5", "6", "7"]
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U"]
    transition_order = ["A>G", "U>C", "G>A"]
    type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]

    g1 = sns.catplot("label", "frac_and_weight", data=data_filter, hue="Mutation", order=label_order, palette="tab20",
                        kind="point", dodge=True, hue_order=mutation_order, join=False, estimator=weighted_varaint,
                     orient="v")
    g1.set_axis_labels("", "Variant Frequency")
    g1.set_xticklabels(fontsize=9, rotation=45)
    g1.set(yscale='log')
    g1.set(ylim=(10**-7, 10**-3))

    # plt.show()
    g1.savefig(output_dir + "/All_Mutations_point_plot", dpi=300)
    plt.close()

    data_filter["passage"] = data_filter["passage"].astype(str)
    g2 = sns.catplot("passage", "frac_and_weight", data=data_filter, hue="Mutation", order=passage_order, palette="tab20"
                        ,kind="point", dodge=True, hue_order=transition_order, join=False, estimator=weighted_varaint,
                     orient="v")
    g2.set_axis_labels("Passage", "Variant Frequency")
    g2.set(yscale='log')
    g2.set(ylim=(10 ** -6, 10 ** -2))
    # g2.set_xticklabels(fontsize=10, rotation=45)
    g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/Prgress reports/20200913 Final report/plots" +
                      "/Transition_Mutations_point_plot_OPV", dpi=300)
    g2.savefig(output_dir + "/Transition_Mutations_point_plot_OPV", dpi=300)
    plt.close()
    data_filter["passage"] = data_filter["passage"].astype(int)
    # data_pmsc = data_filter[data_filter["Type"] == "Premature Stop Codon"]
    # data_pmsc["mutation_type"] = data_pmsc.Mutation.str.contains("A>G") | data_pmsc.Mutation.str.contains("U>C") | \
    #                              data_pmsc.Mutation.str.contains("C>U") | data_pmsc.Mutation.str.contains("G>A")
    # data_pmsc_transition = data_pmsc[data_pmsc["mutation_type"] == True]
    # g3 = sns.catplot("label", "frac_and_weight", data=data_pmsc_transition, hue="Mutation", order=label_order,
    #                  palette="tab20", estimator=weighted_varaint, orient="v", dodge=True, kind="point", col="Type",
    #                  join=False)
    # g3.set_axis_labels("Sample", "Variant Frequency")
    # g3.set(yscale='log')
    # g3.axes.flat[0].set_yscale('symlog', linthreshy=10**-5)
    # g3.set_xticklabels(rotation=45)
    # # plt.show()
    # g3.savefig(output_dir + "/Transitions_PMSC_Mutations_point_plot", dpi=300)
    # plt.close()
    #
    # g4 = sns.relplot("passage", "frac_and_weight", data=data_filter, hue="Mutation", palette="tab20",
    #                     hue_order=transition_order, estimator=weighted_varaint, col="Type", kind="line",
    #                  col_order=type_order)
    #
    # g4.axes.flat[0].set_yscale('symlog', linthreshy=10**-5)
    # g4.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g4.savefig(output_dir + "/Time_Transition_Mutations_line_plot", dpi=300)
    # plt.close()

    # # A>G Prev Context
    #
    data_filter_ag = data_filter[data_filter["Mutation"] == "A>G"]

    data_filter_ag['Prev'].replace('AA', 'ApA', inplace=True)
    data_filter_ag['Prev'].replace('UA', 'UpA', inplace=True)
    data_filter_ag['Prev'].replace('CA', 'CpA', inplace=True)
    data_filter_ag['Prev'].replace('GA', 'GpA', inplace=True)

    context_order = ["UpA", "ApA", "CpA", "GpA"]
    type_order = ["Synonymous", "Non-Synonymous"]
    data_filter_ag["ADAR_like"] = data_filter_ag.Prev.str.contains('UpA') | data_filter_ag.Prev.str.contains('ApA')
    g5 = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=label_order, palette="tab20",
                kind = "point", dodge=True, hue_order=[True, False], estimator=weighted_varaint, orient="v",
                     col="Type", join=False, col_order=type_order)
    g5.set_axis_labels("", "Variant Frequency")
    g5.set(yscale='log')
    g5.set(ylim=(7*10**-7, 4*10**-3))
    g5.set_xticklabels(rotation=45)
    # plt.show()
    g5.savefig(output_dir + "/Context_point_plot", dpi=300)
    plt.close()
    
    data_filter_ag_pass5 = data_filter_ag.loc[data_filter_ag.passage == 5]
    data_filter_ag_pass5 = data_filter_ag_pass5.loc[data_filter_ag_pass5.Type == "Synonymous"]
    print(data_filter_ag_pass5.to_string())

    ax = sns.boxplot("ADAR_like", "Frequency", data=data_filter_ag_pass5, palette=flatui, order=[True, False])
    ax = sns.stripplot("ADAR_like", "Frequency", data=data_filter_ag_pass5, color=".2", order=[True, False])
    old_statannot.add_stat_annotation(ax, data=data_filter_ag_pass5, x="ADAR_like", y="Frequency",
                        boxPairList=[(True, False)], test='Mann-Whitney', textFormat='star', loc='inside', verbose=1)
    ax.set_yscale('log')
    ax.set_xlabel("ADAR-like\nContext")
    ax.set_ylabel("Variant Frequency")
    ax.set(ylim=(10 ** -5, 10 ** -2))
    sns.despine()
    # plt.tight_layout()
    plt.savefig(output_dir + "/context_p5_point_plot_OPV.png", dpi=300)
    plt.close()
    #
    # print(data_filter_ag.to_string())
    # data_codons = data_filter_ag[data_filter_ag["Prob"] > 0.95]
    # data_codons = data_codons[data_codons["ADAR_like"] == True]
    # data_codons = data_codons[data_codons["Type"] == "Synonymous"]
    # g_codons = sns.catplot("label", "frac_and_weight", data=data_codons, hue="Consensus>Mutated_codon",
    #                           order=label_order, palette="tab20", kind="point",
    #                           join=False, estimator=weighted_varaint, orient="v", dodge=True)
    # g_codons.set_axis_labels("", "Variant Frequency")
    # g_codons.set_xticklabels(fontsize=5, rotation=45)
    # g_codons.set(yscale='log')
    # # g_codons.set(ylim=(10**-7, 10**-3))
    #
    # # plt.show()
    # g_codons.savefig(output_dir + "/codons_point_plot", dpi=300)
    # plt.close()
    #
    # g6 = sns.relplot("passage", "frac_and_weight", data=data_filter_ag, hue="Prev", palette="tab20",
    #                     hue_order=context_order, estimator=weighted_varaint, col="Type", kind="line",
    #                  col_order=type_order)
    #
    # g6.axes.flat[0].set_yscale('symlog', linthreshy=10**-5)
    # g6.set(ylim=(0, 10**-2))
    # yaxis = plt.gca().yaxis
    # yaxis.set_minor_locator(MinorSymLogLocator(1e-1))
    # g6.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g6.savefig(output_dir + "/Context_miseq_sample_time_line_plot", dpi=300)
    # plt.close()
    #
    # data_filter_ag.to_csv(input_dir + "/data_mutation_AG_trajectories.csv", sep=',', encoding='utf-8')
    #
    # g7 = sns.relplot(x="passage", y="frac_and_weight", hue="Prev", data=data_filter_ag, palette="Paired", kind="line",
    #                 style="Type", style_order=type_order, hue_order=context_order, estimator=weighted_varaint)
    # g7.set(yscale="log")
    # g7.fig.suptitle("A>G Mutation trajectories in OPV", y=0.99)
    # g7.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g7.savefig(output_dir + "/ADAR_like_AG_Mutation_Context_trajectories_Context_.png", dpi=300)
    # plt.close()
    #
    #
    #
    # g8 = sns.relplot(x="passage", y="frac_and_weight", hue="ADAR_like", data=data_filter_ag, palette="Paired",
    #                  kind="line", style="Type", style_order=type_order, estimator=weighted_varaint)#, hue_order=context_order)
    # g8.set(yscale="log")
    # g8.fig.suptitle("A>G Mutation trajectories in OPV", y=0.99)
    # g8.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g8.savefig(output_dir + "/ADAR_like_AG_Mutation_ADAR_trajectories_Context.png", dpi=300)
    # plt.close()
    #
    # data_filter_uc = data_filter[data_filter["Mutation"] == "U>C"]
    #
    # data_filter_uc['Next'].replace('UA', 'UpA', inplace=True)
    # data_filter_uc['Next'].replace('UU', 'UpU', inplace=True)
    # data_filter_uc['Next'].replace('UC', 'UpC', inplace=True)
    # data_filter_uc['Next'].replace('UG', 'UpG', inplace=True)
    # data_filter_uc = data_filter_uc[data_filter_uc["passage"] != 0]
    #
    # data_filter_uc.to_csv(input_dir + "/data_filter_uc.csv", sep=',', encoding='utf-8')
    # context_order_uc = ["UpA", "UpU", "UpG",  "UpC"]
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
    #                  hue_order=context_order_uc, estimator=weighted_varaint, col="Type", kind="line",
    #                  col_order=type_order)
    #
    # # g8.set(yscale="log")
    # g10.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -4)
    # # g8.set(ylim=(0, 10 ** -6))
    # # yaxis = plt.gca().yaxis
    # # yaxis.set_minor_locator(MinorSymLogLocator(1e-1))
    # g10.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g10.savefig(output_dir + "/UC_Context_miseq_sample_time_line_plot", dpi=300)
    # plt.close()
    #
    #
    # data_filter_ag_grouped = data_filter_ag.groupby(["ADAR_like", "label", "Type"])["frac_and_weight"].agg(lambda x: weighted_varaint(x))
    # data_filter_ag_grouped = data_filter_ag_grouped.reset_index()
    # data_filter_ag_grouped["passage"] = data_filter_ag_grouped["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    # print(data_filter_ag_grouped.to_string())
    #
    # # data_filter_ag_grouped["frac_and_weight"] = np.log10(data_filter_ag_grouped["frac_and_weight"])
    #
    # data_filter_ag_grouped = data_filter_ag_grouped.rename(columns={"frac_and_weight": "Frequency"})
    # data_filter_ag_grouped["Frequency"] = data_filter_ag_grouped["Frequency"].astype(float)
    # data_filter_ag_grouped["passage"] = np.where(data_filter_ag_grouped["passage"] == "RNA Control", 0, data_filter_ag_grouped["passage"])
    # data_filter_ag_grouped["passage"] = data_filter_ag_grouped["passage"].astype(int)
    #
    #
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
    #                                                                     data_reg_adar_syn['Frequency'])
    # print("stats - ADAR SYN slope: %s" % stat_slope1)
    # print("stats - ADAR SYN intercept: %s" % stat_intercept1)
    #
    # data_reg_adar_nonsyn
    #
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
    #                  line_kws={'label':"y={0:.3g}x+{1:.3g}".format(slope1, intercept1)})
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
    # axes[0, 0].set_xlim(0, 14)
    # axes[0, 0].set_ylim(10**-6, 10**-2)
    # axes[0, 1].legend()
    # axes[0, 1].legend(loc=2)
    # axes[0, 1].set_yscale("log")
    # axes[0, 1].set_xlabel('')
    # axes[0, 1].set_ylabel('')
    # axes[0, 1].set_xlim(0, 14)
    # axes[0, 1].set_ylim(10 ** -6, 10 ** -2)
    # axes[1, 0].legend()
    # axes[1, 0].legend(loc=2)
    # axes[1, 0].set_yscale("log")
    # axes[1, 0].set_xlabel('')
    # axes[1, 0].set_ylabel('')
    # axes[1, 0].set_xlim(0, 14)
    # axes[1, 0].set_ylim(10 ** -6, 10 ** -2)
    # axes[1, 1].legend()
    # axes[1, 1].legend(loc=2)
    # axes[1, 1].set_yscale("log")
    # axes[1, 1].set_xlabel('')
    # axes[1, 1].set_ylabel('')
    # axes[1, 1].set_xlim(0, 14)
    # axes[1, 1].set_ylim(10 ** -6, 10 ** -2)
    # #
    # # for ax in axes:
    # #     ax.set_yscale("log")
    # #     ax.set_ylim(10**-4,10**-2)
    # #     # # l = ax.get_ylabel()
    # #     # # ax.set_ylabel(l, fontsize=8)
    # #     ax.set_xlabel('')
    # #     ax.set_ylabel('')
    # fig.text(0.5, 0.01, 'Passage', ha='center')
    # fig.text(0.001, 0.5, 'Frequency', va='center', rotation='vertical')
    # # fig.suptitle("Mutation trajectories", fontsize=14)
    # fig.tight_layout()
    # # plt.show()
    # plt.savefig(output_dir + "/regplot_AG_Mutation_Context_trajectories.png", dpi=300)
    # plt.close()
    #
    #
    #
    # X = data_reg_adar_nonsyn.passage.values.reshape(-1, 1)
    # y = data_reg_adar_nonsyn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # print("ADAR NON-SYN slope: %s" % str(regr.coef_[0]))
    # print("ADAR NON-SYN intercept: %s" % str(regr.intercept_))
    #
    # g11 = sns.lmplot(x="passage", y="Frequency", data=data_filter_ag_grouped, hue="ADAR_like", markers=["o", "x"],
    #                   hue_order=[True, False], fit_reg=True, col="Type", col_order=type_order)
    # g11.set(xlim=(0, 13))
    # g11.set(yscale="log")
    # g11.set(ylim=(10**-6, 10**-2))
    # # props = dict(boxstyle='round', alpha=0.5, color=sns.color_palette()[0])
    # # textstr = "ADAR-like SYN: y={0:.3g}x+{1:.3g}\nnon ADAR-like SYN: y={2:.3g}x+{3:.3g}".format(slope1, intercept1, slope2, intercept2)
    # # g11.ax.text(0.7, 0.9, textstr, transform=g11.axes.transAxes, fontsize=14, bbox=props)
    # # plt.show()
    # g11.savefig(output_dir + "/lmplot_ADAR_Context", dpi=300)


if __name__ == "__main__":
    main()

