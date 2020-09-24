
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
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/"
    output_dir = input_dir + "/20200924_plots"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)


    data_mutations = pd.read_csv(input_dir + "q38_data_mutation.csv")

    columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
               "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts", "Consensus>Mutated_codon"]
    data_filter = pd.DataFrame(data_mutations, columns=columns)
    data_filter["pval"] = data_filter["pval"].fillna(1)
    data_filter["no_variants"] = data_filter["Frequency"] * data_filter["Read_count"]
    # filter based on pval<0.01 and Prob>0.95
    data_filter["no_variants"] = np.where(data_filter["pval"] > 0.01, 0, data_filter["no_variants"])
    data_filter["no_variants"] = np.where(data_filter["Prob"] < 0.95, 0, data_filter["no_variants"])
    region_lst = [629, 835, 1621, 2329, 3196, 3634, 3925, 4915, 5170, 5239, 5785, 7165]
    data_filter = add_Protein_to_pd_df.add_Protein_to_pd_df_func(data_filter, region_lst)
    data_filter["label"] = np.where(data_filter["label"] == "RVB14-RNA Control", "RVB14\nRNA Control", data_filter["label"])

    data_filter["frac_and_weight"] = list(zip(data_filter.no_variants, data_filter.Read_count))
    #p5-p12
    # data_filter = data_filter[data_filter["label"] != "RVB14-RNA Control"]
    # data_filter = data_filter[data_filter["label"] != "RVB14-p2"]

    data_miseq = data_filter[data_filter["label"] != "RVB14-Next-RNA Control"]
    data_miseq = data_miseq[data_filter["label"] != "RVB14-p1"]
    data_miseq["passage"] = data_miseq["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    data_miseq["passage"] = np.where(data_miseq["passage"] == "RNA Control", 0, data_miseq["passage"])
    data_miseq["passage"] = data_miseq["passage"].astype(int)
    data_miseq["Type"] = data_miseq["Type"].fillna("NonCodingRegion")
    data_miseq.to_csv(input_dir + "/data_miseq.csv", sep=',', encoding='utf-8')

    label_order = ["RVB14\nRNA Control", "RVB14-p2", "RVB14-p5", "RVB14-p7", "RVB14-p8", "RVB14-p10", "RVB14-p12"] # , #"RVB14-Next-RNA Control", "RVB14-p1",
    miseq_order = ["RVB14\nRNA Control", "RVB14-p2", "RVB14-p5", "RVB14-p7", "RVB14-p8", "RVB14-p10", "RVB14-p12"] #
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U"]
    transition_order = ["A>G", "U>C", "C>U", "G>A"]
    type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]

    # A>G Prev Context
    data_filter_ag = data_filter[data_filter["Mutation"] == "A>G"]
    data_filter_ag = data_filter_ag.rename(columns={"Prev": "Context"})

    data_filter_ag['Context'].replace('AA', 'ApA', inplace=True)
    data_filter_ag['Context'].replace('UA', 'UpA', inplace=True)
    data_filter_ag['Context'].replace('CA', 'CpA', inplace=True)
    data_filter_ag['Context'].replace('GA', 'GpA', inplace=True)

    context_order = ["UpA", "ApA", "CpA", "GpA"]
    type_order = ["Synonymous", "Non-Synonymous"]

    data_miseq_ag = data_miseq[data_miseq["Mutation"] == "A>G"]
    data_miseq_ag = data_miseq_ag.rename(columns={"Prev": "Context"})

    data_miseq_ag['Context'].replace('AA', 'ApA', inplace=True)
    data_miseq_ag['Context'].replace('UA', 'UpA', inplace=True)
    data_miseq_ag['Context'].replace('CA', 'CpA', inplace=True)
    data_miseq_ag['Context'].replace('GA', 'GpA', inplace=True)


    data_miseq_ag["ADAR_like"] = data_miseq_ag.Context.str.contains('UpA') | data_miseq_ag.Context.str.contains('ApA')
    print(data_miseq_ag.to_string())
    data_miseq_ag.to_csv(input_dir + "/data_mutation_AG_trajectories.csv", sep=',', encoding='utf-8')
    data_miseq_ag_pass5 = data_miseq_ag.loc[data_miseq_ag.label == "RVB14-p5"]
    data_miseq_ag_pass5 = data_miseq_ag_pass5[data_miseq_ag_pass5["pval"] < 0.01]
    # data_miseq_ag_pass5 = data_miseq_ag_pass5.loc[data_miseq_ag_pass5.Type == "Synonymous"]

    # U>C Prev Context
    data_miseq_uc = data_miseq[data_miseq["Mutation"] == "U>C"]

    data_miseq_uc['Next'].replace('UA', 'UpA', inplace=True)
    data_miseq_uc['Next'].replace('UU', 'UpU', inplace=True)
    data_miseq_uc['Next'].replace('UC', 'UpC', inplace=True)
    data_miseq_uc['Next'].replace('UG', 'UpG', inplace=True)
    data_miseq_uc = data_miseq_uc[data_miseq_uc["passage"] != 0]

    data_miseq_uc.to_csv(input_dir + "/data_miseq_uc.csv", sep=',', encoding='utf-8')
    context_order_uc = ["UpA", "UpU", "UpG",  "UpC"]

    #Plots

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

    g2 = sns.catplot("label", "frac_and_weight", data=data_filter, hue="Mutation", order=label_order, palette="tab20"
                        ,kind="point", dodge=True, hue_order=transition_order, join=False, estimator=weighted_varaint,
                     orient="v", legend=False)
    g2.set_axis_labels("", "Variant Frequency")
    g2.set(yscale='log')
    g2.set(ylim=(10 ** -6, 10 ** -3))
    # g2.set_yticklabels(fontsize=12)
    g2.set_xticklabels(fontsize=10, rotation=45)
    # plt.show()
    g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/Transition_Mutations_point_plot_RV", dpi=300)
    # g2.savefig(output_dir + "/Transition_Mutations_point_plot", dpi=300)
    plt.close()

    data_pmsc = data_miseq[data_miseq["Type"] == "Premature Stop Codon"]
    data_pmsc["mutation_type"] = data_pmsc.Mutation.str.contains("A>G") | data_pmsc.Mutation.str.contains("U>C") | \
                                 data_pmsc.Mutation.str.contains("C>U") | data_pmsc.Mutation.str.contains("G>A")
    data_pmsc_transition = data_pmsc[data_pmsc["mutation_type"] == True]
    g3 = sns.catplot("label", "frac_and_weight", data=data_pmsc_transition, hue="Mutation", order=miseq_order,
                     palette="tab20", estimator=weighted_varaint, orient="v", dodge=True, kind="point", col="Type",
                     join=False)
    g3.set_axis_labels("Sample", "Variant Frequency")
    g3.set(yscale='log')
    g3.axes.flat[0].set_yscale('symlog', linthreshy=10**-5)
    g3.set_xticklabels(rotation=45)
    # plt.show()
    g3.savefig(output_dir + "/Transitions_PMSC_Mutations_point_plot", dpi=300)
    plt.close()

    g4 = sns.relplot("passage", "frac_and_weight", data=data_miseq, hue="Mutation", palette="tab20",
                        hue_order=transition_order, estimator=weighted_varaint, col="Type", kind="line",
                     col_order=type_order)

    g4.axes.flat[0].set_yscale('symlog', linthreshy=10**-5)
    g4.set_axis_labels("Passage", "Variant Frequency")
    # plt.show()
    g4.savefig(output_dir + "/Time_Transition_Mutations_line_plot", dpi=300)
    plt.close()

    # A>G Prev Context
    g5 = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="Context", order=label_order, palette="Set2",
                kind = "point", dodge=True, hue_order=context_order, estimator=weighted_varaint, orient="v",
                     col="Type", join=False, col_order=type_order)
    g5.set_axis_labels("", "Variant Frequency")
    g5.set(yscale='log')
    g5.set(ylim=(7*10**-7, 4*10**-3))
    g5.set_xticklabels(rotation=45)
    # plt.show()
    g5.savefig(output_dir + "/Context_point_plot", dpi=300)
    plt.close()

    g_context = sns.catplot("ADAR_like", "frac_and_weight", data=data_miseq_ag_pass5,
                              order=[True, False], palette="Set2", kind="point",
                              join=False, estimator=weighted_varaint, orient="v", dodge=True)
    # add_stat_annotation(g_context, data=data_miseq_ag_pass5, x="ADAR_like", y="Frequency", order=[True, False],
    #                     boxPairList=[(True, False)], test='Mann-Whitney', textFormat='star', loc='inside', verbose=2)
    g_context.set_axis_labels("ADAR-like\nContext", "Variant Frequency")
    g_context.set_xticklabels(rotation=45)
    g_context.set(yscale='log')
    g_context.set(ylim=(5*10**-5, 5*10**-3))
    plt.tight_layout()
    plt.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/context_p5_point_plot.png", dpi=300)
    plt.savefig(output_dir + "/context_p5_point_plot.png", dpi=300)
    plt.close()

    flatui = ["#9b59b6", "#3498db"]
    ax = sns.boxplot("ADAR_like", "Frequency", data=data_miseq_ag_pass5, palette=flatui)
    ax = sns.stripplot("ADAR_like", "Frequency", data=data_miseq_ag_pass5, color=".2")
    add_stat_annotation(ax, data=data_miseq_ag_pass5, x="ADAR_like", y="Frequency",
                        boxPairList=[(True, False)], test='Mann-Whitney', textFormat='star', loc='inside', verbose=2)
    ax.set_yscale('log')
    ax.set_xlabel("ADAR-like\nContext")
    ax.set_ylabel("Variant Frequency")
    ax.set(ylim=(10**-4, 10**-2))
    sns.despine()
    plt.tight_layout()
    plt.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/context_p5_point_plot_v2.png", dpi=300)
    plt.savefig(output_dir + "/context_p5_point_plot_v2.png", dpi=300)
    plt.close()

    data_codons = data_miseq_ag[data_miseq_ag["pval"] < 0.01]
    data_codons = data_codons[data_codons["ADAR_like"] == True]
    data_codons = data_codons[data_codons["Type"] == "Synonymous"]
    g_codons = sns.catplot("label", "frac_and_weight", data=data_codons, hue="Consensus>Mutated_codon",
                              order=miseq_order, palette="tab20", kind="point",
                              join=False, estimator=weighted_varaint, orient="v", dodge=True)
    g_codons.set_axis_labels("", "Variant Frequency")
    g_codons.set_xticklabels(fontsize=5, rotation=45)
    g_codons.set(yscale='log')
    # g_codons.set(ylim=(10**-7, 10**-3))

    # plt.show()
    g_codons.savefig(output_dir + "/codons_point_plot", dpi=300)
    plt.close()

    g_codons_3c = sns.catplot("label", "frac_and_weight", data=data_codons[data_codons["Protein"]=="3C"], hue="Consensus>Mutated_codon",
                              order=miseq_order, palette="tab20", kind="point",
                              join=False, estimator=weighted_varaint, orient="v", dodge=True)
    g_codons_3c.set_axis_labels("", "Variant Frequency")
    g_codons_3c.set_xticklabels(fontsize=5, rotation=45)
    g_codons_3c.set(yscale='log')
    # g_codons.set(ylim=(10**-7, 10**-3))

    # plt.show()
    g_codons_3c.savefig(output_dir + "/codons_point_plot_3c", dpi=300)
    plt.close()

    g6 = sns.relplot("passage", "frac_and_weight", data=data_miseq_ag, hue="Context", palette="tab20",
                        hue_order=context_order, estimator=weighted_varaint, col="Type", kind="line",
                     col_order=type_order)

    g6.axes.flat[0].set_yscale('symlog', linthreshy=10**-5)
    g6.set(ylim=(0, 10**-2))
    yaxis = plt.gca().yaxis
    yaxis.set_minor_locator(fits_new_plotter.MinorSymLogLocator(1e-1))
    g6.set_axis_labels("Passage", "Variant Frequency")
    # plt.show()
    g6.savefig(output_dir + "/Context_miseq_sample_time_line_plot", dpi=300)
    plt.close()

    g7 = sns.relplot(x="passage", y="frac_and_weight", hue="Context", data=data_miseq_ag, palette="Paired", kind="line",
                    style="Type", style_order=type_order, hue_order=context_order, estimator=weighted_varaint)
    g7.set(yscale="log")
    g7.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    g7.set_axis_labels("Passage", "Variant Frequency")
    # plt.show()
    g7.savefig(output_dir + "/ADAR_like_AG_Mutation_Context_trajectories_Context_.png", dpi=300)
    plt.close()

    g8 = sns.relplot(x="passage", y="frac_and_weight", hue="ADAR_like", data=data_miseq_ag, palette="Paired",
                     kind="line", style="Type", style_order=type_order, estimator=weighted_varaint)#, hue_order=context_order)
    g8.set(yscale="log")
    g8.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    g8.set_axis_labels("Passage", "Variant Frequency")
    # plt.show()
    g8.savefig(output_dir + "/ADAR_like_AG_Mutation_ADAR_trajectories_Context.png", dpi=300)
    plt.close()

    g9 = sns.catplot("label", "frac_and_weight", data=data_miseq_uc, hue="Next", order=miseq_order, palette="tab20",
                     hue_order=context_order_uc, estimator=weighted_varaint, orient="v", dodge=True, kind="point",
                     col="Type", join=False, col_order=type_order)
    g9.set_axis_labels("", "Variant Frequency")
    g9.set(yscale='log')
    g9.set(ylim=(10 ** -7, 10 ** -3))
    g9.set_xticklabels(rotation=45)
    # plt.show()
    g9.savefig(output_dir + "/UC_Context_point_plot", dpi=300)
    plt.close()

    g10 = sns.relplot("passage", "frac_and_weight", data=data_miseq_uc, hue="Next", palette="tab20",
                     hue_order=context_order_uc, estimator=weighted_varaint, col="Type", kind="line",
                     col_order=type_order)

    # g8.set(yscale="log")
    g10.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -4)
    # g8.set(ylim=(0, 10 ** -6))
    # yaxis = plt.gca().yaxis
    yaxis.set_minor_locator(fits_new_plotter.MinorSymLogLocator(1e-1))
    g10.set_axis_labels("Passage", "Variant Frequency")
    # plt.show()
    g10.savefig(output_dir + "/UC_Context_miseq_sample_time_line_plot", dpi=300)
    plt.close()


    data_miseq_ag_grouped = data_miseq_ag.groupby(["ADAR_like", "label", "Type"])["frac_and_weight"].agg(lambda x: weighted_varaint(x))
    data_miseq_ag_grouped = data_miseq_ag_grouped.reset_index()
    data_miseq_ag_grouped["passage"] = data_miseq_ag_grouped["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    print(data_miseq_ag_grouped.to_string())

    # data_miseq_ag_grouped["frac_and_weight"] = np.log10(data_miseq_ag_grouped["frac_and_weight"])

    data_miseq_ag_grouped = data_miseq_ag_grouped.rename(columns={"frac_and_weight": "Frequency"})
    data_miseq_ag_grouped["Frequency"] = data_miseq_ag_grouped["Frequency"].astype(float)
    data_miseq_ag_grouped["passage"] = np.where(data_miseq_ag_grouped["passage"] == "RNA Control", 0, data_miseq_ag_grouped["passage"])
    data_miseq_ag_grouped["passage"] = data_miseq_ag_grouped["passage"].astype(int)

    data_reg_adar = data_miseq_ag_grouped[data_miseq_ag_grouped["ADAR_like"] == True]
    data_reg_nonadar = data_miseq_ag_grouped[data_miseq_ag_grouped["ADAR_like"] == False]

    data_reg_adar_syn = data_reg_adar[data_reg_adar["Type"] == "Synonymous"]
    data_reg_nonadar_syn = data_reg_nonadar[data_reg_nonadar["Type"] == "Synonymous"]

    data_reg_adar_nonsyn = data_reg_adar[data_reg_adar["Type"] == "Non-Synonymous"]
    data_reg_nonadar_non_syn = data_reg_nonadar[data_reg_nonadar["Type"] == "Non-Synonymous"]

    fig, axes = plt.subplots(2, 2, sharey=True, sharex=True)
    stat_slope1, stat_intercept1, r_value1, p_value1, std_err1 = stats.linregress(data_reg_adar_syn['passage'],
                                                                        data_reg_adar_syn['Frequency'])
    print("stats - ADAR SYN slope: %s" % stat_slope1)
    print("stats - ADAR SYN intercept: %s" % stat_intercept1)

    # data_reg_adar_nonsyn
    stat_slope2, stat_intercept2, r_value1, p_value1, std_err1 = stats.linregress(data_reg_adar_nonsyn['passage'],
                                                                                  data_reg_adar_nonsyn['Frequency'])

    print("stats - ADAR NON-SYN slope: %s" % stat_slope2)
    print("stats - ADAR NON-SYN intercept: %s" % stat_intercept2)

    from sklearn import linear_model
    regr = linear_model.LinearRegression()
    X = data_reg_adar_syn.passage.values.reshape(-1, 1)
    y = data_reg_adar_syn.Frequency.values.reshape(-1, 1)
    regr.fit(X, y)
    slope1 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    intercept1 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    print("ADAR SYN slope: %s" % str(regr.coef_[0]).split("[")[1].split("]")[0])
    print("ADAR SYN intercept: %s" % str(regr.intercept_).split("[")[1].split("]")[0])

    g1 = sns.regplot(x="passage", y="Frequency", data=data_reg_adar_syn, ax=axes[0, 0],
                     line_kws={'label':"y={0:.3g}x+{1:.3g}".format(slope1, intercept1)})
    g1.set(title="ADAR-like Synonymous")

    X = data_reg_nonadar_syn.passage.values.reshape(-1, 1)
    y = data_reg_nonadar_syn.Frequency.values.reshape(-1, 1)
    regr.fit(X, y)
    slope2 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    intercept2 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    print("NON-ADAR SYN slope: %s" % str(regr.coef_[0]).split("[")[1].split("]")[0])
    print("NON-ADAR SYN intercept: %s" % str(regr.intercept_).split("[")[1].split("]")[0])


    g2 = sns.regplot(x="passage", y="Frequency", data=data_reg_nonadar_syn, ax=axes[0, 1],
                     line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope2, intercept2)})
    g2.set(title="Non ADAR-like Synonymous")

    X = data_reg_adar_nonsyn.passage.values.reshape(-1, 1)
    y = data_reg_adar_nonsyn.Frequency.values.reshape(-1, 1)
    regr.fit(X, y)
    slope3 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    intercept3 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    print("ADAR NON-SYN slope: %s" % str(regr.coef_[0]).split("[")[1].split("]")[0])
    print("ADAR NON-SYN intercept: %s" % str(regr.intercept_).split("[")[1].split("]")[0])

    g3 = sns.regplot(x="passage", y="Frequency", data=data_reg_adar_nonsyn, ax=axes[1, 0],
                     line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope3, intercept3)})
    g3.set(title="ADAR-like Non Synonymous")

    X = data_reg_nonadar_non_syn.passage.values.reshape(-1, 1)
    y = data_reg_nonadar_non_syn.Frequency.values.reshape(-1, 1)
    regr.fit(X, y)
    slope4 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    intercept4 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    print("NON-ADAR NON-SYN slope: %s" % str(regr.coef_[0]).split("[")[1].split("]")[0])
    print("NON-ADAR NON-SYN intercept: %s" % str(regr.intercept_).split("[")[1].split("]")[0])

    g4 = sns.regplot(x="passage", y="Frequency", data=data_reg_nonadar_non_syn, ax=axes[1, 1],
                     line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope4, intercept4)})
    g4.set(title="Non ADAR-like Non Synonymous")

    axes[0, 0].legend()
    axes[0, 0].legend(loc=2)
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_xlabel('')
    axes[0, 0].set_ylabel('')
    axes[0, 0].set_xlim(0, 14, 2)
    axes[0, 0].set_ylim(10**-6, 10**-2)
    axes[0, 1].legend()
    axes[0, 1].legend(loc=2)
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_xlabel('')
    axes[0, 1].set_ylabel('')
    axes[0, 1].set_xlim(0, 14, 2)
    axes[0, 1].set_ylim(10 ** -6, 10 ** -2)
    axes[1, 0].legend()
    axes[1, 0].legend(loc=2)
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_xlabel('')
    axes[1, 0].set_ylabel('')
    axes[1, 0].set_xlim(0, 14, 2)
    axes[1, 0].set_ylim(10 ** -6, 10 ** -2)
    axes[1, 1].legend()
    axes[1, 1].legend(loc=2)
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xlabel('')
    axes[1, 1].set_ylabel('')
    axes[1, 1].set_xlim(0, 14, 2)
    axes[1, 1].set_ylim(10 ** -6, 10 ** -2)
    #
    # for ax in axes:
    #     ax.set_yscale("log")
    #     ax.set_ylim(10**-4,10**-2)
    #     # # l = ax.get_ylabel()
    #     # # ax.set_ylabel(l, fontsize=8)
    #     ax.set_xlabel('')
    #     ax.set_ylabel('')
    fig.text(0.5, 0.01, 'Passage', ha='center')
    fig.text(0.001, 0.5, 'Frequency', va='center', rotation='vertical')
    # fig.suptitle("Mutation trajectories", fontsize=14)
    fig.tight_layout()
    # plt.show()
    plt.savefig(output_dir + "/regplot_AG_Mutation_Context_trajectories.png", dpi=300)
    plt.close()

    X = data_reg_adar_nonsyn.passage.values.reshape(-1, 1)
    y = data_reg_adar_nonsyn.Frequency.values.reshape(-1, 1)
    regr.fit(X, y)
    print("ADAR NON-SYN slope: %s" % str(regr.coef_[0]))
    print("ADAR NON-SYN intercept: %s" % str(regr.intercept_))

    g11 = sns.lmplot(x="passage", y="Frequency", data=data_miseq_ag_grouped, hue="ADAR_like", markers=["o", "x"],
                      hue_order=[True, False], fit_reg=True, col="Type", col_order=type_order)
    g11.set(xlim=(0, 13))
    g11.set(yscale="log")
    g11.set(ylim=(10**-6, 10**-2))
    # props = dict(boxstyle='round', alpha=0.5, color=sns.color_palette()[0])
    # textstr = "ADAR-like SYN: y={0:.3g}x+{1:.3g}\nnon ADAR-like SYN: y={2:.3g}x+{3:.3g}".format(slope1, intercept1, slope2, intercept2)
    # g11.ax.text(0.7, 0.9, textstr, transform=g11.axes.transAxes, fontsize=14, bbox=props)
    # plt.show()
    g11.savefig(output_dir + "/lmplot_ADAR_Context", dpi=300)


if __name__ == "__main__":
    main()

