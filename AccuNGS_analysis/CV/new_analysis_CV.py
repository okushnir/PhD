
import pandas as pd
import os
import matplotlib.pyplot as plt
from FITS_analysis import fits_plotter
import numpy as np
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker
import seaborn as sns
from FITS_analysis import fits_new_plotter
from AccuNGS_analysis.adar_mutation_palette import mutation_palette
from AccuNGS_analysis.Linear_regression import linear_reg
from scipy import stats

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
    date = "20201130"
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/CVB3"
    output_dir = input_dir + "/plots_q38_filtered/%s" % date
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)


    data_mutations = pd.read_csv(input_dir + "/q38_data_mutation.csv")

    columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
               "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts"]
    data_filter = pd.DataFrame(data_mutations, columns=columns)
    data_filter["pval"] = data_filter["pval"].fillna(1)
    data_filter["no_variants"] = data_filter["Frequency"] * data_filter["Read_count"]
    """filter based on pval<0.01 and Prob>0.95"""
    # data_filter["no_variants"] = np.where(data_filter["pval"] > 0.01, 0, data_filter["no_variants"])
    # data_filter["no_variants"] = np.where(data_filter["Prob"] < 0.95, 0, data_filter["no_variants"])
    # data_filter["Read_count"] = data_filter[data_filter["Read_count"] > 10000]
    data_filter["label"] = np.where(data_filter["label"] == "CVB3-RNA Control", "CVB3\nRNA Control", data_filter["label"])

    data_filter["frac_and_weight"] = list(zip(data_filter.no_variants, data_filter.Read_count))

    data_filter["passage"] = data_filter["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    data_filter["passage"] = np.where(data_filter["passage"] == "CVB3\nRNA Control", 0, data_filter["passage"])
    data_filter["passage"] = data_filter["passage"].astype(int)
    data_filter["Type"] = data_filter["Type"].fillna("NonCodingRegion")
    # data_filter.to_csv(input_dir + "/data_filter.csv", sep=',', encoding='utf-8')

    label_order = ["CVB3\nRNA Control", "CVB3-p2", "CVB3-p5", "CVB3-p8", "CVB3-p10", "CVB3-p12"]
    passage_order = ["0", "2", "5", "8", "10", "12"]
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U"]
    transition_order = ["A>G", "U>C", "G>A", "C>U"]
    type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]

    g1 = sns.catplot("label", "frac_and_weight", data=data_filter, hue="Mutation", order=label_order, palette="tab20",
                        kind="point", hue_order=mutation_order, join=False, estimator=weighted_varaint, orient="v", dodge=True)
    g1.set_axis_labels("", "Variant Frequency")
    g1.set_xticklabels(fontsize=9, rotation=45)
    g1.set(yscale='log')
    g1.set(ylim=(10 ** -7, 10 ** -3))
    # plt.show()
    g1.savefig(output_dir + "/All_Mutations_point_plot", dpi=300)
    plt.close()
    data_filter["passage"] = data_filter["passage"].astype(str)
    g2 = sns.catplot("passage", "frac_and_weight", data=data_filter, hue="Mutation", order=passage_order, palette=mutation_palette(4),
                        kind="point", hue_order=transition_order, join=False, estimator=weighted_varaint, orient="v",
                     dodge=True, legend=True)
    g2.set_axis_labels("Passage", "Variant Frequency")
    g2.set(yscale='log')
    g2.set(ylim=(10 ** -6, 10 ** -2))
    # g2.set_yticklabels(fontsize=12)
    # g2.set_xticklabels(fontsize=10, rotation=45)
    # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/Prgress reports/20200913 Final report/plots" +
    #                   "/Transition_Mutations_point_plot_CVB3", dpi=300)
    g2.savefig(output_dir + "/Transition_Mutations_point_plot_CVB3", dpi=300)
    # g2.savefig(output_dir + "/Transition_Mutations_point_plot", dpi=300)
    plt.close()
    data_filter["passage"] = data_filter["passage"].astype(int)
    data_filter["replica"] = 1
    linear_reg(data_filter, output_dir, transition_order, type_order, virus="CVB3", replica=1)

    # data_filter_grouped = data_filter.groupby(["label", "passage", "Type", "Mutation"])[
    #     "frac_and_weight"].agg(
    #     lambda x: weighted_varaint(x))
    # data_filter_grouped = data_filter_grouped.reset_index()
    # data_filter_grouped = data_filter_grouped.rename(columns={"frac_and_weight": "Frequency"})
    # data_filter_grouped["Frequency"] = data_filter_grouped["Frequency"].astype(float)
    # data_filter_grouped = data_filter_grouped[data_filter_grouped["label"] != "RNA Control\nPrimer ID"]
    # data_filter_grouped = data_filter_grouped[data_filter_grouped["label"] != "RNA Control\nRND"]
    # # data_filter_grouped = data_filter_grouped[data_filter_grouped["replica"] == replica]
    # print(data_filter_grouped.to_string())
    #
    # data_reg_ag = data_filter_grouped[data_filter_grouped["Mutation"] == "A>G"]
    # data_reg_uc = data_filter_grouped[data_filter_grouped["Mutation"] == "U>C"]
    # data_reg_ga = data_filter_grouped[data_filter_grouped["Mutation"] == "G>A"]
    # data_reg_cu = data_filter_grouped[data_filter_grouped["Mutation"] == "C>U"]
    #
    # data_reg_ag_syn = data_reg_ag[data_reg_ag["Type"] == "Synonymous"]
    # data_reg_uc_syn = data_reg_uc[data_reg_uc["Type"] == "Synonymous"]
    # data_reg_ga_syn = data_reg_ga[data_reg_ga["Type"] == "Synonymous"]
    # data_reg_cu_syn = data_reg_cu[data_reg_cu["Type"] == "Synonymous"]
    #
    # data_reg_ag_non_syn = data_reg_ag[data_reg_ag["Type"] == "Non-Synonymous"]
    # data_reg_uc_non_syn = data_reg_uc[data_reg_uc["Type"] == "Non-Synonymous"]
    # data_reg_ga_non_syn = data_reg_ga[data_reg_ga["Type"] == "Non-Synonymous"]
    # data_reg_cu_non_syn = data_reg_cu[data_reg_cu["Type"] == "Non-Synonymous"]
    #
    # data_reg_ga_pmsc = data_reg_ga[data_reg_ga["Type"] == "Premature Stop Codon"]
    # data_reg_cu_pmsc = data_reg_cu[data_reg_cu["Type"] == "Premature Stop Codon"]
    #
    # stat_slope1, stat_intercept1, r_value1, p_value1, std_err1 = stats.linregress(data_reg_ag_syn['passage'],
    #                                                                               data_reg_ag_syn
    #                                                                               ['Frequency'])
    # stat_slope2, stat_intercept2, r_value2, p_value2, std_err2 = stats.linregress(data_reg_uc_syn['passage'],
    #                                                                               data_reg_uc_syn[
    #                                                                                   'Frequency'])
    # stat_slope3, stat_intercept3, r_value3, p_value3, std_err3 = stats.linregress(data_reg_ga_syn['passage'],
    #                                                                               data_reg_ga_syn[
    #                                                                                   'Frequency'])
    # stat_slope4, stat_intercept4, r_value4, p_value4, std_err4 = stats.linregress(data_reg_cu_syn['passage'],
    #                                                                               data_reg_cu_syn[
    #                                                                                   'Frequency'])
    # # data_reg_adar_nonsyn
    # stat_slope5, stat_intercept5, r_value5, p_value5, std_err5 = stats.linregress(data_reg_ag_non_syn['passage'],
    #                                                                               data_reg_ag_non_syn[
    #                                                                                   'Frequency'])
    # stat_slope6, stat_intercept6, r_value6, p_value6, std_err6 = stats.linregress(data_reg_uc_non_syn['passage'],
    #                                                                               data_reg_uc_non_syn[
    #                                                                                   'Frequency'])
    # stat_slope7, stat_intercept7, r_value7, p_value7, std_err7 = stats.linregress(data_reg_ga_non_syn['passage'],
    #                                                                               data_reg_ga_non_syn[
    #                                                                                   'Frequency'])
    # stat_slope8, stat_intercept8, r_value8, p_value8, std_err8 = stats.linregress(data_reg_cu_non_syn['passage'],
    #                                                                               data_reg_cu_non_syn[
    #                                                                                   'Frequency'])
    # # pmsc
    # stat_slope11, stat_intercept11, r_value11, p_value11, std_err11 = stats.linregress(data_reg_ga_pmsc['passage'],
    #                                                                                    data_reg_ga_pmsc[
    #                                                                                        'Frequency'])
    # stat_slope12, stat_intercept12, r_value12, p_value12, std_err12 = stats.linregress(data_reg_cu_pmsc['passage'],
    #                                                                                    data_reg_cu_pmsc[
    #                                                                                        'Frequency'])
    # data_filter_grouped = data_filter_grouped.rename(columns={"passage": "Passage"})
    # reg_plot = sns.lmplot(x="Passage", y="Frequency", data=data_filter_grouped, hue="Mutation",
    #                       hue_order=transition_order, fit_reg=True, col="Mutation",
    #                       col_order=transition_order, row="Type", row_order=type_order, palette=mutation_palette(4),
    #                       line_kws={'label': "Linear Reg"}, legend=True, height=6)  # markers=["o", "v", "x"]
    # reg_plot.fig.subplots_adjust(wspace=.02)
    # ax = reg_plot.axes[0, 0]
    # ax.legend()
    # leg = ax.get_legend()
    # leg._loc = 2
    # L_labels = leg.get_texts()
    # label_line_1 = "y={0:.3g}x+{1:.3g}".format(stat_slope1, stat_intercept1)
    # label_line_2 = "y={0:.3g}x+{1:.3g}".format(stat_slope2, stat_intercept2)
    # label_line_3 = "y={0:.3g}x+{1:.3g}".format(stat_slope3, stat_intercept3)
    # label_line_4 = "y={0:.3g}x+{1:.3g}".format(stat_slope4, stat_intercept4)
    # label_line_5 = "y={0:.3g}x+{1:.3g}".format(stat_slope5, stat_intercept5)
    # label_line_6 = "y={0:.3g}x+{1:.3g}".format(stat_slope6, stat_intercept6)
    # label_line_7 = "y={0:.3g}x+{1:.3g}".format(stat_slope7, stat_intercept7)
    # label_line_8 = "y={0:.3g}x+{1:.3g}".format(stat_slope8, stat_intercept8)
    # label_line_11 = "y={0:.3g}x+{1:.3g}".format(stat_slope11, stat_intercept11)
    # label_line_12 = "y={0:.3g}x+{1:.3g}".format(stat_slope12, stat_intercept12)
    # L_labels[0].set_text(label_line_1)
    # ax = reg_plot.axes[0, 1]
    # ax.legend()
    # leg = ax.get_legend()
    # leg._loc = 2
    # L_labels = leg.get_texts()
    # L_labels[0].set_text(label_line_2)
    # ax = reg_plot.axes[0, 2]
    # ax.legend()
    # leg = ax.get_legend()
    # leg._loc = 2
    # L_labels = leg.get_texts()
    # L_labels[0].set_text(label_line_3)
    # ax = reg_plot.axes[0, 3]
    # ax.legend()
    # leg = ax.get_legend()
    # leg._loc = 2
    # L_labels = leg.get_texts()
    # L_labels[0].set_text(label_line_4)
    # ax = reg_plot.axes[1, 0]
    # ax.legend()
    # leg = ax.get_legend()
    # leg._loc = 2
    # L_labels = leg.get_texts()
    # L_labels[0].set_text(label_line_5)
    # ax = reg_plot.axes[1, 1]
    # ax.legend()
    # leg = ax.get_legend()
    # leg._loc = 2
    # L_labels = leg.get_texts()
    # L_labels[0].set_text(label_line_6)
    # ax = reg_plot.axes[1, 2]
    # ax.legend()
    # leg = ax.get_legend()
    # leg._loc = 2
    # L_labels = leg.get_texts()
    # L_labels[0].set_text(label_line_7)
    # ax = reg_plot.axes[1, 3]
    # ax.legend()
    # leg = ax.get_legend()
    # leg._loc = 2
    # L_labels = leg.get_texts()
    # L_labels[0].set_text(label_line_8)
    # ax = reg_plot.axes[2, 2]
    # ax.legend()
    # leg = ax.get_legend()
    # leg._loc = 2
    # L_labels = leg.get_texts()
    # L_labels[0].set_text(label_line_11)
    # ax = reg_plot.axes[2, 3]
    # ax.legend()
    # leg = ax.get_legend()
    # leg._loc = 2
    # L_labels = leg.get_texts()
    # L_labels[0].set_text(label_line_12)
    #
    # reg_plot.set(xlim=(0, 13))
    # reg_plot.set(ylim=(0.000, 0.001))
    # # reg_plot.fig.suptitle("RV #%s" % str(replica), y=0.99)
    # # plt.tight_layout()
    # reg_plot.savefig(output_dir + "/transition_lmplot.png", dpi=300)
    # plt.close()
    #
    # columns = ["Mutation", "Type", "Slope", "Intercept"]
    # mutation_rate_df = pd.DataFrame(columns=columns)
    # mutation_rate_df.loc[0] = ["A>G", "Synonymous", stat_slope1, stat_intercept1]
    # mutation_rate_df.loc[1] = ["U>C", "Synonymous", stat_slope2, stat_intercept2]
    # mutation_rate_df.loc[2] = ["G>A", "Synonymous", stat_slope3, stat_intercept3]
    # mutation_rate_df.loc[3] = ["C>U", "Synonymous", stat_slope4, stat_intercept4]
    # mutation_rate_df.loc[4] = ["A>G", "Non-Synonymous", stat_slope5, stat_intercept5]
    # mutation_rate_df.loc[5] = ["U>C", "Non-Synonymous", stat_slope6, stat_intercept6]
    # mutation_rate_df.loc[6] = ["G>A", "Non-Synonymous", stat_slope7, stat_intercept7]
    # mutation_rate_df.loc[7] = ["C>U", "Non-Synonymous", stat_slope8, stat_intercept8]
    # mutation_rate_df.loc[8] = ["G>A",  "Pre Mature Stop Codon", stat_slope11, stat_intercept11]
    # mutation_rate_df.loc[9] = ["C>U",  "Pre Mature Stop Codon", stat_slope12, stat_intercept12]
    # mutation_rate_df["Virus"] = "CVB3"
    # mutation_rate_df["Replica"] = 1
    # mutation_rate_df.to_csv(output_dir + "/mutation_rate.csv", sep=',', encoding='utf-8')

    data_filter["Transition"] = data_filter.Mutation.str.contains("A>G") | data_filter.Mutation.str.contains("U>C") \
                                | data_filter.Mutation.str.contains("C>U") | data_filter.Mutation.str.contains("G>A")
    data_transition = data_filter[data_filter["Transition"] == True]
    position_mutation = sns.relplot(x="Pos", y="Frequency", data=data_transition, hue="Mutation",
                                    col="passage", hue_order=transition_order, col_wrap=3,
                                    palette=mutation_palette(4), height=4)
    #
    # data_pmsc = data_filter[data_filter["Type"] == "Premature Stop Codon"]
    # data_pmsc["mutation_type"] = data_pmsc.Mutation.str.contains("A>G") | data_pmsc.Mutation.str.contains("U>C") | \
    #                              data_pmsc.Mutation.str.contains("C>U") | data_pmsc.Mutation.str.contains("G>A")
    # data_pmsc_transition = data_pmsc[data_pmsc["mutation_type"] == True]
    # g3 = sns.catplot("label", "frac_and_weight", data=data_pmsc_transition, hue="Mutation", order=label_order, palette="tab20",
    #                     estimator=weighted_varaint, orient="v", dodge=True, kind="point",
    #                  col="Type", join=False)
    # g3.set_axis_labels("Sample", "Variant Frequency")
    # g3.set(yscale='log')
    # g3.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -5)
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
    #
    # # A>G Prev Context
    #
    data_filter_ag = data_filter[data_filter["Mutation"] == "A>G"]
    data_filter_ag['Prev'].replace('AA', 'ApA', inplace=True)
    data_filter_ag['Prev'].replace('UA', 'UpA', inplace=True)
    data_filter_ag['Prev'].replace('CA', 'CpA', inplace=True)
    data_filter_ag['Prev'].replace('GA', 'GpA', inplace=True)
    data_filter_ag = data_filter_ag.rename(columns={"Prev": "Context"})
    data_filter_ag["ADAR_like"] = data_filter_ag.Context.str.contains('UpA') | data_filter_ag.Context.str.contains(
        'ApA')
    data_filter_ag.to_pickle(output_dir + "/data_filter_ag.pkl")

    type_order = ["Synonymous", "Non-Synonymous"]

    g5 = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=label_order,
                     palette=mutation_palette(2), kind="point", dodge=True, hue_order=[True, False], estimator=weighted_varaint,
                     orient="v", col="Type", join=False, col_order=type_order)
    g5.set_axis_labels("", "Variant Frequency")
    g5.set(yscale='log')
    g5.set(ylim=(7 * 10 ** -7, 4 * 10 ** -3))
    g5.set_xticklabels(fontsize=9, rotation=90)
    # plt.show()
    g5.savefig(output_dir + "/Context_point_plot", dpi=300)
    plt.close()

    data_filter_ag = data_filter_ag[data_filter_ag["passage"] != 0]
    g6 = sns.relplot("passage", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", palette=mutation_palette(2),
                        hue_order=[True, False], estimator=weighted_varaint, col="Type", kind="line",
                     col_order=type_order)

    g6.axes.flat[0].set_yscale('symlog', linthreshy=10**-5)
    g6.set(ylim=(0, 10 ** -2))
    yaxis = plt.gca().yaxis
    yaxis.set_minor_locator(fits_new_plotter.MinorSymLogLocator(1e-1))
    g6.set_axis_labels("Passage", "Variant Frequency")
    # plt.show()
    g6.savefig(output_dir + "/Context_miseq_sample_time_line_plot", dpi=300)
    plt.close()


    # data_filter_uc = data_filter[data_filter["Mutation"] == "U>C"]
    #
    # data_filter_uc['Next'].replace('UA', 'UpA', inplace=True)
    # data_filter_uc['Next'].replace('UU', 'UpU', inplace=True)
    # data_filter_uc['Next'].replace('UC', 'UpC', inplace=True)
    # data_filter_uc['Next'].replace('UG', 'UpG', inplace=True)
    #
    #
    # data_filter_uc.to_csv(input_dir + "/data_filter_uc.csv", sep=',', encoding='utf-8')
    # context_order_uc = ["UpA", "UpU", "UpG",  "UpC"]
    #
    #
    # g7 = sns.catplot("label", "frac_and_weight", data=data_filter_uc, hue="Next", order=label_order, palette="tab20",
    #                     hue_order=context_order_uc, estimator=weighted_varaint, orient="v", dodge=True, kind="point",
    #                  col="Type", join=False, col_order=type_order)
    # g7.set_axis_labels("", "Variant Frequency")
    # g7.set(yscale='log')
    # g7.set(ylim=(10**-7, 10**-3))
    # g7.set_xticklabels(rotation=45)
    # # plt.show()
    # g7.savefig(output_dir + "/UC_Context_point_plot", dpi=300)
    # plt.close()
    #
    # data_filter_uc = data_filter_uc[data_filter_uc["passage"] != 0]
    # g8 = sns.relplot("passage", "frac_and_weight", data=data_filter_uc, hue="Next", palette="tab20",
    #                     hue_order=context_order_uc, estimator=weighted_varaint, col="Type", kind="line",
    #                  col_order=type_order)
    #
    # # g8.set(yscale="log")
    # g8.axes.flat[0].set_yscale('symlog', linthreshy=10**-4)
    # # g8.set(ylim=(0, 10 ** -6))
    # yaxis = plt.gca().yaxis
    # yaxis.set_minor_locator(MinorSymLogLocator(1e-1))
    # g8.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g8.savefig(output_dir + "/UC_Context_miseq_sample_time_line_plot", dpi=300)
    # plt.close()

if __name__ == "__main__":
    main()
