
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from AccuNGS_analysis.adar_mutation_palette import mutation_palette
from scipy import stats

def weighted_varaint(x, **kws):
    var, count = map(np.asarray, zip(*x))
    return var.sum() / count.sum()


def linear_reg(data_filter, output_dir, transition_order, type_order, virus, replica, cu=True, ag=True, uc=True,
               ga= True, output_file="/mutation_rate"):
    """

    :param data_filter:
    :param output_dir:
    :param transition_order:
    :param type_order:
    :param virus:
    :param replica:
    :param cu:
    :param output_file:
    :return:
    """
    data_filter_grouped = data_filter.groupby(["label", "passage", "Type", "Mutation", "replica"])[
        "frac_and_weight"].agg(
        lambda x: weighted_varaint(x))
    data_filter_grouped = data_filter_grouped.reset_index()

    data_filter_grouped = data_filter_grouped.rename(columns={"frac_and_weight": "Frequency"})
    data_filter_grouped["Frequency"] = data_filter_grouped["Frequency"].astype(float)
    data_filter_grouped = data_filter_grouped[data_filter_grouped["label"] != "RNA Control\nPrimer ID"]
    data_filter_grouped = data_filter_grouped[data_filter_grouped["label"] != "RNA Control\nRND"]
    data_filter_grouped = data_filter_grouped[data_filter_grouped["replica"] == replica]
    print(data_filter_grouped.to_string())

    data_reg_ag = data_filter_grouped[data_filter_grouped["Mutation"] == "A>G"]
    data_reg_uc = data_filter_grouped[data_filter_grouped["Mutation"] == "U>C"]
    data_reg_ga = data_filter_grouped[data_filter_grouped["Mutation"] == "G>A"]
    data_reg_cu = data_filter_grouped[data_filter_grouped["Mutation"] == "C>U"]

    data_reg_ag_syn = data_reg_ag[data_reg_ag["Type"] == "Synonymous"]
    data_reg_uc_syn = data_reg_uc[data_reg_uc["Type"] == "Synonymous"]
    data_reg_ga_syn = data_reg_ga[data_reg_ga["Type"] == "Synonymous"]
    data_reg_cu_syn = data_reg_cu[data_reg_cu["Type"] == "Synonymous"]

    data_reg_ag_non_syn = data_reg_ag[data_reg_ag["Type"] == "Non-Synonymous"]
    data_reg_uc_non_syn = data_reg_uc[data_reg_uc["Type"] == "Non-Synonymous"]
    data_reg_ga_non_syn = data_reg_ga[data_reg_ga["Type"] == "Non-Synonymous"]
    data_reg_cu_non_syn = data_reg_cu[data_reg_cu["Type"] == "Non-Synonymous"]

    data_reg_ga_pmsc = data_reg_ga[data_reg_ga["Type"] == "Premature Stop Codon"]
    data_reg_cu_pmsc = data_reg_cu[data_reg_cu["Type"] == "Premature Stop Codon"]

    if ag==True:
        slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(data_reg_ag_syn['passage'],
                                                                                      data_reg_ag_syn
                                                                                      ['Frequency'])
        slope5, intercept5, r_value5, p_value5, std_err5 = stats.linregress(data_reg_ag_non_syn['passage'],
                                                                                      data_reg_ag_non_syn[
                                                                                          'Frequency'])
    if uc==True:
        slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(data_reg_uc_syn['passage'],
                                                                                      data_reg_uc_syn[
                                                                                          'Frequency'])
        slope6, intercept6, r_value6, p_value6, std_err6 = stats.linregress(data_reg_uc_non_syn['passage'],
                                                                            data_reg_uc_non_syn[
                                                                                'Frequency'])
    if ga == True:
        slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(data_reg_ga_syn['passage'],
                                                                                      data_reg_ga_syn[
                                                                                          'Frequency'])
        slope7, intercept7, r_value7, p_value7, std_err7 = stats.linregress(data_reg_ga_non_syn['passage'],
                                                                            data_reg_ga_non_syn[
                                                                                'Frequency'])
        slope11, intercept11, r_value11, p_value11, std_err11 = stats.linregress(data_reg_ga_pmsc['passage'],
                                                                                 data_reg_ga_pmsc[
                                                                                     'Frequency'])
    if cu == True:
        slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(data_reg_cu_syn['passage'],
                                                                                      data_reg_cu_syn[
                                                                                          'Frequency'])
        slope8, intercept8, r_value8, p_value8, std_err8 = stats.linregress(data_reg_cu_non_syn['passage'],
                                                                                      data_reg_cu_non_syn[
                                                                                          'Frequency'])
        slope12, intercept12, r_value12, p_value12, std_err12 = stats.linregress(data_reg_cu_pmsc['passage'],
                                                                                           data_reg_cu_pmsc[
                                                                                               'Frequency'])
    data_filter_grouped = data_filter_grouped.rename(columns={"passage": "Passage"})
    reg_plot = sns.lmplot(x="Passage", y="Frequency", data=data_filter_grouped, hue="Mutation",
                          hue_order=transition_order, fit_reg=True, col="Mutation",
                          col_order=transition_order, row="Type", row_order=type_order, palette=mutation_palette(4),
                          line_kws={'label': "Linear Reg"}, legend=True, height=6)  # markers=["o", "v", "x"]

    if ag == None:
        slope1, intercept1, r_value1, p_value1, std_err1 = 0, 0, 0, 0, 0
        slope5, intercept5, r_value5, p_value5, std_err5 = 0, 0, 0, 0, 0
    if uc == None:
        slope2, intercept2, r_value2, p_value2, std_err2 = 0, 0, 0, 0, 0
        slope6, intercept6, r_value6, p_value6, std_err6 = 0, 0, 0, 0, 0
    if ga == None:
        slope3, intercept3, r_value3, p_value3, std_err3 = 0, 0, 0, 0, 0
        slope7, intercept7, r_value7, p_value7, std_err7 = 0, 0, 0, 0, 0
        slope11, intercept11, r_value11, p_value11, std_err11 = 0, 0, 0, 0, 0
    if cu == None:
        slope4, intercept4, r_value4, p_value4, std_err4 = 0, 0, 0, 0, 0
        slope8, intercept8, r_value8, p_value8, std_err8 = 0, 0, 0, 0, 0
        slope12, intercept12, r_value12, p_value12, std_err12 = 0, 0, 0, 0, 0

    label_line_1 = "y={0:.3g}x+{1:.3g} pval={2:.3g} r-squared={3:.3g}".format(slope1, intercept1, p_value1,
                                                                              r_value1 ** 2)
    label_line_2 = "y={0:.3g}x+{1:.3g} pval={2:.3g} r-squared={3:.3g}".format(slope2, intercept2, p_value2,
                                                                              r_value2 ** 2)
    label_line_3 = "y={0:.3g}x+{1:.3g} pval={2:.3g} r-squared={3:.3g}".format(slope3, intercept3, p_value3,
                                                                              r_value3 ** 2)
    label_line_4 = "y={0:.3g}x+{1:.3g} pval={2:.3g} r-squared={3:.3g}".format(slope4, intercept4, p_value4,
                                                                              r_value4 ** 2)
    label_line_5 = "y={0:.3g}x+{1:.3g} pval={2:.3g} r-squared={3:.3g}".format(slope5, intercept5, p_value5,
                                                                              r_value5 ** 2)
    label_line_6 = "y={0:.3g}x+{1:.3g} pval={2:.3g} r-squared={3:.3g}".format(slope6, intercept6, p_value6,
                                                                              r_value6 ** 2)
    label_line_7 = "y={0:.3g}x+{1:.3g} pval={2:.3g} r-squared={3:.3g}".format(slope7, intercept7, p_value7,
                                                                              r_value7 ** 2)
    label_line_8 = "y={0:.3g}x+{1:.3g} pval={2:.3g} r-squared={3:.3g}".format(slope8, intercept8, p_value8,
                                                                              r_value8 ** 2)
    label_line_11 = "y={0:.3g}x+{1:.3g} pval={2:.3g} r-squared={3:.3g}".format(slope11, intercept11, p_value11,
                                                                              r_value11 ** 2)
    label_line_12 = "y={0:.3g}x+{1:.3g} pval={2:.3g} r-squared={3:.3g}".format(slope12, intercept12, p_value12,
                                                                              r_value12 ** 2)
    reg_plot.fig.subplots_adjust(wspace=.02)
    if ag == True:
        ax = reg_plot.axes[0, 0]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_1)
    if uc == True:
        ax = reg_plot.axes[0, 1]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_2)
    if ga == True:
        ax = reg_plot.axes[0, 2]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_3)
    if cu == True:
        ax = reg_plot.axes[0, 3]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_4)
    if ag == True:
        ax = reg_plot.axes[1, 0]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_5)
    if uc == True:
        ax = reg_plot.axes[1, 1]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_6)
    if ga == True:
        ax = reg_plot.axes[1, 2]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_7)
    if cu == True:
        ax = reg_plot.axes[1, 3]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_8)
    if ga == True:
        ax = reg_plot.axes[2, 2]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_11)
    if cu == True:
        ax = reg_plot.axes[2, 3]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_12)

    # reg_plot.set(xlim=(0, 13))
    # reg_plot.set(ylim=(0.000, 0.001))
    # reg_plot.fig.suptitle("RV #%s" % str(replica), y=0.99)
    # plt.tight_layout()
    reg_plot.savefig(output_dir + output_file + "_lmplot_%s.png" % replica, dpi=300)
    plt.close()

    columns = ["Mutation", "Type", "Slope", "Intercept", "pval"]
    mutation_rate_df = pd.DataFrame(columns=columns)
    mutation_rate_df.loc[0] = ["A>G", "Synonymous", slope1, intercept1, p_value1]
    mutation_rate_df.loc[1] = ["U>C", "Synonymous", slope2, intercept2, p_value2]
    mutation_rate_df.loc[2] = ["G>A", "Synonymous", slope3, intercept3, p_value3]
    mutation_rate_df.loc[3] = ["C>U", "Synonymous", slope4, intercept4, p_value3]
    mutation_rate_df.loc[4] = ["A>G", "Non-Synonymous", slope5, intercept5, p_value5]
    mutation_rate_df.loc[5] = ["U>C", "Non-Synonymous", slope6, intercept6, p_value6]
    mutation_rate_df.loc[6] = ["G>A", "Non-Synonymous", slope7, intercept7, p_value7]
    mutation_rate_df.loc[7] = ["C>U", "Non-Synonymous", slope8, intercept8, p_value8]
    mutation_rate_df.loc[8] = ["G>A", "Pre Mature Stop Codon", slope11, intercept11, p_value11]
    mutation_rate_df.loc[9] = ["C>U", "Pre Mature Stop Codon", slope12, intercept12, p_value11]
    mutation_rate_df["Virus"] = virus
    mutation_rate_df["Replica"] = replica
    mutation_rate_df.to_csv(output_dir + output_file + ".csv", sep=',', encoding='utf-8')
    mutation_rate_df.to_pickle(output_dir + output_file + str(replica) + ".pkl")

    return data_filter_grouped
