
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from AccuNGS_analysis.adar_mutation_palette import mutation_palette
from scipy import stats

def weighted_varaint(x, **kws):
    var, count = map(np.asarray, zip(*x))
    return var.sum() / count.sum()

def linear_reg(data_filter, output_dir, transition_order, type_order, virus, replica, cu=True):
    data_filter_grouped = data_filter.groupby(["label", "passage", "Type", "Mutation"])[
        "frac_and_weight"].agg(
        lambda x: weighted_varaint(x))
    data_filter_grouped = data_filter_grouped.reset_index()

    data_filter_grouped = data_filter_grouped.rename(columns={"frac_and_weight": "Frequency"})
    data_filter_grouped["Frequency"] = data_filter_grouped["Frequency"].astype(float)
    data_filter_grouped = data_filter_grouped[data_filter_grouped["label"] != "RNA Control\nPrimer ID"]
    data_filter_grouped = data_filter_grouped[data_filter_grouped["label"] != "RNA Control\nRND"]
    # data_filter_grouped = data_filter_grouped[data_filter_grouped["replica"] == replica]
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

    stat_slope1, stat_intercept1, r_value1, p_value1, std_err1 = stats.linregress(data_reg_ag_syn['passage'],
                                                                                  data_reg_ag_syn
                                                                                  ['Frequency'])
    stat_slope2, stat_intercept2, r_value2, p_value2, std_err2 = stats.linregress(data_reg_uc_syn['passage'],
                                                                                  data_reg_uc_syn[
                                                                                      'Frequency'])
    stat_slope3, stat_intercept3, r_value3, p_value3, std_err3 = stats.linregress(data_reg_ga_syn['passage'],
                                                                                  data_reg_ga_syn[
                                                                                      'Frequency'])
    if cu == True:
        stat_slope4, stat_intercept4, r_value4, p_value4, std_err4 = stats.linregress(data_reg_cu_syn['passage'],
                                                                                      data_reg_cu_syn[
                                                                                          'Frequency'])
        stat_slope8, stat_intercept8, r_value8, p_value8, std_err8 = stats.linregress(data_reg_cu_non_syn['passage'],
                                                                                      data_reg_cu_non_syn[
                                                                                          'Frequency'])
        stat_slope12, stat_intercept12, r_value12, p_value12, std_err12 = stats.linregress(data_reg_cu_pmsc['passage'],
                                                                                           data_reg_cu_pmsc[
                                                                                               'Frequency'])
    # data_reg_adar_nonsyn
    stat_slope5, stat_intercept5, r_value5, p_value5, std_err5 = stats.linregress(data_reg_ag_non_syn['passage'],
                                                                                  data_reg_ag_non_syn[
                                                                                      'Frequency'])
    stat_slope6, stat_intercept6, r_value6, p_value6, std_err6 = stats.linregress(data_reg_uc_non_syn['passage'],
                                                                                  data_reg_uc_non_syn[
                                                                                      'Frequency'])
    stat_slope7, stat_intercept7, r_value7, p_value7, std_err7 = stats.linregress(data_reg_ga_non_syn['passage'],
                                                                                  data_reg_ga_non_syn[
                                                                                      'Frequency'])

    # pmsc
    stat_slope11, stat_intercept11, r_value11, p_value11, std_err11 = stats.linregress(data_reg_ga_pmsc['passage'],
                                                                                       data_reg_ga_pmsc[
                                                                                           'Frequency'])

    data_filter_grouped = data_filter_grouped.rename(columns={"passage": "Passage"})
    reg_plot = sns.lmplot(x="Passage", y="Frequency", data=data_filter_grouped, hue="Mutation",
                          hue_order=transition_order, fit_reg=True, col="Mutation",
                          col_order=transition_order, row="Type", row_order=type_order, palette=mutation_palette(4),
                          line_kws={'label': "Linear Reg"}, legend=True, height=6)  # markers=["o", "v", "x"]
    reg_plot.fig.subplots_adjust(wspace=.02)
    ax = reg_plot.axes[0, 0]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    if cu == None:
        stat_slope4 = 0
        stat_intercept4 = 0
        stat_slope8 = 0
        stat_intercept8 = 0
        stat_slope12 = 0
        stat_intercept12 = 0

    label_line_1 = "y={0:.3g}x+{1:.3g}".format(stat_slope1, stat_intercept1)
    label_line_2 = "y={0:.3g}x+{1:.3g}".format(stat_slope2, stat_intercept2)
    label_line_3 = "y={0:.3g}x+{1:.3g}".format(stat_slope3, stat_intercept3)
    label_line_4 = "y={0:.3g}x+{1:.3g}".format(stat_slope4, stat_intercept4)
    label_line_5 = "y={0:.3g}x+{1:.3g}".format(stat_slope5, stat_intercept5)
    label_line_6 = "y={0:.3g}x+{1:.3g}".format(stat_slope6, stat_intercept6)
    label_line_7 = "y={0:.3g}x+{1:.3g}".format(stat_slope7, stat_intercept7)
    label_line_8 = "y={0:.3g}x+{1:.3g}".format(stat_slope8, stat_intercept8)
    label_line_11 = "y={0:.3g}x+{1:.3g}".format(stat_slope11, stat_intercept11)
    label_line_12 = "y={0:.3g}x+{1:.3g}".format(stat_slope12, stat_intercept12)



    L_labels[0].set_text(label_line_1)
    ax = reg_plot.axes[0, 1]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    L_labels[0].set_text(label_line_2)
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
    ax = reg_plot.axes[1, 0]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    L_labels[0].set_text(label_line_5)
    ax = reg_plot.axes[1, 1]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    L_labels[0].set_text(label_line_6)
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
    reg_plot.savefig(output_dir + "/transition_lmplot.png", dpi=300)
    plt.close()

    columns = ["Mutation", "Type", "Slope", "Intercept"]
    mutation_rate_df = pd.DataFrame(columns=columns)
    mutation_rate_df.loc[0] = ["A>G", "Synonymous", stat_slope1, stat_intercept1]
    mutation_rate_df.loc[1] = ["U>C", "Synonymous", stat_slope2, stat_intercept2]
    mutation_rate_df.loc[2] = ["G>A", "Synonymous", stat_slope3, stat_intercept3]
    mutation_rate_df.loc[3] = ["C>U", "Synonymous", stat_slope4, stat_intercept4]
    mutation_rate_df.loc[4] = ["A>G", "Non-Synonymous", stat_slope5, stat_intercept5]
    mutation_rate_df.loc[5] = ["U>C", "Non-Synonymous", stat_slope6, stat_intercept6]
    mutation_rate_df.loc[6] = ["G>A", "Non-Synonymous", stat_slope7, stat_intercept7]
    mutation_rate_df.loc[7] = ["C>U", "Non-Synonymous", stat_slope8, stat_intercept8]
    mutation_rate_df.loc[8] = ["G>A", "Pre Mature Stop Codon", stat_slope11, stat_intercept11]
    mutation_rate_df.loc[9] = ["C>U", "Pre Mature Stop Codon", stat_slope12, stat_intercept12]
    mutation_rate_df["Virus"] = virus
    mutation_rate_df["Replica"] = replica
    mutation_rate_df.to_csv(output_dir + "/mutation_rate.csv", sep=',', encoding='utf-8')
    mutation_rate_df.to_pickle(output_dir + "/mutation_rate.pkl")

    return data_filter_grouped
