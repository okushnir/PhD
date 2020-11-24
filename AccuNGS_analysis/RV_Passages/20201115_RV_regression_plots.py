
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
from AccuNGS_analysis.adar_mutation_palette import mutation_palette
from AccuNGS_analysis.Linear_regression import linear_reg

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
    replica_lst = [1, 2, 3]
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/"
    prefix = "inosine_predict_context"
    output_dir = input_dir + "20201124_10000coverage_%s" % prefix
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)
    for replica in replica_lst:
        data_filter = pd.read_pickle(input_dir + prefix + "/data_filter.pkl")
        data_filter_ag = pd.read_pickle(input_dir + prefix + "/data_filter_ag.pkl")
        data_filter_uc = pd.read_pickle(input_dir + prefix +"/data_filter_uc.pkl")
        data_filter["passage"] = data_filter["passage"].astype(int)
        data_filter_ag["passage"] = data_filter_ag["passage"].astype(int)
        data_filter_uc["passage"] = data_filter_uc["passage"].astype(int)
        data_filter["replica"] = data_filter["replica"].astype(int)
        data_filter_ag["replica"] = data_filter_ag["replica"].astype(int)
        data_filter_uc["replica"] = data_filter_uc["replica"].astype(int)
        data_filter = data_filter[(data_filter['Read_count'] > 10000)]
        data_filter_ag = data_filter_ag[(data_filter_ag['Read_count'] > 10000)]
        data_filter_uc = data_filter_uc[(data_filter_uc['Read_count'] > 10000)]

        # data_filter = data_filter[data_filter["label"] != "p10-1"]
        # data_filter_ag = data_filter_ag[data_filter_ag["label"] != "p10-1"]
        # data_filter_uc = data_filter_uc[data_filter_uc["label"] != "p10-1"]

        data_filter["replica"] = np.where(data_filter["label"] == "RNA Control\nPrimer ID", replica,
                                          data_filter["replica"])
        data_filter_ag["replica"] = np.where(data_filter_ag["label"] == "RNA Control\nPrimer ID", replica,
                                          data_filter_ag["replica"])
        data_filter_uc["replica"] = np.where(data_filter_uc["label"] == "RNA Control\nPrimer ID", replica,
                                          data_filter_uc["replica"])

        data_filter = data_filter[data_filter["replica"] == replica]
        data_filter_ag = data_filter_ag[data_filter_ag["replica"] == replica]
        data_filter_uc = data_filter_uc[data_filter_uc["replica"] == replica]

        """Plots"""
        transition_order = ["A>G", "U>C", "G>A", "C>U"]
        type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]
        type_order_ag = ["Synonymous", "Non-Synonymous"]
        adar_preference = ["High", "Intermediate", "Low"]

        linear_reg(data_filter, output_dir, transition_order, type_order, replica=replica, virus="RVB14")

        """A>G Prev Context"""
        mutation_ag = sns.catplot("passage", "frac_and_weight", data=data_filter_ag, hue="5`_ADAR_Preference",
                         palette=mutation_palette(3, adar=True, ag=True), kind="point", dodge=True, estimator=weighted_varaint,
                         orient="v", col="Type", join=False, col_order=type_order_ag, hue_order=adar_preference)
        mutation_ag.set(yscale="log")
        mutation_ag.set(ylim=(1*10**-5, 1*10**-2))
        mutation_ag.fig.suptitle("A>G ADAR_like Mutation in RV #%s" % replica, y=0.99)
        plt.subplots_adjust(top=0.85)
        mutation_ag.set_axis_labels("Passage", "Variant Frequency")
        mutation_ag.savefig(output_dir + "/ag_ADAR_like_Mutation_col_replica%s.png" % replica, dpi=300)
        plt.close()

        # data_filter_ag_grouped = data_filter_ag.groupby(["label", "Type", "passage", "replica",
        #                                                  "ADAR_grade_five", "5`_ADAR_Preference"])[
        #     "frac_and_weight"].agg(
        #     lambda x: weighted_varaint(x))
        # data_filter_ag_grouped = data_filter_ag_grouped.reset_index()

        # print(data_filter_ag_grouped.to_string())

        # data_filter_ag = data_filter_ag_grouped.rename(columns={"frac_and_weight": "Frequency"})
        # data_filter_ag_grouped["Frequency"] = data_filter_ag_grouped["Frequency"].astype(float)
        data_filter_ag = data_filter_ag[data_filter_ag["label"] != "RNA Control\nPrimer ID"]
        data_filter_ag = data_filter_ag[data_filter_ag["label"] != "RNA Control\nRND"]
        data_filter_ag = data_filter_ag[data_filter_ag["replica"] == replica]

        data_reg_full_adar = data_filter_ag[data_filter_ag["ADAR_grade_five"] == 1]
        data_reg_semi_adar = data_filter_ag[data_filter_ag["ADAR_grade_five"] == 0.5]
        data_reg_nonadar = data_filter_ag[data_filter_ag["ADAR_grade_five"] == 0]

        data_reg_full_adar_syn = data_reg_full_adar[data_reg_full_adar["Type"] == "Synonymous"]
        data_reg_semi_adar_syn = data_reg_semi_adar[data_reg_semi_adar["Type"] == "Synonymous"]
        data_reg_nonadar_syn = data_reg_nonadar[data_reg_nonadar["Type"] == "Synonymous"]

        data_reg_full_adar_non_syn = data_reg_full_adar[data_reg_full_adar["Type"] == "Non-Synonymous"]
        data_reg_semi_adar_non_syn = data_reg_semi_adar[data_reg_semi_adar["Type"] == "Non-Synonymous"]
        data_reg_nonadar_non_syn = data_reg_nonadar[data_reg_nonadar["Type"] == "Non-Synonymous"]

        slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(data_reg_full_adar_syn['passage'],
                                                                                      data_reg_full_adar_syn[
                                                                                          'Frequency'])
        slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(data_reg_semi_adar_syn['passage'],
                                                                                      data_reg_semi_adar_syn[
                                                                                          'Frequency'])
        slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(data_reg_nonadar_syn['passage'],
                                                                                      data_reg_nonadar_syn[
                                                                                          'Frequency'])
        # data_reg_adar_nonsyn
        slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(
            data_reg_full_adar_non_syn['passage'],
            data_reg_full_adar_non_syn[
                'Frequency'])
        slope5, intercept5, r_value5, p_value5, std_err5 = stats.linregress(
            data_reg_semi_adar_non_syn['passage'],
            data_reg_semi_adar_non_syn[
                'Frequency'])
        slope6, intercept6, r_value6, p_value6, std_err6 = stats.linregress(
            data_reg_nonadar_non_syn['passage'],
            data_reg_nonadar_non_syn[
                'Frequency'])
        data_filter_ag = data_filter_ag.rename(columns={"passage": "Passage"})
        ag_reg_plot = sns.lmplot(x="Passage", y="Frequency", data=data_filter_ag, hue="5`_ADAR_Preference",
                                 hue_order=adar_preference, markers=["o", "v", "x"], fit_reg=True,
                                 col="5`_ADAR_Preference",
                                 col_order=adar_preference, row="Type", row_order=type_order_ag,
                                 palette=mutation_palette(3, adar=True, ag=True),
                                 line_kws={'label': "Linear Reg"}, legend=True, height=6)
        ag_reg_plot.fig.subplots_adjust(wspace=.02)
        ax = ag_reg_plot.axes[0, 0]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        label_line_1 = "y={0:.3g}x+{1:.3g} pval={2:.3g}".format(slope1, intercept1, p_value1)
        label_line_2 = "y={0:.3g}x+{1:.3g} pval={2:.3g}".format(slope2, intercept2, p_value2)
        label_line_3 = "y={0:.3g}x+{1:.3g} pval={2:.3g}".format(slope3, intercept3, p_value3)
        label_line_4 = "y={0:.3g}x+{1:.3g} pval={2:.3g}".format(slope4, intercept4, p_value4)
        label_line_5 = "y={0:.3g}x+{1:.3g} pval={2:.3g}".format(slope5, intercept5, p_value5)
        label_line_6 = "y={0:.3g}x+{1:.3g} pval={2:.3g}".format(slope6, intercept6, p_value6)
        L_labels[0].set_text(label_line_1)
        ax = ag_reg_plot.axes[0, 1]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_2)
        ax = ag_reg_plot.axes[0, 2]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_3)
        ax = ag_reg_plot.axes[1, 0]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_4)
        ax = ag_reg_plot.axes[1, 1]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_5)
        ax = ag_reg_plot.axes[1, 2]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        L_labels[0].set_text(label_line_6)

        ag_reg_plot.set(xlim=(0, 13))
        ag_reg_plot.set(ylim=(0.000, 0.003))
        # ag_reg_plot.fig.suptitle("RV #%s" % str(replica), y=0.99)
        # plt.tight_layout()
        ag_reg_plot.savefig(output_dir + "/ag_lmplot_ADAR_Context_replica%s" % replica, dpi=300)
        plt.close()

        data_filter_ag = data_filter_ag[data_filter_ag["Protein"] != "2A"]
        data_filter_ag = data_filter_ag[data_filter_ag["Protein"] != "3'UTR"]
        data_filter_ag = data_filter_ag[data_filter_ag["Type"] == "Synonymous"]

        position_mutation_ag = sns.relplot(x="Pos", y="Frequency", data=data_filter_ag, hue="5`_ADAR_Preference",
                                           col="Passage", col_wrap=3, palette=mutation_palette(3, adar=True, ag=True),
                                           hue_order=adar_preference, height=4,
                                           style="5`_ADAR_Preference", style_order=["High", "Low", "Intermediate"])

        position_mutation_ag.set_axis_labels("", "Variant Frequency")
        position_mutation_ag.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -4)
        position_mutation_ag.axes.flat[0].set_ylim(10 ** -5, 10 ** -2)
        plt.savefig(output_dir + "/ag_position_mutation_replica%s.png" % replica, dpi=300)
        plt.close()

        columns = ["Mutation", "Replica", "Type", "Slope", "Intercept", "p-val"]
        mutation_rate_ag_df = pd.DataFrame(columns=columns)
        mutation_rate_ag_df.loc[0] = ["High ADAR-like A>G", replica, "Synonymous", slope1, intercept1, p_value1]
        mutation_rate_ag_df.loc[1] = ["Intermediate ADAR-like A>G", replica, "Synonymous", slope2, intercept2, p_value2]
        mutation_rate_ag_df.loc[2] = ["Low ADAR-like A>G", replica, "Synonymous", slope3, intercept3, p_value3]
        mutation_rate_ag_df.loc[3] = ["High ADAR-like A>G", replica, "Non-Synonymous", slope4, intercept4, p_value4]
        mutation_rate_ag_df.loc[4] = ["Intermediate ADAR-like A>G", replica, "Non-Synonymous", slope5,
                                      intercept5, p_value5]
        mutation_rate_ag_df.loc[5] = ["Low ADAR-like A>G", replica, "Non-Synonymous", slope6, intercept6, p_value6]

        mutation_rate_ag_df.to_csv(output_dir + "/mutation_rate_ag%s.csv" % replica, sep=',', encoding='utf-8')
        mutation_rate_ag_df.to_pickle(output_dir + "/mutation_rate_ag%s.pkl" % replica)

        """U>C Context"""
        mutation_uc = sns.catplot("passage", "frac_and_weight", data=data_filter_uc, hue="3`_ADAR_Preference",
                                  palette=mutation_palette(3, adar=True, uc=True), kind="point", dodge=True, estimator=weighted_varaint,
                                  orient="v", col="Type", join=False, hue_order=adar_preference,
                                  col_order=type_order_ag)
        mutation_uc.set(yscale="log")
        mutation_uc.set(ylim=(1 * 10 ** -5, 1 * 10 ** -2))
        mutation_uc.set_axis_labels("Passage", "Variant Frequency")
        mutation_uc.savefig(output_dir + "/uc_ADAR_like_Mutation_col_replica%s.png" % replica, dpi=300)
        plt.close()

        # data_filter_uc_grouped = data_filter_uc.groupby(["label", "Type", "passage", "replica",
        #                                                  "ADAR_grade_three", "3`_ADAR_Preference"])["frac_and_weight"].agg(
        #     lambda x: weighted_varaint(x))
        #
        # data_filter_uc_grouped = data_filter_uc_grouped.reset_index()
        # data_filter_uc_grouped = data_filter_uc_grouped.rename(columns={"frac_and_weight": "Frequency"})
        # data_filter_uc_grouped["Frequency"] = data_filter_uc_grouped["Frequency"].astype(float)
        data_filter_uc = data_filter_uc[data_filter_uc["label"] != "RNA Control\nPrimer ID"]
        data_filter_uc = data_filter_uc[data_filter_uc["label"] != "RNA Control\nRND"]
        data_filter_uc = data_filter_uc[data_filter_uc["replica"] == replica]

        data_reg_full_adar_uc = data_filter_uc[data_filter_uc["ADAR_grade_three"] == 1]
        data_reg_semi_adar_uc = data_filter_uc[data_filter_uc["ADAR_grade_three"] == 0.5]
        data_reg_nonadar_uc = data_filter_uc[data_filter_uc["ADAR_grade_three"] == 0]

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
        data_filter_uc = data_filter_uc.rename(columns={"passage": "Passage"})
        uc_lmplot = sns.lmplot(x="Passage", y="Frequency", data=data_filter_uc, hue="3`_ADAR_Preference",
                               markers=["o", "v", "x"], hue_order=adar_preference, fit_reg=True, col="3`_ADAR_Preference",
                               col_order=adar_preference, row="Type", row_order=type_order_ag,
                               palette=mutation_palette(3, adar=True, uc=True),
                               line_kws={'label': "Linear Reg"}, legend=True, height=6)
        uc_lmplot.fig.subplots_adjust(wspace=.02)
        ax = uc_lmplot.axes[0, 0]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        label_line_7 = "y={0:.3g}x+{1:.3g} pval={2:.3g}".format(stat_slope7, stat_intercept7, p_value7)
        label_line_8 = "y={0:.3g}x+{1:.3g} pval={2:.3g}".format(stat_slope8, stat_intercept8, p_value8)
        label_line_9 = "y={0:.3g}x+{1:.3g} pval={2:.3g}".format(stat_slope9, stat_intercept9, p_value9)
        label_line_10 = "y={0:.3g}x+{1:.3g} pval={2:.3g}".format(stat_slope10, stat_intercept10, p_value10)
        label_line_11 = "y={0:.3g}x+{1:.3g} pval={2:.3g}".format(stat_slope11, stat_intercept11, p_value11)
        label_line_12 = "y={0:.3g}x+{1:.3g} pval={2:.3g}".format(stat_slope12, stat_intercept12, p_value12)
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
        uc_lmplot.set(ylim=(0.000, 0.003))
        # uc_lmplot.fig.suptitle("RV #%s" % str(replica), y=0.99)
        # plt.tight_layout()
        uc_lmplot.savefig(output_dir + "/uc_lmplot_ADAR_Context_replica%s" % replica, dpi=300)
        plt.close()

        data_filter_uc = data_filter_uc[data_filter_uc["Protein"] != "2A"]
        data_filter_uc = data_filter_uc[data_filter_uc["Protein"] != "3'UTR"]
        data_filter_uc = data_filter_uc[data_filter_uc["Type"] == "Synonymous"]

        position_mutation_uc = sns.relplot(x="Pos", y="Frequency", data=data_filter_uc, hue="3`_ADAR_Preference",
                                        col="Passage", col_wrap=3, palette=mutation_palette(3, adar=True, uc=True),
                                        hue_order=adar_preference, height=4,
                                           style="3`_ADAR_Preference", style_order=["High", "Low", "Intermediate"])

    mutation_rate_df1 = pd.read_pickle(output_dir + "/mutation_rate1.pkl")
    mutation_rate_df2 = pd.read_pickle(output_dir + "/mutation_rate2.pkl")
    mutation_rate_df3 = pd.read_pickle(output_dir + "/mutation_rate3.pkl")
    mutation_rate_df_all = pd.concat([mutation_rate_df1, mutation_rate_df2, mutation_rate_df3], sort=False)
    mutation_rate_df_all.to_csv(output_dir + "/mutation_rate_all.csv", sep=',', encoding='utf-8')
    mutation_rate_df_all_grouped = mutation_rate_df_all.groupby(["Mutation", "Type"])["Slope", "Intercept"].agg(np.median)
    mutation_rate_df_all_grouped = mutation_rate_df_all_grouped.reset_index()
    mutation_rate_df_all_grouped.to_csv(output_dir + "/mutation_rate_median.csv", sep=',', encoding='utf-8')

    mutation_rate_ag_df1 = pd.read_pickle(output_dir + "/mutation_rate_ag1.pkl")
    mutation_rate_ag_df2 = pd.read_pickle(output_dir + "/mutation_rate_ag2.pkl")
    mutation_rate_ag_df3 = pd.read_pickle(output_dir + "/mutation_rate_ag3.pkl")
    mutation_rate_ag_df_all = pd.concat([mutation_rate_ag_df1, mutation_rate_ag_df2, mutation_rate_ag_df3], sort=False)
    mutation_rate_ag_df_all.to_csv(output_dir + "/mutation_rate_ag_all.csv", sep=',', encoding='utf-8')
    mutation_rate_ag_df_all.to_pickle(output_dir + "/mutation_rate_ag_all.pkl")

    # mutation_rate_ag_df = pd.read_csv(output_dir + "/mutation_rate_ag_all.csv", sep=',', encoding='utf-8')
    mutation_rate_ag_df_all_grouped = mutation_rate_ag_df_all.groupby(["Mutation", "Type"])["Slope", "Intercept"].agg(np.median)
    mutation_rate_ag_df_all_grouped = mutation_rate_ag_df_all_grouped.reset_index()
    mutation_rate_ag_df_all_grouped.to_csv(output_dir + "/mutation_rate_ag_median.csv", sep=',', encoding='utf-8')


if __name__ == "__main__":
    main()
