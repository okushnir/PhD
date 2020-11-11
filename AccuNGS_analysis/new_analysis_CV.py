
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
    date = "20201109"
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
    # filter based on pval<0.01 and Prob>0.95
    data_filter["no_variants"] = np.where(data_filter["pval"] > 0.01, 0, data_filter["no_variants"])
    data_filter["no_variants"] = np.where(data_filter["Prob"] < 0.95, 0, data_filter["no_variants"])
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
