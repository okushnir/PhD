
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

def checkKey(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        raise Exception()


def main():
    input_dir = "/Users/odedkushnir/Projects/fitness/CirSeq/PV/Mahoney/"
    output_dir = input_dir + "20201112_plots"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)


    data_mutations = pd.read_csv(input_dir + "q30_data_mutation.csv")

    columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
               "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts", "Consensus>Mutated_codon"]
    data_filter = pd.DataFrame(data_mutations, columns=columns)
    data_filter["pval"] = data_filter["pval"].fillna(1)
    data_filter["no_variants"] = data_filter["Frequency"] * data_filter["Read_count"]

    data_filter["frac_and_weight"] = list(zip(data_filter.no_variants, data_filter.Read_count))
    data_filter["passage"] = data_filter["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    data_filter["passage"] = np.where(data_filter["passage"] == "RNA Control", 0, data_filter["passage"])
    data_filter["passage"] = data_filter["passage"].astype(int)
    data_filter["Type"] = data_filter["Type"].fillna("NonCodingRegion")
    data_filter = data_filter[data_filter["Mutation"] != "C>U"]
    data_filter.to_csv(output_dir + "/data_filter.csv", sep=',', encoding='utf-8')

    label_order = ["PV-p3", "PV-p4", "PV-p5", "PV-p6", "PV-p7", "PV-p8"]
    passage_order = ["3", "4", "5", "6", "7", "8"]
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U"]
    transition_order = ["A>G", "U>C", "G>A"]
    type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]

    g1 = sns.catplot("label", "frac_and_weight", data=data_filter, hue="Mutation", order=label_order, palette="tab20",
                        kind="point", dodge=True, hue_order=mutation_order, join=False, estimator=weighted_varaint,
                     orient="v")
    g1.set_axis_labels("", "Variant Frequency")
    g1.set_xticklabels(fontsize=9, rotation=45)
    g1.set(yscale='log')
    g1.set(ylim=(10**-5, 10**-1))

    # plt.show()
    g1.savefig(output_dir + "/All_Mutations_point_plot", dpi=300)
    plt.close()

    data_filter["passage"] = data_filter["passage"].astype(str)
    g2 = sns.catplot("passage", "frac_and_weight", data=data_filter, hue="Mutation", order=passage_order, palette=mutation_palette(4)
                        ,kind="point", dodge=True, hue_order=transition_order, join=False, estimator=weighted_varaint,
                     orient="v")
    g2.set_axis_labels("Passage", "Variant Frequency")
    g2.set(yscale='log')
    g2.set(ylim=(10 ** -6, 10 ** -2))
    # g2.set_xticklabels(fontsize=10, rotation=45)
    # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/Prgress reports/20200913 Final report/plots" +
    #                   "/Transition_Mutations_point_plot_Mahoney", dpi=300)
    g2.savefig(output_dir + "/Transition_Mutations_point_plot_Mahoney", dpi=300)
    plt.close()

    data_filter["passage"] = data_filter["passage"].astype(int)

    linear_reg(data_filter, output_dir, transition_order, type_order, virus="PV1", replica=1, cu=None)

    # A>G Prev Context
    data_filter_ag = data_filter[data_filter["Mutation"] == "A>G"]

    data_filter_ag['Prev'].replace('AA', 'ApA', inplace=True)
    data_filter_ag['Prev'].replace('UA', 'UpA', inplace=True)
    data_filter_ag['Prev'].replace('CA', 'CpA', inplace=True)
    data_filter_ag['Prev'].replace('GA', 'GpA', inplace=True)

    context_order = ["UpA", "ApA", "CpA", "GpA"]
    type_order = ["Synonymous", "Non-Synonymous"]
    data_filter_ag["ADAR_like"] = data_filter_ag.Prev.str.contains('UpA') | data_filter_ag.Prev.str.contains('ApA')
    data_filter_ag_pass8 = data_filter_ag.loc[data_filter_ag.passage == 8]
    data_filter_ag_pass8 = data_filter_ag_pass8.loc[data_filter_ag_pass8.Type == "Synonymous"]
    print(data_filter_ag_pass8.to_string())

    g5 = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=label_order, palette=mutation_palette(2),
                kind = "point", dodge=True, hue_order=[True, False], estimator=weighted_varaint, orient="v",
                     col="Type", join=False, col_order=type_order)
    g5.set_axis_labels("", "Variant Frequency")
    g5.set(yscale='log')
    g5.set(ylim=(7*10**-7, 4*10**-3))
    g5.set_xticklabels(rotation=45)
    # plt.show()
    g5.savefig(output_dir + "/Context_point_plot", dpi=300)
    plt.close()
    #

    ax = sns.boxplot("ADAR_like", "Frequency", data=data_filter_ag_pass8, palette=mutation_palette(2), order=[True, False])
    ax = sns.stripplot("ADAR_like", "Frequency", data=data_filter_ag_pass8, color=".2", order=[True, False])
    old_statannot.add_stat_annotation(ax, data=data_filter_ag_pass8, x="ADAR_like", y="Frequency",
                        boxPairList=[(True, False)], test='Mann-Whitney', textFormat='star', loc='inside', verbose=1)
    ax.set_yscale('log')
    ax.set_xlabel("ADAR-like\nContext")
    ax.set_ylabel("Variant Frequency")
    ax.set(ylim=(10 ** -4, 10 ** -2))
    # sns.despine()
    # plt.tight_layout()
    plt.savefig(output_dir + "/context_p8_point_plot_v2.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
