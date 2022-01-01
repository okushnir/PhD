import datetime

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
from datetime import datetime

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
    input_dir = "/Users/odedkushnir/PhD_Projects/After_review/CirSeq/PV/Mahoney/"
    date = datetime.today().strftime("%Y%m%d")
    print(date)
    prefix = "inosine_predict_context"
    output_dir = input_dir + prefix
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)


    data_mutations = pd.read_csv(input_dir + "q30_data_mutation.csv")
    data_adar = pd.read_csv("/Users/odedkushnir/PhD_Projects/fitness/CirSeq/PV/Mahoney/inosine_results/Output/PV1_adar1_trans.csv")
    data_mutations = data_mutations.merge(data_adar, on="Pos", how="inner")

    columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
               "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts",
               "Consensus>Mutated_codon", "fiveGrade", "threeGrade"]
    data_filter = pd.DataFrame(data_mutations, columns=columns)
    data_filter["pval"] = data_filter["pval"].fillna(1)
    data_filter["no_variants"] = data_filter["Frequency"] * data_filter["Read_count"]
    # filter based on pval<0.01 and Prob>0.95
    # data_filter["no_variants"] = np.where(data_filter["pval"] > 0.01, 0, data_filter["no_variants"])
    data_filter["no_variants"] = np.where(data_filter["Prob"] < 0.95, 0, data_filter["no_variants"])

    data_filter["frac_and_weight"] = list(zip(data_filter.no_variants, data_filter.Read_count))
    data_filter["passage"] = data_filter["label"].apply(lambda x: x.split("-")[-1][1])
    data_filter["passage"] = data_filter["passage"].astype(int)
    data_filter["Type"] = data_filter["Type"].fillna("NonCodingRegion")
    data_filter = data_filter.loc[data_filter.Mutation != "C>U"]

    data_filter_ag = data_filter[data_filter["Mutation"] == "A>G"]
    data_filter_uc = data_filter[data_filter["Mutation"] == "U>C"]

    print("25_quantile_ag %s" % str(data_filter_ag["fiveGrade"].quantile(0.25)))
    print("75_quantile_ag %s" % str(data_filter_ag["fiveGrade"].quantile(0.75)))
    print("25_quantile_uc %s" % str(data_filter_uc["threeGrade"].quantile(0.25)))
    print("75_quantile_uc %s" % str(data_filter_uc["threeGrade"].quantile(0.75)))

    data_filter["ADAR_grade_five"] = np.where(data_filter["fiveGrade"] < data_filter_ag["fiveGrade"].quantile(0.25), 0,
                                              np.where(data_filter["fiveGrade"] <= data_filter_ag["fiveGrade"].
                                                       quantile(0.75), 0.5, 1))
    data_filter["5`_ADAR_Preference"] = np.where(data_filter["fiveGrade"] < data_filter_ag["fiveGrade"].quantile(0.25),
                                                 "Low", np.where(data_filter["fiveGrade"] <=
                                                                 data_filter_ag["fiveGrade"].quantile(0.75),
                                                                 "Intermediate", "High"))
    data_filter["ADAR_grade_five"] = np.where(data_filter["fiveGrade"] == 0, "0", data_filter["ADAR_grade_five"])
    data_filter["5`_ADAR_Preference"] = np.where(data_filter["fiveGrade"] == 0, "Low", data_filter["5`_ADAR_Preference"])

    data_filter["ADAR_grade_three"] = np.where(data_filter["threeGrade"] < data_filter_uc["threeGrade"].quantile(0.25), 0,
                                               np.where(data_filter["threeGrade"] <= data_filter_uc["threeGrade"].
                                                        quantile(0.75), 0.5, 1))
    data_filter["3`_ADAR_Preference"] = np.where(data_filter["threeGrade"] < data_filter_uc["threeGrade"].quantile(0.25),
                                                 "Low", np.where(data_filter["threeGrade"] <=
                                                                 data_filter_uc["threeGrade"].quantile(0.75),
                                                                 "Intermediate", "High"))
    data_filter["ADAR_grade_three"] = np.where(data_filter["threeGrade"] == 0, "0", data_filter["ADAR_grade_five"])
    data_filter["3`_ADAR_Preference"] = np.where(data_filter["threeGrade"] == 0, "Low", data_filter["3`_ADAR_Preference"])

    data_filter_ag = data_filter[data_filter["Mutation"] == "A>G"]
    data_filter_ag['Prev'].replace('AA', 'ApA', inplace=True)
    data_filter_ag['Prev'].replace('UA', 'UpA', inplace=True)
    data_filter_ag['Prev'].replace('CA', 'CpA', inplace=True)
    data_filter_ag['Prev'].replace('GA', 'GpA', inplace=True)
    data_filter_ag["ADAR_like"] = data_filter_ag.Prev.str.contains('UpA') | data_filter_ag.Prev.str.contains('ApA')

    data_filter_uc = data_filter[data_filter["Mutation"] == "U>C"]
    data_filter_uc['Next'].replace('AA', 'ApA', inplace=True)
    data_filter_uc['Next'].replace('UA', 'UpA', inplace=True)
    data_filter_uc['Next'].replace('CA', 'CpA', inplace=True)
    data_filter_uc['Next'].replace('GA', 'GpA', inplace=True)
    data_filter_uc["ADAR_like"] = data_filter_uc.Next.str.contains('UpA') | data_filter_uc.Next.str.contains('ApA')

    data_filter.to_csv(output_dir + "/data_filter.csv", sep=',', encoding='utf-8')
    data_filter_ag.to_csv(output_dir + "/data_filter_ag.csv", sep=',', encoding='utf-8')
    data_filter_uc.to_csv(output_dir + "/data_filter_uc.csv", sep=',', encoding='utf-8')
    data_filter.to_pickle(output_dir + "/data_filter.pkl")
    data_filter_ag.to_pickle(output_dir + "/data_filter_ag.pkl")
    data_filter_uc.to_pickle(output_dir + "/data_filter_uc.pkl")

    # columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
    #            "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts", "Consensus>Mutated_codon"]
    # data_filter = pd.DataFrame(data_mutations, columns=columns)
    # data_filter["pval"] = data_filter["pval"].fillna(1)
    # data_filter["no_variants"] = data_filter["Frequency"] * data_filter["Read_count"]
    #
    # data_filter["frac_and_weight"] = list(zip(data_filter.no_variants, data_filter.Read_count))
    # data_filter["passage"] = data_filter["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    # data_filter["passage"] = np.where(data_filter["passage"] == "RNA Control", 0, data_filter["passage"])
    # data_filter["passage"] = data_filter["passage"].astype(int)
    # data_filter["Type"] = data_filter["Type"].fillna("NonCodingRegion")
    # data_filter = data_filter[data_filter["Mutation"] != "C>U"]
    # data_filter.to_csv(output_dir + "/data_filter.csv", sep=',', encoding='utf-8')
    """Plots"""
    output_dir = input_dir + date + "_plots"
    plus_minus = u"\u00B1"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)
    label_order = ["PV-p3", "PV-p4", "PV-p5", "PV-p6", "PV-p7", "PV-p8"]
    passage_order = ["p3", "p4", "p5", "p6", "p7", "p8"]
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U"]
    transition_order = ["A>G", "U>C", "G>A"]
    type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]

    g1 = sns.catplot("label", "frac_and_weight", data=data_filter, hue="Mutation", order=label_order, palette="tab20",
                        kind="point", dodge=True, hue_order=mutation_order, join=False, estimator=weighted_varaint,
                     orient="v")
    g1.set_axis_labels("Passage", "Variant Frequency {} CI=95%".format(plus_minus))
    g1.set_xticklabels(fontsize=9, rotation=45)
    g1.set(yscale='log')
    g1.set(ylim=(10**-5, 10**-1))

    # plt.show()
    g1.savefig(output_dir + "/All_Mutations_point_plot", dpi=300)
    plt.close()

    data_filter["passage"] = data_filter["passage"].astype(str)
    data_filter["passage"] = "p" + data_filter["passage"]
    g2 = sns.catplot("passage", "frac_and_weight", data=data_filter, hue="Mutation", order=passage_order, palette=mutation_palette(4)
                        ,kind="point", dodge=True, hue_order=transition_order, join=False, estimator=weighted_varaint,
                     orient="v")
    g2.set_axis_labels("Passage", "Variant Frequency {} CI=95%".format(plus_minus))
    g2.set(yscale='log')
    g2.set(ylim=(10 ** -6, 10 ** -2))
    # g2.set_xticklabels(fontsize=10, rotation=45)
    # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/Prgress reports/20200913 Final report/plots" +
    #                   "/Transition_Mutations_point_plot_Mahoney", dpi=300)
    g2.savefig(output_dir + "/Transition_Mutations_point_plot_Mahoney", dpi=300)
    plt.close()

    # data_filter["passage"] = data_filter["passage"].astype(int)
    #
    # linear_reg(data_filter, output_dir, transition_order, type_order, virus="PV1", replica=1, cu=None)
    #
    # # A>G Prev Context
    # data_filter_ag = data_filter[data_filter["Mutation"] == "A>G"]
    #
    # data_filter_ag['Prev'].replace('AA', 'ApA', inplace=True)
    # data_filter_ag['Prev'].replace('UA', 'UpA', inplace=True)
    # data_filter_ag['Prev'].replace('CA', 'CpA', inplace=True)
    # data_filter_ag['Prev'].replace('GA', 'GpA', inplace=True)
    #
    # context_order = ["UpA", "ApA", "CpA", "GpA"]
    # type_order = ["Synonymous", "Non-Synonymous"]
    # data_filter_ag["ADAR_like"] = data_filter_ag.Prev.str.contains('UpA') | data_filter_ag.Prev.str.contains('ApA')
    # data_filter_ag_pass8 = data_filter_ag.loc[data_filter_ag.passage == 8]
    # data_filter_ag_pass8 = data_filter_ag_pass8.loc[data_filter_ag_pass8.Type == "Synonymous"]
    # print(data_filter_ag_pass8.to_string())

    # g5 = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=label_order, palette=mutation_palette(2),
    #             kind = "point", dodge=True, hue_order=[True, False], estimator=weighted_varaint, orient="v",
    #                  col="Type", join=False, col_order=type_order)
    # g5.set_axis_labels("", "Variant Frequency")
    # g5.set(yscale='log')
    # g5.set(ylim=(7*10**-7, 4*10**-3))
    # g5.set_xticklabels(rotation=45)
    # # plt.show()
    # g5.savefig(output_dir + "/Context_point_plot", dpi=300)
    # plt.close()
    #

    # ax = sns.boxplot("ADAR_like", "Frequency", data=data_filter_ag_pass8, palette=mutation_palette(2), order=[True, False])
    # ax = sns.stripplot("ADAR_like", "Frequency", data=data_filter_ag_pass8, color=".2", order=[True, False])
    # old_statannot.add_stat_annotation(ax, data=data_filter_ag_pass8, x="ADAR_like", y="Frequency",
    #                     boxPairList=[(True, False)], test='Mann-Whitney', textFormat='star', loc='inside', verbose=1)
    # ax.set_yscale('log')
    # ax.set_xlabel("ADAR-like\nContext")
    # ax.set_ylabel("Variant Frequency")
    # ax.set(ylim=(10 ** -4, 10 ** -2))
    # # sns.despine()
    # # plt.tight_layout()
    # plt.savefig(output_dir + "/context_p8_point_plot_v2.png", dpi=300)
    # plt.close()

    data_filter_synonymous = data_filter.loc[data_filter.Type == "Synonymous"]
    # data_filter_synonymous["ADAR_like"] = (data_filter_synonymous.Prev.str.contains('UpA') | data_filter_synonymous.Prev.str.contains('ApA'))
    data_filter_synonymous["Mutation"] = np.where(((data_filter_synonymous["Mutation"] == "A>G") &
                                                   (data_filter_synonymous["5`_ADAR_Preference"] == "High")),
                                                  "High\nADAR-like\nA>G", np.where(((data_filter_synonymous["Mutation"] == "A>G")
                                                                                    & (data_filter_synonymous["5`_ADAR_Preference"] == "Intermediate")),
                                                                                   "Intermediate\nADAR-like\nA>G",
                                                                                   np.where(((data_filter_synonymous["Mutation"] == "A>G") &
                                                                                             (data_filter_synonymous["5`_ADAR_Preference"] == "Low")),
                                                                                            "Low\nADAR-like\nA>G",
                                                                                            data_filter_synonymous["Mutation"])))
    data_filter_synonymous["Mutation_adar"] = np.where(((data_filter_synonymous["Mutation"] == "U>C") &
                                                        (data_filter_synonymous["3`_ADAR_Preference"] == "High")),
                                                       "High\nADAR-like\nU>C", np.where(((data_filter_synonymous["Mutation"] == "U>C")
                                                                                         & (data_filter_synonymous["3`_ADAR_Preference"] == "Intermediate")),
                                                                                        "Intermediate\nADAR-like\nU>C",
                                                                                        np.where(((data_filter_synonymous["Mutation"] == "U>C") &
                                                                                                  (data_filter_synonymous["3`_ADAR_Preference"] == "Low")),
                                                                                                 "Low\nADAR-like\nU>C",
                                                                                                 data_filter_synonymous["Mutation"])))
    mutation_adar_order = ["High\nADAR-like\nA>G", "Low\nADAR-like\nA>G",
                           "High\nADAR-like\nU>C", "Low\nADAR-like\nU>C"]

    data_filter_synonymous["passage"] = data_filter_synonymous["passage"].astype(str)
    catplot_adar = sns.catplot(x="passage", y="frac_and_weight", data=data_filter_synonymous, hue="Mutation_adar",
                               order=passage_order, palette=mutation_palette(4, adar=True), kind="point", dodge=True,
                               hue_order=mutation_adar_order, join=False, estimator=weighted_varaint, orient="v",
                               legend=True)
    catplot_adar.set_axis_labels("Passage", "Variant Frequency {} CI=95%".format(plus_minus))
    catplot_adar.set(yscale='log')
    catplot_adar.set(ylim=(10 ** -6, 10 ** -2))
    # catplot_adar.set_xticklabels(fontsize=8)
    # plt.tight_layout()
    plt.savefig(output_dir + "/adar_pref_mutation_point_plot_PV1.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
