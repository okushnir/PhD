
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from AccuNGS_analysis.adar_mutation_palette import mutation_palette
from statannotations.Annotator import Annotator
import contextlib


sns.set(font_scale=1.2)
sns.set_style("ticks")
sns.despine()


def weighted_varaint(x, **kws):
    var, count = map(np.asarray, zip(*x))
    return var.sum() / count.sum()


def analysis(input_dir, output_dir, q_file_name, data_adar, columns, removed_mutation=None, replica=None, filter_reads=None):
    data_mutations = pd.read_csv(input_dir + q_file_name)
    data_mutations = data_mutations[data_mutations["Rank"] != 0]
    data_mutations = data_mutations.merge(data_adar, on="Pos", how="inner")
    data_filter = pd.DataFrame(data_mutations, columns=columns)
    data_filter["pval"] = data_filter["pval"].fillna(1)
    data_filter["no_variants"] = data_filter["Frequency"] * data_filter["Read_count"]
    if filter_reads is True:
        # filter based on pval<0.01 and Prob>0.95
        # data_filter["no_variants"] = np.where(data_filter["pval"] > 0.01, 0, data_filter["no_variants"])
        data_filter["no_variants"] = np.where(data_filter["Prob"] < 0.95, 0, data_filter["no_variants"])
        data_filter["Read_count"] = data_filter[data_filter["Read_count"] > 10000]

    data_filter["frac_and_weight"] = list(zip(data_filter.no_variants, data_filter.Read_count))
    data_filter["passage"] = data_filter["label"].apply(lambda x: x.split("-")[-1][1])
    data_filter["passage"] = data_filter["passage"].astype(int)
    data_filter["Type"] = data_filter["Type"].fillna("NonCodingRegion")
    if removed_mutation is not None:
        data_filter = data_filter.loc[data_filter.Mutation != removed_mutation]
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
    return data_filter


def plots(input_dir, date, data_filter, virus, passage_order, transition_order, pairs, label_order, pairs_adar, filter_reads=None):
    output_dir = input_dir + date + "_plots"
    plus_minus = u"\u00B1"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)
    if filter_reads is True:
        data_filter["no_variants"] = np.where(data_filter["Prob"] < 0.95, 0, data_filter["no_variants"])
        data_filter["Read_count"] = data_filter[data_filter["Read_count"] > 10000]
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U"]
    type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]
    # g1 = sns.catplot("label", "frac_and_weight", data=data_filter, hue="Mutation", order=label_order, palette="tab20",
    #                     kind="point", dodge=True, hue_order=mutation_order, join=False, estimator=weighted_varaint,
    #                  orient="v")
    # g1.set_axis_labels("Passage", "Variant Frequency {} CI=95%".format(plus_minus))
    # g1.set_xticklabels(fontsize=9, rotation=45)
    # g1.set(yscale='log')
    # g1.set(ylim=(10**-5, 10**-1))
    #
    # # plt.show()
    # g1.savefig(output_dir + "/All_Mutations_point_plot", dpi=300)
    # plt.close()
    data_filter["passage"] = data_filter["passage"].astype(str)
    data_filter["passage"] = "p" + data_filter["passage"]
    g2 = sns.catplot("passage", "frac_and_weight", data=data_filter, hue="Mutation", order=passage_order,
                     palette=mutation_palette(4)
                     , kind="point", dodge=0.5, hue_order=transition_order, join=False, estimator=weighted_varaint,
                     orient="v")
    g2.set_axis_labels("Passage", "Variant Frequency {} CI=95%".format(plus_minus))
    g2.set(yscale='log')
    g2.set(ylim=(10 ** -6, 10 ** -2))
    # g2.set_xticklabels(fontsize=10, rotation=45)
    # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/Prgress reports/20200913 Final report/plots" +
    #                   "/Transition_Mutations_point_plot_Mahoney", dpi=300)
    g2.savefig(output_dir + "/Transition_Mutations_point_plot_{0}".format(virus), dpi=300)
    plt.close()

    passage_g = sns.boxplot(x="passage", y="Frequency", data=data_filter, hue="Mutation", order=passage_order,
                            palette=mutation_palette(4), dodge=True, hue_order=transition_order)
    passage_g.set_yscale('log')
    passage_g.set_ylim(10 ** -6, 10 ** -2)

    annot = Annotator(passage_g, pairs, x="passage", y="Frequency", hue="Mutation", data=data_filter,
                      order=passage_order, hue_order=transition_order)
    annot.configure(test='t-test_welch', text_format='star', loc='inside', verbose=2,
                    comparisons_correction="Bonferroni")
    annot.apply_test()
    file_path = output_dir + "/sts.csv"
    with open(file_path, "w") as o:
        with contextlib.redirect_stdout(o):
            passage_g, test_results = annot.annotate()
    plt.legend(bbox_to_anchor=(1.05, 0.5), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(output_dir + "/Transition_Mutations_box_stat_plot_{0}".format(virus), dpi=300)
    plt.close()

    data_filter_synonymous = data_filter.loc[data_filter.Type == "Synonymous"]
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
                               order=passage_order, palette=mutation_palette(4, adar=True), kind="point", dodge=0.5,
                               hue_order=mutation_adar_order, join=False, estimator=weighted_varaint, orient="v",
                               legend=True)
    catplot_adar.set_axis_labels("Passage", "Variant Frequency {0} CI=95%".format(plus_minus))
    catplot_adar.set(yscale='log')
    catplot_adar.set(ylim=(10 ** -6, 10 ** -2))
    plt.savefig(output_dir + "/adar_pref_mutation_point_plot_{0}.png".format(virus), dpi=300)
    plt.close()

    adar_g = sns.boxplot(x="passage", y="Frequency", data=data_filter_synonymous, hue="Mutation_adar",
                         order=passage_order, palette=mutation_palette(4, adar=True), dodge=True,
                         hue_order=mutation_adar_order)
    adar_g.set_yscale('log')
    adar_g.set_ylim(10 ** -6, 10 ** -1)
    adar_g.set(xlabel="Passage", ylabel="Variant Frequency")
    annot = Annotator(adar_g, pairs_adar, x="passage", y="Frequency", hue="Mutation_adar",
                      data=data_filter_synonymous, hue_order=mutation_adar_order, order=passage_order)
    annot.configure(test='t-test_welch', text_format='star', loc='outside', verbose=2,
                    comparisons_correction="Bonferroni")
    annot.apply_test()
    file_path = output_dir + "/sts_adar.csv"
    with open(file_path, "w") as o:
        with contextlib.redirect_stdout(o):
            adar_g, test_results = annot.annotate()
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(output_dir + "/adar_pref_mutation_box_plot_{0}.png".format(virus), dpi=300)
    plt.close()
