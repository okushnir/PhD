
from datetime import datetime
from AccuNGS_analysis.new_analysis_fuctions import *


def main():
    input_dir = "/Users/odedkushnir/PhD_Projects/After_review/AccuNGS/CVB3/"
    date = datetime.today().strftime("%Y%m%d")
    prefix = "inosine_predict_context"
    output_dir = input_dir + prefix
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    q_file_name = "/q38_data_mutation.csv"
    data_adar = pd.read_csv("/Users/odedkushnir/PhD_Projects/After_review/AccuNGS/CVB3/InosinePredict_reults/Output/CVB3_adar1_trans.csv")
    columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
               "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts",
               "Consensus>Mutated_codon", "fiveGrade", "threeGrade"]
    removed_mutation = None
    virus = "CVB3"
    data_filter = analysis(input_dir, output_dir, q_file_name, data_adar, columns, virus, removed_mutation, filter_reads=True)

    """Plots"""
    passage_order = ["RNA\nControl", "p2", "p5", "p8", "p10", "p12"]
    transition_order = ["A>G", "U>C", "G>A", "C>U"]
    pairs = [(("RNA\nControl", "A>G"), ("RNA\nControl", "G>A")), (("p2", "A>G"), ("p2", "G>A")),
             (("p5", "A>G"), ("p5", "G>A")), (("p8", "A>G"), ("p8", "G>A")),
             (("p10", "A>G"), ("p10", "G>A")), (("p12", "A>G"), ("p12", "G>A")),
             (("RNA\nControl", "A>G"), ("RNA\nControl", "U>C")), (("p2", "A>G"), ("p2", "U>C")),
             (("p5", "A>G"), ("p5", "U>C")), (("p8", "A>G"), ("p8", "U>C")),
             (("p10", "A>G"), ("p10", "U>C")), (("p12", "A>G"), ("p12", "U>C")),
             (("RNA\nControl", "A>G"), ("RNA\nControl", "C>U")), (("p2", "A>G"), ("p2", "C>U")),
             (("p5", "A>G"), ("p5", "C>U")), (("p8", "A>G"), ("p8", "C>U")),
             (("p10", "A>G"), ("p10", "C>U")), (("p12", "A>G"), ("p12", "C>U"))]

    pairs_adar = [(("RNA\nControl", "High\nADAR-like\nA>G"), ("RNA\nControl", "Low\nADAR-like\nA>G")),
                  (("p2", "High\nADAR-like\nA>G"), ("p2", "Low\nADAR-like\nA>G")),
                  (("p5", "High\nADAR-like\nA>G"), ("p5", "Low\nADAR-like\nA>G")),
                  (("p8", "High\nADAR-like\nA>G"), ("p8", "Low\nADAR-like\nA>G")),
                  (("p10", "High\nADAR-like\nA>G"), ("p10", "Low\nADAR-like\nA>G")),
                  (("p12", "High\nADAR-like\nA>G"), ("p12", "Low\nADAR-like\nA>G")),
                  (("RNA\nControl", "High\nADAR-like\nU>C"), ("RNA\nControl", "Low\nADAR-like\nU>C")),
                  (("p2", "High\nADAR-like\nU>C"), ("p2", "Low\nADAR-like\nU>C")),
                  (("p5", "High\nADAR-like\nU>C"), ("p5", "Low\nADAR-like\nU>C")),
                  (("p8", "High\nADAR-like\nU>C"), ("p8", "Low\nADAR-like\nU>C")),
                  (("p10", "High\nADAR-like\nU>C"), ("p10", "Low\nADAR-like\nU>C")),
                  (("p12", "High\nADAR-like\nU>C"), ("p12", "Low\nADAR-like\nU>C"))]
    label_order = ["CVB3\nRNA Control", "CVB3-p2", "CVB3-p5", "CVB3-p8", "CVB3-p10", "CVB3-p12"]
    plots(input_dir, date, data_filter, virus, passage_order, transition_order, pairs, label_order, pairs_adar
          , filter_reads=True)


if __name__ == "__main__":
    main()

    # import pandas as pd
    # import os
    # import matplotlib.pyplot as plt
    # from FITS_analysis import fits_plotter
    # import numpy as np
    # from matplotlib.ticker import ScalarFormatter
    # import matplotlib.ticker as ticker
    # import seaborn as sns
    # from FITS_analysis import fits_new_plotter
    # from AccuNGS_analysis.adar_mutation_palette import mutation_palette
    # from AccuNGS_analysis.Linear_regression import linear_reg
    # from scipy import stats
    # from statannotations.Annotator import Annotator
    # import contextlib
    #
    # sns.set(font_scale=1.2)
    # sns.set_style("ticks")
    # sns.despine()
    #
    #
    # def weighted_varaint(x, **kws):
    #     var, count = map(np.asarray, zip(*x))
    #     return var.sum() / count.sum()
    #
    # # flatui = ["#3498db", "#9b59b6"]
    # # date = "20211211"
    # # input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/CVB3"
    # input_dir = "/Users/odedkushnir/PhD_Projects/After_review/AccuNGS/CVB3/"
    # date = datetime.today().strftime("%Y%m%d")
    # print(date)
    # prefix = "inosine_predict_context"
    # output_dir = input_dir + prefix
    # try:
    #     os.mkdir(output_dir)
    # except OSError:
    #     print("Creation of the directory %s failed" % output_dir)
    # else:
    #     print("Successfully created the directory %s " % output_dir)
    #
    #
    # data_mutations = pd.read_csv(input_dir + "/q38_data_mutation.csv")
    # data_adar = pd.read_csv("/Users/odedkushnir/PhD_Projects/After_review/AccuNGS/CVB3/InosinePredict_reults/Output/CVB3_adar1_trans.csv")
    # data_mutations = data_mutations.merge(data_adar, on="Pos", how="inner")
    #
    # columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
    #            "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts",
    #            "Consensus>Mutated_codon", "fiveGrade", "threeGrade"]
    # data_filter = pd.DataFrame(data_mutations, columns=columns)
    # data_filter["pval"] = data_filter["pval"].fillna(1)
    # data_filter["no_variants"] = data_filter["Frequency"] * data_filter["Read_count"]
    # # filter based on pval<0.01 and Prob>0.95
    # # data_filter["no_variants"] = np.where(data_filter["pval"] > 0.01, 0, data_filter["no_variants"])
    # data_filter["no_variants"] = np.where(data_filter["Prob"] < 0.95, 0, data_filter["no_variants"])
    #
    # data_filter["frac_and_weight"] = list(zip(data_filter.no_variants, data_filter.Read_count))
    # data_filter["passage"] = data_filter["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    # data_filter["passage"] = np.where(data_filter["passage"] == "RNA Control", 0, data_filter["passage"])
    # data_filter["passage"] = data_filter["passage"].astype(int)
    # data_filter["Type"] = data_filter["Type"].fillna("NonCodingRegion")
    # # data_filter = data_filter.loc[data_filter.Mutation != "C>U"]
    #
    # data_filter_ag = data_filter[data_filter["Mutation"] == "A>G"]
    # data_filter_uc = data_filter[data_filter["Mutation"] == "U>C"]
    #
    # print("25_quantile_ag %s" % str(data_filter_ag["fiveGrade"].quantile(0.25)))
    # print("75_quantile_ag %s" % str(data_filter_ag["fiveGrade"].quantile(0.75)))
    # print("25_quantile_uc %s" % str(data_filter_uc["threeGrade"].quantile(0.25)))
    # print("75_quantile_uc %s" % str(data_filter_uc["threeGrade"].quantile(0.75)))
    #
    # data_filter["ADAR_grade_five"] = np.where(data_filter["fiveGrade"] < data_filter_ag["fiveGrade"].quantile(0.25), 0,
    #                                           np.where(data_filter["fiveGrade"] <= data_filter_ag["fiveGrade"].
    #                                                    quantile(0.75), 0.5, 1))
    # data_filter["5`_ADAR_Preference"] = np.where(data_filter["fiveGrade"] < data_filter_ag["fiveGrade"].quantile(0.25),
    #                                              "Low", np.where(data_filter["fiveGrade"] <=
    #                                                              data_filter_ag["fiveGrade"].quantile(0.75),
    #                                                              "Intermediate", "High"))
    # data_filter["ADAR_grade_five"] = np.where(data_filter["fiveGrade"] == 0, "0", data_filter["ADAR_grade_five"])
    # data_filter["5`_ADAR_Preference"] = np.where(data_filter["fiveGrade"] == 0, "Low", data_filter["5`_ADAR_Preference"])
    #
    # data_filter["ADAR_grade_three"] = np.where(data_filter["threeGrade"] < data_filter_uc["threeGrade"].quantile(0.25), 0,
    #                                            np.where(data_filter["threeGrade"] <= data_filter_uc["threeGrade"].
    #                                                     quantile(0.75), 0.5, 1))
    # data_filter["3`_ADAR_Preference"] = np.where(data_filter["threeGrade"] < data_filter_uc["threeGrade"].quantile(0.25),
    #                                              "Low", np.where(data_filter["threeGrade"] <=
    #                                                              data_filter_uc["threeGrade"].quantile(0.75),
    #                                                              "Intermediate", "High"))
    # data_filter["ADAR_grade_three"] = np.where(data_filter["threeGrade"] == 0, "0", data_filter["ADAR_grade_five"])
    # data_filter["3`_ADAR_Preference"] = np.where(data_filter["threeGrade"] == 0, "Low", data_filter["3`_ADAR_Preference"])
    #
    # data_filter_ag = data_filter[data_filter["Mutation"] == "A>G"]
    # data_filter_ag['Prev'].replace('AA', 'ApA', inplace=True)
    # data_filter_ag['Prev'].replace('UA', 'UpA', inplace=True)
    # data_filter_ag['Prev'].replace('CA', 'CpA', inplace=True)
    # data_filter_ag['Prev'].replace('GA', 'GpA', inplace=True)
    # data_filter_ag["ADAR_like"] = data_filter_ag.Prev.str.contains('UpA') | data_filter_ag.Prev.str.contains('ApA')
    #
    # data_filter_uc = data_filter[data_filter["Mutation"] == "U>C"]
    # data_filter_uc['Next'].replace('AA', 'ApA', inplace=True)
    # data_filter_uc['Next'].replace('UA', 'UpA', inplace=True)
    # data_filter_uc['Next'].replace('CA', 'CpA', inplace=True)
    # data_filter_uc['Next'].replace('GA', 'GpA', inplace=True)
    # data_filter_uc["ADAR_like"] = data_filter_uc.Next.str.contains('UpA') | data_filter_uc.Next.str.contains('ApA')
    #
    # data_filter.to_csv(output_dir + "/data_filter.csv", sep=',', encoding='utf-8')
    # data_filter_ag.to_csv(output_dir + "/data_filter_ag.csv", sep=',', encoding='utf-8')
    # data_filter_uc.to_csv(output_dir + "/data_filter_uc.csv", sep=',', encoding='utf-8')
    # data_filter.to_pickle(output_dir + "/data_filter.pkl")
    # data_filter_ag.to_pickle(output_dir + "/data_filter_ag.pkl")
    # data_filter_uc.to_pickle(output_dir + "/data_filter_uc.pkl")
    #
    # # columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
    # #            "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts"]
    # # data_filter = pd.DataFrame(data_mutations, columns=columns)
    # # data_filter["pval"] = data_filter["pval"].fillna(1)
    # # data_filter["no_variants"] = data_filter["Frequency"] * data_filter["Read_count"]
    # # """filter based on pval<0.01 and Prob>0.95"""
    # # # data_filter["no_variants"] = np.where(data_filter["pval"] > 0.01, 0, data_filter["no_variants"])
    # # # data_filter["no_variants"] = np.where(data_filter["Prob"] < 0.95, 0, data_filter["no_variants"])
    # # # data_filter["Read_count"] = data_filter[data_filter["Read_count"] > 10000]
    # # data_filter["label"] = np.where(data_filter["label"] == "CVB3-RNA Control", "CVB3\nRNA Control", data_filter["label"])
    # #
    # # data_filter["frac_and_weight"] = list(zip(data_filter.no_variants, data_filter.Read_count))
    # #
    # # data_filter["passage"] = data_filter["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    # # data_filter["passage"] = np.where(data_filter["passage"] == "CVB3\nRNA Control", 0, data_filter["passage"])
    # # data_filter["passage"] = data_filter["passage"].astype(int)
    # # data_filter["Type"] = data_filter["Type"].fillna("NonCodingRegion")
    # # data_filter.to_csv(input_dir + "/data_filter.csv", sep=',', encoding='utf-8')
    #
    #
    # """Plots"""
    # output_dir = input_dir + date + "_plots"
    # plus_minus = u"\u00B1"
    # try:
    #     os.mkdir(output_dir)
    # except OSError:
    #     print("Creation of the directory %s failed" % output_dir)
    # else:
    #     print("Successfully created the directory %s " % output_dir)
    # label_order = ["CVB3\nRNA Control", "CVB3-p2", "CVB3-p5", "CVB3-p8", "CVB3-p10", "CVB3-p12"]
    # passage_order = ["RNA\nControl", "p2", "p5", "p8", "p10", "p12"]
    # mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U"]
    # transition_order = ["A>G", "U>C", "G>A", "C>U"]
    # type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]
    #
    # # g1 = sns.catplot("label", "frac_and_weight", data=data_filter, hue="Mutation", order=label_order, palette="tab20",
    # #                     kind="point", hue_order=mutation_order, join=False, estimator=weighted_varaint, orient="v", dodge=True)
    # # g1.set_axis_labels("", "Variant Frequency {} CI=95%".format(plus_minus))
    # # g1.set_xticklabels(fontsize=9, rotation=45)
    # # g1.set(yscale='log')
    # # g1.set(ylim=(10 ** -7, 10 ** -3))
    # # # plt.show()
    # # g1.savefig(output_dir + "/All_Mutations_point_plot", dpi=300)
    # # plt.close()
    # data_filter["passage"] = data_filter["passage"].astype(str)
    # data_filter["passage"] = "p" + data_filter["passage"]
    # data_filter["passage"] = np.where(data_filter["passage"] == "p0", "RNA\nControl", data_filter["passage"])
    # g2 = sns.catplot("passage", "frac_and_weight", data=data_filter, hue="Mutation", order=passage_order, palette=mutation_palette(4),
    #                     kind="point", hue_order=transition_order, join=False, estimator=weighted_varaint, orient="v",
    #                  dodge=True, legend=True)
    # g2.set_axis_labels("Passage", "Variant Frequency {} CI=95%".format(plus_minus))
    # g2.set(yscale='log')
    # g2.set(ylim=(10 ** -6, 10 ** -2))
    # # g2.set_yticklabels(fontsize=12)
    # # g2.set_xticklabels(fontsize=10, rotation=45)
    # # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/Prgress reports/20200913 Final report/plots" +
    # #                   "/Transition_Mutations_point_plot_CVB3", dpi=300)
    # g2.savefig(output_dir + "/Transition_Mutations_point_plot_CVB3", dpi=300)
    # # g2.savefig(output_dir + "/Transition_Mutations_point_plot", dpi=300)
    # plt.close()
    #
    # pairs = [(("RNA\nControl", "A>G"), ("RNA\nControl", "G>A")), (("p2", "A>G"), ("p2", "G>A")),
    #          (("p5", "A>G"), ("p5", "G>A")), (("p8", "A>G"), ("p8", "G>A")),
    #          (("p10", "A>G"), ("p10", "G>A")), (("p12", "A>G"), ("p12", "G>A")),
    #          (("RNA\nControl", "A>G"), ("RNA\nControl", "U>C")), (("p2", "A>G"), ("p2", "U>C")),
    #          (("p5", "A>G"), ("p5", "U>C")), (("p8", "A>G"), ("p8", "U>C")),
    #          (("p10", "A>G"), ("p10", "U>C")), (("p12", "A>G"), ("p12", "U>C")),
    #          (("RNA\nControl", "A>G"), ("RNA\nControl", "C>U")), (("p2", "A>G"), ("p2", "C>U")),
    #          (("p5", "A>G"), ("p5", "C>U")), (("p8", "A>G"), ("p8", "C>U")),
    #          (("p10", "A>G"), ("p10", "C>U")), (("p12", "A>G"), ("p12", "C>U"))]
    # passage_g = sns.boxplot(x="passage", y="Frequency", data=data_filter, hue="Mutation", order=passage_order,
    #                         palette=mutation_palette(4), dodge=True, hue_order=transition_order)
    # passage_g.set_yscale('log')
    # passage_g.set_ylim(10 ** -6, 10 ** -2)
    #
    # annot = Annotator(passage_g, pairs, x="passage", y="Frequency", hue="Mutation", data=data_filter, hue_order=transition_order)
    # annot.configure(test='t-test_welch', text_format='star', loc='outside', verbose=2, comparisons_correction="Bonferroni") #"Wilcoxon test"
    # annot.apply_test()#alternative="less"
    # file_path = output_dir + "/sts.csv"
    # with open(file_path, "w") as o:
    #     with contextlib.redirect_stdout(o):
    #         passage_g, test_results = annot.annotate()
    # plt.legend(bbox_to_anchor=(1.05, 0.5), loc=2, borderaxespad=0.)
    # plt.tight_layout()
    # plt.savefig(output_dir + "/Transition_Mutations_box_stat_plot_CVB3", dpi=300)
    # plt.close()
    #
    # data_filter_synonymous = data_filter.loc[data_filter.Type == "Synonymous"]
    # # data_filter_synonymous["ADAR_like"] = (data_filter_synonymous.Prev.str.contains('UpA') | data_filter_synonymous.Prev.str.contains('ApA'))
    # data_filter_synonymous["Mutation"] = np.where(((data_filter_synonymous["Mutation"] == "A>G") &
    #                                                (data_filter_synonymous["5`_ADAR_Preference"] == "High")),
    #                                               "High\nADAR-like\nA>G", np.where(((data_filter_synonymous["Mutation"] == "A>G")
    #                                                                                 & (data_filter_synonymous["5`_ADAR_Preference"] == "Intermediate")),
    #                                                                                "Intermediate\nADAR-like\nA>G",
    #                                                                                np.where(((data_filter_synonymous["Mutation"] == "A>G") &
    #                                                                                          (data_filter_synonymous["5`_ADAR_Preference"] == "Low")),
    #                                                                                         "Low\nADAR-like\nA>G",
    #                                                                                         data_filter_synonymous["Mutation"])))
    # data_filter_synonymous["Mutation_adar"] = np.where(((data_filter_synonymous["Mutation"] == "U>C") &
    #                                                     (data_filter_synonymous["3`_ADAR_Preference"] == "High")),
    #                                                    "High\nADAR-like\nU>C", np.where(((data_filter_synonymous["Mutation"] == "U>C")
    #                                                                                      & (data_filter_synonymous["3`_ADAR_Preference"] == "Intermediate")),
    #                                                                                     "Intermediate\nADAR-like\nU>C",
    #                                                                                     np.where(((data_filter_synonymous["Mutation"] == "U>C") &
    #                                                                                               (data_filter_synonymous["3`_ADAR_Preference"] == "Low")),
    #                                                                                              "Low\nADAR-like\nU>C",
    #                                                                                              data_filter_synonymous["Mutation"])))
    # mutation_adar_order = ["High\nADAR-like\nA>G", "Low\nADAR-like\nA>G",
    #                        "High\nADAR-like\nU>C", "Low\nADAR-like\nU>C"]
    #
    # data_filter_synonymous["passage"] = data_filter_synonymous["passage"].astype(str)
    # catplot_adar = sns.catplot(x="passage", y="frac_and_weight", data=data_filter_synonymous, hue="Mutation_adar",
    #                            order=passage_order, palette=mutation_palette(4, adar=True), kind="point", dodge=True,
    #                            hue_order=mutation_adar_order, join=False, estimator=weighted_varaint, orient="v",
    #                            legend=True)
    # catplot_adar.set_axis_labels("Passage", "Variant Frequency {} CI=95%".format(plus_minus))
    # catplot_adar.set(yscale='log')
    # catplot_adar.set(ylim=(10 ** -6, 10 ** -2))
    # # catplot_adar.set_xticklabels(fontsize=8)
    # # plt.tight_layout()
    # plt.savefig(output_dir + "/adar_pref_mutation_point_plot_CVB3.png", dpi=300)
    # plt.close()
    #
    # adar_g = sns.boxplot(x="passage", y="Frequency", data=data_filter_synonymous, hue="Mutation_adar",
    #                      order=passage_order, palette=mutation_palette(4, adar=True), dodge=True,
    #                      hue_order=mutation_adar_order)
    # adar_g.set_yscale('log')
    # adar_g.set_ylim(10 ** -6, 10 ** -2)
    # annot = Annotator(adar_g, pairs, x="passage", y="Frequency", hue="Mutation_adar",
    #                   data=data_filter_synonymous, hue_order=mutation_adar_order, order=passage_order)
    # annot.configure(test='t-test_welch', text_format='star', loc='outside', verbose=2,
    #                 comparisons_correction="Bonferroni")
    # annot.apply_test()
    # file_path = output_dir + "/sts_adar.csv"
    # with open(file_path, "w") as o:
    #     with contextlib.redirect_stdout(o):
    #         adar_g, test_results = annot.annotate()
    # plt.legend(bbox_to_anchor=(1.05, 0.5), loc=2, borderaxespad=0.)
    # plt.tight_layout()
    # plt.savefig(output_dir + "/adar_pref_mutation_box_plot_CVB3.png", dpi=300)
    # plt.close()
    #
    # data_filter["passage"] = data_filter["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    # data_filter["passage"] = np.where(data_filter["passage"] == "RNA Control", 0, data_filter["passage"])
    # data_filter["passage"] = data_filter["passage"].astype(int)
    # data_filter["replica"] = 1
    # linear_reg(data_filter, output_dir, transition_order, type_order, virus="CVB3", replica=1)
    #
    # # data_filter_grouped = data_filter.groupby(["label", "passage", "Type", "Mutation"])[
    # #     "frac_and_weight"].agg(
    # #     lambda x: weighted_varaint(x))
    # # data_filter_grouped = data_filter_grouped.reset_index()
    # # data_filter_grouped = data_filter_grouped.rename(columns={"frac_and_weight": "Frequency"})
    # # data_filter_grouped["Frequency"] = data_filter_grouped["Frequency"].astype(float)
    # # data_filter_grouped = data_filter_grouped[data_filter_grouped["label"] != "RNA Control\nPrimer ID"]
    # # data_filter_grouped = data_filter_grouped[data_filter_grouped["label"] != "RNA Control\nRND"]
    # # # data_filter_grouped = data_filter_grouped[data_filter_grouped["replica"] == replica]
    # # print(data_filter_grouped.to_string())
    # #
    # # data_reg_ag = data_filter_grouped[data_filter_grouped["Mutation"] == "A>G"]
    # # data_reg_uc = data_filter_grouped[data_filter_grouped["Mutation"] == "U>C"]
    # # data_reg_ga = data_filter_grouped[data_filter_grouped["Mutation"] == "G>A"]
    # # data_reg_cu = data_filter_grouped[data_filter_grouped["Mutation"] == "C>U"]
    # #
    # # data_reg_ag_syn = data_reg_ag[data_reg_ag["Type"] == "Synonymous"]
    # # data_reg_uc_syn = data_reg_uc[data_reg_uc["Type"] == "Synonymous"]
    # # data_reg_ga_syn = data_reg_ga[data_reg_ga["Type"] == "Synonymous"]
    # # data_reg_cu_syn = data_reg_cu[data_reg_cu["Type"] == "Synonymous"]
    # #
    # # data_reg_ag_non_syn = data_reg_ag[data_reg_ag["Type"] == "Non-Synonymous"]
    # # data_reg_uc_non_syn = data_reg_uc[data_reg_uc["Type"] == "Non-Synonymous"]
    # # data_reg_ga_non_syn = data_reg_ga[data_reg_ga["Type"] == "Non-Synonymous"]
    # # data_reg_cu_non_syn = data_reg_cu[data_reg_cu["Type"] == "Non-Synonymous"]
    # #
    # # data_reg_ga_pmsc = data_reg_ga[data_reg_ga["Type"] == "Premature Stop Codon"]
    # # data_reg_cu_pmsc = data_reg_cu[data_reg_cu["Type"] == "Premature Stop Codon"]
    # #
    # # stat_slope1, stat_intercept1, r_value1, p_value1, std_err1 = stats.linregress(data_reg_ag_syn['passage'],
    # #                                                                               data_reg_ag_syn
    # #                                                                               ['Frequency'])
    # # stat_slope2, stat_intercept2, r_value2, p_value2, std_err2 = stats.linregress(data_reg_uc_syn['passage'],
    # #                                                                               data_reg_uc_syn[
    # #                                                                                   'Frequency'])
    # # stat_slope3, stat_intercept3, r_value3, p_value3, std_err3 = stats.linregress(data_reg_ga_syn['passage'],
    # #                                                                               data_reg_ga_syn[
    # #                                                                                   'Frequency'])
    # # stat_slope4, stat_intercept4, r_value4, p_value4, std_err4 = stats.linregress(data_reg_cu_syn['passage'],
    # #                                                                               data_reg_cu_syn[
    # #                                                                                   'Frequency'])
    # # # data_reg_adar_nonsyn
    # # stat_slope5, stat_intercept5, r_value5, p_value5, std_err5 = stats.linregress(data_reg_ag_non_syn['passage'],
    # #                                                                               data_reg_ag_non_syn[
    # #                                                                                   'Frequency'])
    # # stat_slope6, stat_intercept6, r_value6, p_value6, std_err6 = stats.linregress(data_reg_uc_non_syn['passage'],
    # #                                                                               data_reg_uc_non_syn[
    # #                                                                                   'Frequency'])
    # # stat_slope7, stat_intercept7, r_value7, p_value7, std_err7 = stats.linregress(data_reg_ga_non_syn['passage'],
    # #                                                                               data_reg_ga_non_syn[
    # #                                                                                   'Frequency'])
    # # stat_slope8, stat_intercept8, r_value8, p_value8, std_err8 = stats.linregress(data_reg_cu_non_syn['passage'],
    # #                                                                               data_reg_cu_non_syn[
    # #                                                                                   'Frequency'])
    # # # pmsc
    # # stat_slope11, stat_intercept11, r_value11, p_value11, std_err11 = stats.linregress(data_reg_ga_pmsc['passage'],
    # #                                                                                    data_reg_ga_pmsc[
    # #                                                                                        'Frequency'])
    # # stat_slope12, stat_intercept12, r_value12, p_value12, std_err12 = stats.linregress(data_reg_cu_pmsc['passage'],
    # #                                                                                    data_reg_cu_pmsc[
    # #                                                                                        'Frequency'])
    # # data_filter_grouped = data_filter_grouped.rename(columns={"passage": "Passage"})
    # # reg_plot = sns.lmplot(x="Passage", y="Frequency", data=data_filter_grouped, hue="Mutation",
    # #                       hue_order=transition_order, fit_reg=True, col="Mutation",
    # #                       col_order=transition_order, row="Type", row_order=type_order, palette=mutation_palette(4),
    # #                       line_kws={'label': "Linear Reg"}, legend=True, height=6)  # markers=["o", "v", "x"]
    # # reg_plot.fig.subplots_adjust(wspace=.02)
    # # ax = reg_plot.axes[0, 0]
    # # ax.legend()
    # # leg = ax.get_legend()
    # # leg._loc = 2
    # # L_labels = leg.get_texts()
    # # label_line_1 = "y={0:.3g}x+{1:.3g}".format(stat_slope1, stat_intercept1)
    # # label_line_2 = "y={0:.3g}x+{1:.3g}".format(stat_slope2, stat_intercept2)
    # # label_line_3 = "y={0:.3g}x+{1:.3g}".format(stat_slope3, stat_intercept3)
    # # label_line_4 = "y={0:.3g}x+{1:.3g}".format(stat_slope4, stat_intercept4)
    # # label_line_5 = "y={0:.3g}x+{1:.3g}".format(stat_slope5, stat_intercept5)
    # # label_line_6 = "y={0:.3g}x+{1:.3g}".format(stat_slope6, stat_intercept6)
    # # label_line_7 = "y={0:.3g}x+{1:.3g}".format(stat_slope7, stat_intercept7)
    # # label_line_8 = "y={0:.3g}x+{1:.3g}".format(stat_slope8, stat_intercept8)
    # # label_line_11 = "y={0:.3g}x+{1:.3g}".format(stat_slope11, stat_intercept11)
    # # label_line_12 = "y={0:.3g}x+{1:.3g}".format(stat_slope12, stat_intercept12)
    # # L_labels[0].set_text(label_line_1)
    # # ax = reg_plot.axes[0, 1]
    # # ax.legend()
    # # leg = ax.get_legend()
    # # leg._loc = 2
    # # L_labels = leg.get_texts()
    # # L_labels[0].set_text(label_line_2)
    # # ax = reg_plot.axes[0, 2]
    # # ax.legend()
    # # leg = ax.get_legend()
    # # leg._loc = 2
    # # L_labels = leg.get_texts()
    # # L_labels[0].set_text(label_line_3)
    # # ax = reg_plot.axes[0, 3]
    # # ax.legend()
    # # leg = ax.get_legend()
    # # leg._loc = 2
    # # L_labels = leg.get_texts()
    # # L_labels[0].set_text(label_line_4)
    # # ax = reg_plot.axes[1, 0]
    # # ax.legend()
    # # leg = ax.get_legend()
    # # leg._loc = 2
    # # L_labels = leg.get_texts()
    # # L_labels[0].set_text(label_line_5)
    # # ax = reg_plot.axes[1, 1]
    # # ax.legend()
    # # leg = ax.get_legend()
    # # leg._loc = 2
    # # L_labels = leg.get_texts()
    # # L_labels[0].set_text(label_line_6)
    # # ax = reg_plot.axes[1, 2]
    # # ax.legend()
    # # leg = ax.get_legend()
    # # leg._loc = 2
    # # L_labels = leg.get_texts()
    # # L_labels[0].set_text(label_line_7)
    # # ax = reg_plot.axes[1, 3]
    # # ax.legend()
    # # leg = ax.get_legend()
    # # leg._loc = 2
    # # L_labels = leg.get_texts()
    # # L_labels[0].set_text(label_line_8)
    # # ax = reg_plot.axes[2, 2]
    # # ax.legend()
    # # leg = ax.get_legend()
    # # leg._loc = 2
    # # L_labels = leg.get_texts()
    # # L_labels[0].set_text(label_line_11)
    # # ax = reg_plot.axes[2, 3]
    # # ax.legend()
    # # leg = ax.get_legend()
    # # leg._loc = 2
    # # L_labels = leg.get_texts()
    # # L_labels[0].set_text(label_line_12)
    # #
    # # reg_plot.set(xlim=(0, 13))
    # # reg_plot.set(ylim=(0.000, 0.001))
    # # # reg_plot.fig.suptitle("RV #%s" % str(replica), y=0.99)
    # # # plt.tight_layout()
    # # reg_plot.savefig(output_dir + "/transition_lmplot.png", dpi=300)
    # # plt.close()
    # #
    # # columns = ["Mutation", "Type", "Slope", "Intercept"]
    # # mutation_rate_df = pd.DataFrame(columns=columns)
    # # mutation_rate_df.loc[0] = ["A>G", "Synonymous", stat_slope1, stat_intercept1]
    # # mutation_rate_df.loc[1] = ["U>C", "Synonymous", stat_slope2, stat_intercept2]
    # # mutation_rate_df.loc[2] = ["G>A", "Synonymous", stat_slope3, stat_intercept3]
    # # mutation_rate_df.loc[3] = ["C>U", "Synonymous", stat_slope4, stat_intercept4]
    # # mutation_rate_df.loc[4] = ["A>G", "Non-Synonymous", stat_slope5, stat_intercept5]
    # # mutation_rate_df.loc[5] = ["U>C", "Non-Synonymous", stat_slope6, stat_intercept6]
    # # mutation_rate_df.loc[6] = ["G>A", "Non-Synonymous", stat_slope7, stat_intercept7]
    # # mutation_rate_df.loc[7] = ["C>U", "Non-Synonymous", stat_slope8, stat_intercept8]
    # # mutation_rate_df.loc[8] = ["G>A",  "Pre Mature Stop Codon", stat_slope11, stat_intercept11]
    # # mutation_rate_df.loc[9] = ["C>U",  "Pre Mature Stop Codon", stat_slope12, stat_intercept12]
    # # mutation_rate_df["Virus"] = "CVB3"
    # # mutation_rate_df["Replica"] = 1
    # # mutation_rate_df.to_csv(output_dir + "/mutation_rate.csv", sep=',', encoding='utf-8')
    #
    # data_filter["Transition"] = data_filter.Mutation.str.contains("A>G") | data_filter.Mutation.str.contains("U>C") \
    #                             | data_filter.Mutation.str.contains("C>U") | data_filter.Mutation.str.contains("G>A")
    # data_transition = data_filter[data_filter["Transition"] == True]
    # position_mutation = sns.relplot(x="Pos", y="Frequency", data=data_transition, hue="Mutation",
    #                                 col="passage", hue_order=transition_order, col_wrap=3,
    #                                 palette=mutation_palette(4), height=4)
    # #
    # # data_pmsc = data_filter[data_filter["Type"] == "Premature Stop Codon"]
    # # data_pmsc["mutation_type"] = data_pmsc.Mutation.str.contains("A>G") | data_pmsc.Mutation.str.contains("U>C") | \
    # #                              data_pmsc.Mutation.str.contains("C>U") | data_pmsc.Mutation.str.contains("G>A")
    # # data_pmsc_transition = data_pmsc[data_pmsc["mutation_type"] == True]
    # # g3 = sns.catplot("label", "frac_and_weight", data=data_pmsc_transition, hue="Mutation", order=label_order, palette="tab20",
    # #                     estimator=weighted_varaint, orient="v", dodge=True, kind="point",
    # #                  col="Type", join=False)
    # # g3.set_axis_labels("Sample", "Variant Frequency")
    # # g3.set(yscale='log')
    # # g3.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -5)
    # # g3.set_xticklabels(rotation=45)
    # # # plt.show()
    # # g3.savefig(output_dir + "/Transitions_PMSC_Mutations_point_plot", dpi=300)
    # # plt.close()
    # #
    # # g4 = sns.relplot("passage", "frac_and_weight", data=data_filter, hue="Mutation", palette="tab20",
    # #                     hue_order=transition_order, estimator=weighted_varaint, col="Type", kind="line",
    # #                  col_order=type_order)
    # #
    # # g4.axes.flat[0].set_yscale('symlog', linthreshy=10**-5)
    # # g4.set_axis_labels("Passage", "Variant Frequency")
    # # # plt.show()
    # # g4.savefig(output_dir + "/Time_Transition_Mutations_line_plot", dpi=300)
    # # plt.close()
    # #
    # # # A>G Prev Context
    # #
    # data_filter_ag = data_filter[data_filter["Mutation"] == "A>G"]
    # data_filter_ag['Prev'].replace('AA', 'ApA', inplace=True)
    # data_filter_ag['Prev'].replace('UA', 'UpA', inplace=True)
    # data_filter_ag['Prev'].replace('CA', 'CpA', inplace=True)
    # data_filter_ag['Prev'].replace('GA', 'GpA', inplace=True)
    # data_filter_ag = data_filter_ag.rename(columns={"Prev": "Context"})
    # data_filter_ag["ADAR_like"] = data_filter_ag.Context.str.contains('UpA') | data_filter_ag.Context.str.contains(
    #     'ApA')
    # data_filter_ag.to_pickle(output_dir + "/data_filter_ag.pkl")
    #
    # type_order = ["Synonymous", "Non-Synonymous"]
    #
    # g5 = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=label_order,
    #                  palette=mutation_palette(2), kind="point", dodge=True, hue_order=[True, False], estimator=weighted_varaint,
    #                  orient="v", col="Type", join=False, col_order=type_order)
    # g5.set_axis_labels("", "Variant Frequency {} CI=95%".format(plus_minus))
    # g5.set(yscale='log')
    # g5.set(ylim=(7 * 10 ** -7, 4 * 10 ** -3))
    # g5.set_xticklabels(fontsize=9, rotation=90)
    # # plt.show()
    # g5.savefig(output_dir + "/Context_point_plot", dpi=300)
    # plt.close()
    #
    # data_filter_ag = data_filter_ag[data_filter_ag["passage"] != 0]
    # g6 = sns.relplot("passage", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", palette=mutation_palette(2),
    #                     hue_order=[True, False], estimator=weighted_varaint, col="Type", kind="line",
    #                  col_order=type_order)
    #
    # g6.axes.flat[0].set_yscale('symlog', linthreshy=10**-5)
    # g6.set(ylim=(0, 10 ** -2))
    # yaxis = plt.gca().yaxis
    # yaxis.set_minor_locator(fits_new_plotter.MinorSymLogLocator(1e-1))
    # g6.set_axis_labels("Passage", "Variant Frequency {} CI=95%".format(plus_minus))
    # # plt.show()
    # g6.savefig(output_dir + "/Context_miseq_sample_time_line_plot", dpi=300)
    # plt.close()
    #
    #
    # # data_filter_uc = data_filter[data_filter["Mutation"] == "U>C"]
    # #
    # # data_filter_uc['Next'].replace('UA', 'UpA', inplace=True)
    # # data_filter_uc['Next'].replace('UU', 'UpU', inplace=True)
    # # data_filter_uc['Next'].replace('UC', 'UpC', inplace=True)
    # # data_filter_uc['Next'].replace('UG', 'UpG', inplace=True)
    # #
    # #
    # # data_filter_uc.to_csv(input_dir + "/data_filter_uc.csv", sep=',', encoding='utf-8')
    # # context_order_uc = ["UpA", "UpU", "UpG",  "UpC"]
    # #
    # #
    # # g7 = sns.catplot("label", "frac_and_weight", data=data_filter_uc, hue="Next", order=label_order, palette="tab20",
    # #                     hue_order=context_order_uc, estimator=weighted_varaint, orient="v", dodge=True, kind="point",
    # #                  col="Type", join=False, col_order=type_order)
    # # g7.set_axis_labels("", "Variant Frequency")
    # # g7.set(yscale='log')
    # # g7.set(ylim=(10**-7, 10**-3))
    # # g7.set_xticklabels(rotation=45)
    # # # plt.show()
    # # g7.savefig(output_dir + "/UC_Context_point_plot", dpi=300)
    # # plt.close()
    # #
    # # data_filter_uc = data_filter_uc[data_filter_uc["passage"] != 0]
    # # g8 = sns.relplot("passage", "frac_and_weight", data=data_filter_uc, hue="Next", palette="tab20",
    # #                     hue_order=context_order_uc, estimator=weighted_varaint, col="Type", kind="line",
    # #                  col_order=type_order)
    # #
    # # # g8.set(yscale="log")
    # # g8.axes.flat[0].set_yscale('symlog', linthreshy=10**-4)
    # # # g8.set(ylim=(0, 10 ** -6))
    # # yaxis = plt.gca().yaxis
    # # yaxis.set_minor_locator(MinorSymLogLocator(1e-1))
    # # g8.set_axis_labels("Passage", "Variant Frequency")
    # # # plt.show()
    # # g8.savefig(output_dir + "/UC_Context_miseq_sample_time_line_plot", dpi=300)
    # # plt.close()