
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
from datetime import datetime
from statannotations.Annotator import Annotator
import contextlib
from scikit_posthocs import posthoc_dunn


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
    # input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/"
    """Local"""
    # input_dir = "/Users/odedkushnir/PhD_Projects/After_review/AccuNGS/RV/passages/"
    input_dir = "C:/Users/odedku/PhD_Projects/After_review/AccuNGS/RV/passages/"
    prefix = "inosine_predict_context"
    date = datetime.today().strftime("%Y%m%d")
    output_dir = input_dir + "{0}_{1}".format(date, prefix)
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    data_filter = pd.read_pickle(input_dir + prefix + "/data_filter.pkl")
    data_filter["passage"] = data_filter["passage"].astype(int)
    data_filter["no_variants"] = np.where(data_filter["Prob"] < 0.95, 0, data_filter["no_variants"])
    data_filter = data_filter.loc[data_filter["Read_count"] > 10000]

    replica_lst = [1, 2, 3]
    for replica in replica_lst:
        data_filter_replica = data_filter[data_filter["replica"] == replica]
        data_filter_replica["passage_p"] = data_filter_replica["passage"]
        data_filter_replica["passage_p"] = data_filter_replica["passage_p"].astype(str)
        data_filter_replica["passage_p"] = "p" + data_filter_replica["passage_p"]
        if replica == 2:
            data_filter_replica = pd.read_pickle(input_dir + prefix + "/data_filter.pkl")
            data_filter_replica["passage"] = data_filter_replica["passage"].astype(int)
            data_filter_replica["no_variants"] = np.where(data_filter_replica["Prob"] < 0.95, 0, data_filter_replica["no_variants"])
            data_filter_replica = data_filter_replica.loc[data_filter_replica["Read_count"] > 10000]
            data_filter_replica["passage_p"] = data_filter_replica["passage"]
            data_filter_replica["passage_p"] = data_filter_replica["passage_p"].astype(str)
            data_filter_replica["passage_p"] = "p" + data_filter_replica["passage_p"]
            data_filter_replica["replica"] = np.where(((data_filter_replica["replica"] == 1) &
                                                       (data_filter_replica["passage"] == 0)), 2,
                                                      data_filter_replica["replica"])
            data_filter_replica = data_filter_replica[data_filter_replica["replica"] == replica]
            data_filter_replica.to_csv(output_dir + "/data_filter_replica{0}.csv".format(str(replica)), sep=",")
        data_filter_replica["passage_p"] = np.where(data_filter_replica["passage_p"] == "p0", "RNA Control", data_filter_replica["passage_p"])
        data_filter_replica["Mutation Type"] = np.where(data_filter_replica["Mutation"] == "A>G", "Transitions",
                                                np.where(data_filter_replica["Mutation"] == "U>C", "Transitions",
                                                         np.where(data_filter_replica["Mutation"] == "G>A", "Transitions",
                                                                  np.where(data_filter_replica["Mutation"] == "C>U",
                                                                           "Transitions",
                                                                           np.where(data_filter_replica["Mutation"] ==
                                                                                    "G>U", "Oxidations",
                                                                                    np.where(data_filter_replica["Mutation"] == "A>C", "Oxidations",
                                                                           "Transversions"))))))
        
        if replica == 1:
            passage_p_order = ["RNA Control", "p2", "p5", "p8", "p12"]
            pairs = [(("RNA Control", "A>G"), ("RNA Control", "G>A")), (("p2", "A>G"), ("p2", "G>A")),
                     (("p5", "A>G"), ("p5", "G>A")), (("p8", "A>G"), ("p8", "G>A")),
                     (("p12", "A>G"), ("p12", "G>A")),
                     (("RNA Control", "A>G"), ("RNA Control", "U>C")), (("p2", "A>G"), ("p2", "U>C")),
                     (("p5", "A>G"), ("p5", "U>C")), (("p8", "A>G"), ("p8", "U>C")),
                     (("p12", "A>G"), ("p12", "U>C")),
                     (("RNA Control", "A>G"), ("RNA Control", "C>U")), (("p2", "A>G"), ("p2", "C>U")),
                     (("p5", "A>G"), ("p5", "C>U")), (("p8", "A>G"), ("p8", "C>U")),
                     (("p12", "A>G"), ("p12", "C>U"))]
            pairs_adar = [(("RNA Control", "High ADAR-like A>G"), ("RNA Control", "Low ADAR-like A>G")),
                          (("p2", "High ADAR-like A>G"), ("p2", "Low ADAR-like A>G")),
                          (("p5", "High ADAR-like A>G"), ("p5", "Low ADAR-like A>G")),
                          (("p8", "High ADAR-like A>G"), ("p8", "Low ADAR-like A>G")),
                          (("p12", "High ADAR-like A>G"), ("p12", "Low ADAR-like A>G")),
                          (("RNA Control", "High ADAR-like U>C"), ("RNA Control", "Low ADAR-like U>C")),
                          (("p2", "High ADAR-like U>C"), ("p2", "Low ADAR-like U>C")),
                          (("p5", "High ADAR-like U>C"), ("p5", "Low ADAR-like U>C")),
                          (("p8", "High ADAR-like U>C"), ("p8", "Low ADAR-like U>C")),
                          (("p12", "High ADAR-like U>C"), ("p12", "Low ADAR-like U>C"))]
            trans_pairs = [(("RNA Control", "Transitions"), ("RNA Control", "Transversions")),
                           (("RNA Control", "Transitions"), ("RNA Control", "Oxidations")),
                           (("p2", "Transitions"), ("p2", "Transversions")),
                           (("p2", "Transitions"), ("p2", "Oxidations")),
                           (("p5", "Transitions"), ("p5", "Transversions")),
                           (("p5", "Transitions"), ("p5", "Oxidations")),
                           (("p8", "Transitions"), ("p8", "Transversions")),
                           (("p8", "Transitions"), ("p8", "Oxidations")),
                           (("p12", "Transitions"), ("p12", "Transversions")),
                           (("p12", "Transitions"), ("p12", "Oxidations"))]
        else:
            passage_p_order = ["RNA Control", "p2", "p5", "p8", "p10", "p12"]
            pairs = [(("RNA Control", "A>G"), ("RNA Control", "G>A")), (("p2", "A>G"), ("p2", "G>A")),
                     (("p5", "A>G"), ("p5", "G>A")), (("p8", "A>G"), ("p8", "G>A")),
                     (("p10", "A>G"), ("p10", "G>A")), (("p12", "A>G"), ("p12", "G>A")),
                     (("RNA Control", "A>G"), ("RNA Control", "U>C")), (("p2", "A>G"), ("p2", "U>C")),
                     (("p5", "A>G"), ("p5", "U>C")), (("p8", "A>G"), ("p8", "U>C")),
                     (("p10", "A>G"), ("p10", "U>C")), (("p12", "A>G"), ("p12", "U>C")),
                     (("RNA Control", "A>G"), ("RNA Control", "C>U")), (("p2", "A>G"), ("p2", "C>U")),
                     (("p5", "A>G"), ("p5", "C>U")), (("p8", "A>G"), ("p8", "C>U")),
                     (("p10", "A>G"), ("p10", "C>U")), (("p12", "A>G"), ("p12", "C>U"))]
            pairs_adar = [(("RNA Control", "High ADAR-like A>G"), ("RNA Control", "Low ADAR-like A>G")),
                          (("p2", "High ADAR-like A>G"), ("p2", "Low ADAR-like A>G")),
                          (("p5", "High ADAR-like A>G"), ("p5", "Low ADAR-like A>G")),
                          (("p8", "High ADAR-like A>G"), ("p8", "Low ADAR-like A>G")),
                          (("p10", "High ADAR-like A>G"), ("p10", "Low ADAR-like A>G")),
                          (("p12", "High ADAR-like A>G"), ("p12", "Low ADAR-like A>G")),
                          (("RNA Control", "High ADAR-like U>C"), ("RNA Control", "Low ADAR-like U>C")),
                          (("p2", "High ADAR-like U>C"), ("p2", "Low ADAR-like U>C")),
                          (("p5", "High ADAR-like U>C"), ("p5", "Low ADAR-like U>C")),
                          (("p8", "High ADAR-like U>C"), ("p8", "Low ADAR-like U>C")),
                          (("p10", "High ADAR-like U>C"), ("p10", "Low ADAR-like U>C")),
                          (("p12", "High ADAR-like U>C"), ("p12", "Low ADAR-like U>C"))]
            trans_pairs = [(("RNA Control", "Transitions"), ("RNA Control", "Transversions")),
                           (("RNA Control", "Transitions"), ("RNA Control", "Oxidations")),
                           (("p2", "Transitions"), ("p2", "Transversions")),
                           (("p2", "Transitions"), ("p2", "Oxidations")),
                           (("p5", "Transitions"), ("p5", "Transversions")),
                           (("p5", "Transitions"), ("p5", "Oxidations")),
                           (("p8", "Transitions"), ("p8", "Transversions")),
                           (("p8", "Transitions"), ("p8", "Oxidations")),
                           (("p10", "Transitions"), ("p10", "Transversions")),
                           (("p10", "Transitions"), ("p10", "Oxidations")),
                           (("p12", "Transitions"), ("p12", "Transversions")),
                           (("p12", "Transitions"), ("p12", "Oxidations"))]

        """Plots"""
        label_order = ["RNA Control\nRND", "RNA Control\nPrimer ID", "p2-1", "p2-2", "p2-3", "p5-1", "p5-2", "p5-3",
                       "p8-1",
                       "p8-2", "p8-3", "p10-2", "p10-3", "p12-1", "p12-2", "p12-3"]
        type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]
        type_order_ag = ["Synonymous", "Non-Synonymous"]
        context_order = ["UpA", "ApA", "CpA", "GpA"]
        context_order_uc = ["UpU", "UpA", "UpC", "UpG"]
        adar_preference = ["High", "Intermediate", "Low"]
        mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "G>C", "C>G", "A>U", "U>A", "G>U", "C>A"]
        Transitions_order = ["A>G", "U>C", "G>A", "C>U"]
        mutation_type_order = ["Transitions", "Transversions", "Oxidations"]
        plus_minus = u"\u00B1"
        x_order = range(0, 14, 1)
        x_ticks = ["RNA\nControl", "", "2", "", "", "5", "", "", "8", "", "10", "", "12", ""]
        dodge = 0.65
        # g1 = sns.catplot(x="passage", data=data_filter_replica,
        #                  hue="Mutation Type", order=x_order,
        #                  palette="Set2", kind="count", hue_order=["Transitions", "Transversions"], col="Mutation Type",
        #                  col_wrap=2,
        #                  col_order=["Transitions", "Transversions"], dodge=False, saturation=1)
        # g1.set_axis_labels("Passage", "Count")
        # g1.set(xticklabels=x_ticks)  # yscale='log',
        # plt.tight_layout()
        # g1.savefig(output_dir + "/All_Mutations_count_plot_replica{0}.png" .format(str(replica)), dpi=300)
        # plt.close()
        #
        # g1_mutation = sns.catplot(x="passage", data=data_filter_replica,
        #                           hue="Mutation", order=x_order, palette=mutation_palette(12), kind="count",
        #                           hue_order=mutation_order, dodge=False, saturation=1)
        # g1_mutation.set_axis_labels("Passage", "Count")
        # g1_mutation.set(xticklabels=x_ticks)  # yscale='log',
        # # plt.tight_layout()
        # g1_mutation.savefig(output_dir + "/All_Mutations_count_plot2_replica{0}.png".format(str(replica)), dpi=300)
        # plt.close()

        dist_g = sns.distplot(data_filter_replica["Frequency"])
        dist_g.set_xscale('log')
        dist_g.set_xlim(10 ** -5, 10 ** -2)
        plt.axvline(data_filter_replica["Frequency"].median(), color='red')
        plt.axvline(data_filter_replica["Frequency"].mean(), color='blue')

        sns.despine()
        plt.tight_layout()
        plt.savefig(output_dir + "/histplot_frequency_replica{0}".format(replica), dpi=300)
        plt.close()
        
        mutation_type_g = sns.catplot(x="passage", y="frac_and_weight", data=data_filter_replica, hue="Mutation Type",
                                      order=x_order, palette="Set2", kind="point", dodge=dodge,
                                      hue_order=mutation_type_order,
                                      estimator=weighted_varaint, orient="v", legend=True)
        mutation_type_g.set_axis_labels("Passage", "Variant Frequency {} CI=95%".format(plus_minus))
        mutation_type_g.set(yscale='log', ylim=(10 ** -5, 10 ** -2), xticklabels=x_ticks)
        plt.savefig(output_dir + "/All_Mutations_point_plot_replica{0}.png".format(str(replica)), dpi=300)
        plt.close()

        passage_g = sns.catplot(x="passage", y="frac_and_weight", data=data_filter_replica, hue="Mutation",
                                order=x_order, palette=mutation_palette(4), kind="point",
                                dodge=dodge, hue_order=Transitions_order,
                                join=False, estimator=weighted_varaint, orient="v", legend=True)
        passage_g.set_axis_labels("Passage", "Variant Frequency {} CI=95%".format(plus_minus))
        passage_g.set(yscale='log', ylim=(10 ** -5, 10 ** -2),
                      xticklabels=x_ticks)
        plt.savefig(output_dir + "/Transitions_Mutations_point_plot_RVB14_replica%s" % str(replica), dpi=300)
        plt.close()

    # data_filter["passage"] = data_filter["passage"].astype(int)
    #
    #
    # g4 = sns.relplot("passage", "frac_and_weight", data=data_filter, hue="Mutation", palette=mutation_palette(4),
    #                  hue_order=Transitions_order, estimator=weighted_varaint, col="Type", kind="line",
    #                  col_order=type_order)
    #
    # g4.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -5)
    # g4.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g4.savefig(output_dir + "/Time_Transitions_Mutations_line_plot", dpi=300)
    # plt.close()
        """ADAR preferences"""
        data_filter_replica_synonymous = data_filter_replica.loc[data_filter_replica.Type == "Synonymous"]
        #data_filter_replica_synonymous["ADAR_like"] = (data_filter_synonymous.Prev.str.contains('UpA') |data_filter_replica_synonymous.Prev.str.contains('ApA'))
        data_filter_replica_synonymous["Mutation_adar"] =data_filter_replica_synonymous["Mutation"]

        data_filter_replica_synonymous["Mutation_adar"] = np.where(((data_filter_replica_synonymous["Mutation"] == "A>G") &
                                                            (data_filter_replica_synonymous["5`_ADAR_Preference"] == "High")),
                                                           "High\nADAR-like\nA>G",
                                                           np.where(((data_filter_replica_synonymous["Mutation"] == "A>G")
                                                                     & (data_filter_replica_synonymous[
                                                                            "5`_ADAR_Preference"] == "Intermediate")),
                                                                    "Intermediate\nADAR-like\nA>G",
                                                                    np.where(
                                                                        ((data_filter_replica_synonymous["Mutation"] == "A>G") &
                                                                         (data_filter_replica_synonymous[
                                                                              "5`_ADAR_Preference"] == "Low")),
                                                                        "Low\nADAR-like\nA>G",
                                                                       data_filter_replica_synonymous["Mutation_adar"])))
        data_filter_replica_synonymous["Mutation_adar"] = np.where(((data_filter_replica_synonymous["Mutation"] == "U>C") &
                                                            (data_filter_replica_synonymous["3`_ADAR_Preference"] == "High")),
                                                           "High\nADAR-like\nU>C",
                                                           np.where(((data_filter_replica_synonymous["Mutation"] == "U>C")
                                                                     & (data_filter_replica_synonymous[
                                                                            "3`_ADAR_Preference"] == "Intermediate")),
                                                                    "Intermediate\nADAR-like\nU>C",
                                                                    np.where(
                                                                        ((data_filter_replica_synonymous["Mutation"] == "U>C") &
                                                                         (data_filter_replica_synonymous[
                                                                              "3`_ADAR_Preference"] == "Low")),
                                                                        "Low\nADAR-like\nU>C",
                                                                       data_filter_replica_synonymous["Mutation_adar"])))
        mutation_adar_order = ["High\nADAR-like\nA>G", "Low\nADAR-like\nA>G",
                               "High\nADAR-like\nU>C", "Low\nADAR-like\nU>C"]
        # data_filter_replica_synonymous["passage"] = data_filter_replica_synonymous["passage"].astype(str)
        # data_filter_replica_synonymous["passage"] = "p" + data_filter_replica_synonymous["passage"]
        catplot_adar = sns.catplot(x="passage", y="frac_and_weight", data=data_filter_replica_synonymous, hue="Mutation_adar",
                                   order=x_order, palette=mutation_palette(4, adar=True), kind="point", dodge=dodge,
                                   hue_order=mutation_adar_order, join=False, estimator=weighted_varaint, orient="v",
                                   legend=True) #comm
        catplot_adar.set_axis_labels("Passage", "Variant Frequency {} CI=95%".format(plus_minus))
        catplot_adar.set(yscale='log', ylim=(10 ** -6, 10 ** -2),
                         xticklabels=x_ticks)
        # catplot_adar.set_xticklabels(fontsize=8)
        # plt.tight_layout()
        plt.savefig(output_dir + "/adar_pref_mutation_point_plot_RVB14_replica{0}.png".format(replica), dpi=300)
        plt.close()


        """Stat"""
        df_sts = data_filter_replica.copy()
        df_sts["Frequency"] = np.where(df_sts["Prob"] < 0.95, 0, df_sts["Frequency"])
        # df_sts = df_sts.loc[df_sts["Prob"] > 0.95]
        df_sts["passage_p"] = np.where(df_sts["passage_p"] == "RNA\nControl", "RNA Control", df_sts["passage_p"])
        mutation_type_g1 = sns.boxplot(x="passage_p", y="Frequency", data=df_sts, hue="Mutation Type",
                                 order=passage_p_order, palette="Set2", dodge=True,
                                 hue_order=mutation_type_order)
        mutation_type_g1.set_yscale('log')
        mutation_type_g1.set_ylim(10 ** -5, 10 ** -2)
        mutation_type_g1.set(xlabel="Passage", ylabel="Variant Frequency")
        annot = Annotator(mutation_type_g1, trans_pairs, x="passage_p", y="Frequency", hue="Mutation Type",
                          data=df_sts, order=passage_p_order, hue_order=mutation_type_order)
        annot.configure(test='Kruskal', text_format='star', loc='outside', verbose=2,
                        comparisons_correction="Benjamini-Hochberg")#Bonferroni
        annot.apply_test()
        file_path = output_dir + "/sts_trans{0}.csv".format(replica)
        with open(file_path, "w") as o:
            with contextlib.redirect_stdout(o):
                mutation_type_g1, test_results = annot.annotate()
        plt.legend(bbox_to_anchor=(1.05, 0.5), loc=2, borderaxespad=0.)
        plt.tight_layout()
        plt.savefig(output_dir + "/All_Mutations_box_stat_plot_RVB14_replica{0}".format(replica), dpi=300)
        plt.close()

        dunn_df = posthoc_dunn(df_sts, val_col="Frequency", group_col="Mutation Type", p_adjust="fdr_bh")
        dunn_df.to_csv(output_dir + "/dunn_df_All_Mutations_replica{0}.csv".format(str(replica)), sep=",")

        passage_g1 = sns.boxplot(x="passage_p", y="Frequency", data=df_sts, hue="Mutation",
                                 order= passage_p_order, palette=mutation_palette(4), dodge=True,
                                 hue_order=Transitions_order)
        passage_g1.set_yscale('log')
        passage_g1.set_ylim(10 ** -5, 10 ** -2)
        passage_g1.set(xlabel="Passage", ylabel="Variant Frequency")
        annot = Annotator(passage_g1, pairs, x="passage_p", y="Frequency", hue="Mutation", data=df_sts,
                          order= passage_p_order, hue_order=Transitions_order)
        annot.configure(test='Kruskal', text_format='star', loc='outside', verbose=2,
                        comparisons_correction="Benjamini-Hochberg")#Bonferroni
        annot.apply_test()
        file_path = output_dir + "/sts{0}.csv".format(replica)
        with open(file_path, "w") as o:
            with contextlib.redirect_stdout(o):
                passage_g1, test_results = annot.annotate()

        plt.legend(bbox_to_anchor=(1.05, 0.5), loc=2, borderaxespad=0.)
        plt.tight_layout()
        plt.savefig(output_dir + "/Transitions_Mutations_box_stat_plot_RVB14_replica{0}".format(replica), dpi=300)
        plt.close()
        dunn_df = posthoc_dunn(data_filter_replica, val_col="Frequency", group_col="Mutation", p_adjust="fdr_bh")
        dunn_df.to_csv(output_dir + "/dunn_df_Transitions_Mutations_replica{0}.csv".format(str(replica)), sep=",")

        df_sts_syn = data_filter_replica_synonymous.copy()
        # df_sts_syn["Frequency"] = np.where(df_sts_syn["Prob"] < 0.95, 0, df_sts_syn["Frequency"])
        df_sts_syn = df_sts_syn.loc[df_sts_syn["Prob"] > 0.95]
        df_sts_syn["Mutation_adar"] = np.where(df_sts_syn["Mutation_adar"] == "High\nADAR-like\nA>G",
                                               "High ADAR-like A>G",
                                               np.where(df_sts_syn["Mutation_adar"] == "Intermediate\nADAR-like\nA>G",
                                                        "Intermediate ADAR-like A>G",
                                                        np.where(df_sts_syn["Mutation_adar"] == "Low\nADAR-like\nA>G",
                                                                 "Low ADAR-like A>G",
                                                                 np.where(df_sts_syn[
                                                                              "Mutation_adar"] == "High\nADAR-like\nU>C",
                                                                          "High ADAR-like U>C",
                                                                          np.where(df_sts_syn[
                                                                                       "Mutation_adar"] == "Intermediate\nADAR-like\nU>C",
                                                                                   "Intermediate ADAR-like U>C",
                                                                                   np.where(df_sts_syn[
                                                                                                "Mutation_adar"] == "Low\nADAR-like\nU>C",
                                                                                            "Low ADAR-like U>C",
                                                                                            df_sts_syn[
                                                                                                "Mutation_adar"]))))))
        mutation_adar_order = ["High ADAR-like A>G", "Low ADAR-like A>G",
                               "High ADAR-like U>C", "Low ADAR-like U>C"]
        adar_g = sns.boxplot(x="passage_p", y="Frequency", data=df_sts_syn, hue="Mutation_adar",
                             order=passage_p_order, palette=mutation_palette(4, adar=True), dodge=True,
                             hue_order=mutation_adar_order)
        adar_g.set_yscale('log')
        adar_g.set_ylim(10 ** -6, 10 ** -1)
        adar_g.set(xlabel="Passage", ylabel="Variant Frequency")
        annot = Annotator(adar_g, pairs_adar, x="passage_p", y="Frequency", hue="Mutation_adar", data=df_sts_syn,
                          hue_order=mutation_adar_order, order=passage_p_order)
        annot.configure(test='Kruskal', text_format='star', loc='outside', verbose=2,
                        comparisons_correction="Benjamini-Hochberg")#Bonferroni
        annot.apply_test()
        file_path = output_dir + "/sts_adar_{0}.csv".format(replica)
        with open(file_path, "w") as o:
            with contextlib.redirect_stdout(o):
                adar_g, test_results = annot.annotate()
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.tight_layout()
        plt.savefig(output_dir + "/adar_pref_mutation_box_stat_plot_RVB14_replica{0}".format(replica), dpi=300)
        plt.close()
        dunn_df = posthoc_dunn(df_sts_syn, val_col="Frequency", group_col="Mutation_adar", p_adjust="fdr_bh")
        dunn_df.to_csv(output_dir + "/dunn_df_adar_pref_replica{0}.csv".format(str(replica)), sep=",")

    # data_filter_pass5 = data_filter.loc[data_filter.passage == 5]
    # ata_filter_pass5 = data_filter.loc[data_filter.replica == 2]
    # data_filter_pass5 = data_filter_pass5[data_filter_pass5["pval"] < 0.01]
    # data_filter_pass5 = data_filter_pass5.loc[data_filter_pass5.Type == "Synonymous"]
    # # data_filter_pass5["ADAR_like"] = (data_filter_pass5.Prev.str.contains('UpA') | data_filter_pass5.Prev.str.contains('ApA'))
    # data_filter_pass5["Mutation"] = np.where(((data_filter_pass5["Mutation"] == "A>G") &
    #                                          (data_filter_pass5["5`_ADAR_Preference"] == "High")),
    #                                          "High\nADAR-like\nA>G", np.where(((data_filter_pass5["Mutation"] == "A>G")
    #                                                   & (data_filter_pass5["5`_ADAR_Preference"] == "Intermediate")),
    #                                                   "Intermediate\nADAR-like\nA>G",
    #                                                   np.where(((data_filter_pass5["Mutation"] == "A>G") &
    #                                                             (data_filter_pass5["5`_ADAR_Preference"] == "Low")),
    #                                                            "Low\nADAR-like\nA>G",
    #                                                            data_filter_pass5["Mutation"])))
    # data_filter_pass5["Mutation_adar"] = np.where(((data_filter_pass5["Mutation"] == "U>C") &
    #                                          (data_filter_pass5["3`_ADAR_Preference"] == "High")),
    #                                          "High\nADAR-like\nU>C", np.where(((data_filter_pass5["Mutation"] == "U>C")
    #                                                   & (data_filter_pass5["3`_ADAR_Preference"] == "Intermediate")),
    #                                                   "Intermediate\nADAR-like\nU>C",
    #                                                   np.where(((data_filter_pass5["Mutation"] == "U>C") &
    #                                                             (data_filter_pass5["3`_ADAR_Preference"] == "Low")),
    #                                                            "Low\nADAR-like\nU>C",
    #                                                            data_filter_pass5["Mutation"])))
    # mutation_adar_order = ["High\nADAR-like\nA>G", "Intermediate\nADAR-like\nA>G", "Low\nADAR-like\nA>G",
    #                        "High\nADAR-like\nU>C", "Intermediate\nADAR-like\nU>C", "Low\nADAR-like\nU>C", "G>A",
    #                        "C>U"]
    #
    # data_filter_pass5["log10_Frequency"] = data_filter_pass5["Frequency"].apply(lambda x: np.log10(x))
    # print(data_filter_pass5.to_string())
    # ax = sns.violinplot("Mutation_adar", "log10_Frequency", data=data_filter_pass5, palette=mutation_palette(8, gray=True),
    #                     order=mutation_adar_order)
    # # # ax = sns.stripplot("Mutation_adar", "Frequency", data=data_filter_pass5, color=".2", order=mutation_adar_order,
    # # #                    alpha=.25)
    # old_statannot.add_stat_annotation(ax, data=data_filter_pass5, x="Mutation_adar", y="log10_Frequency",
    #                     boxPairList=[("High\nADAR-like\nA>G", "High\nADAR-like\nU>C"),("High\nADAR-like\nA>G", "G>A"),
    #                                  ("High\nADAR-like\nA>G", "C>U")], test='Mann-Whitney', textFormat='star',
    #                                   loc='outside', verbose=2,
    #                                   order=mutation_adar_order)
    # #
    # # ax.set_yscale('log')
    # ax.set_xlabel("Mutation")
    # ax.set_ylabel("Variant Frequency [log10]")
    # ax.set(ylim=(-5, -1))
    # plt.xticks(fontsize=7)
    # sns.despine()
    # plt.tight_layout()
    # # # plt.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/context_p5_point_plot_v2.png", dpi=300)
    # plt.savefig(output_dir + "/mutation_p5_box_plot_RVB14.png", dpi=300)
    # plt.close()

    # # A>G Prev Context
    # data_filter_ag["ADAR_like"] = (data_filter_ag.Prev.str.contains('UpA') | data_filter_ag.Prev.str.contains('ApA'))
    #
    # g5 = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=label_order,
    #                  palette=mutation_palette(2), kind="point", dodge=True, hue_order=[True, False], estimator=weighted_varaint,
    #                  orient="v", col="Type", join=False, col_order=type_order_ag)
    # g5.set_axis_labels("", "Variant Frequency")
    # g5.set(yscale='log')
    # g5.set(ylim=(7 * 10 ** -7, 4 * 10 ** -3))
    # g5.set_xticklabels(fontsize=9, rotation=90)
    # # plt.show()
    # g5.savefig(output_dir + "/Context_point_plot", dpi=300)
    # plt.close()
    #
    # adar_all_g = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=label_order,
    #                          palette=mutation_palette(2), kind="point", dodge=True, hue_order=[True, False],
    #                          estimator=weighted_varaint, orient="v", join=False)
    # adar_all_g.set_axis_labels("", "Variant Frequency")
    # adar_all_g.set(yscale='log')
    # adar_all_g.set(ylim=(7 * 10 ** -7, 4 * 10 ** -3))
    # adar_all_g.set_xticklabels(rotation=45)
    # # plt.show()
    # adar_all_g.savefig(output_dir + "/Context_point_all_mutations_plot", dpi=300)
    # plt.close()
    # data_filter_ag["passage"] = data_filter_ag["passage"].astype(int)
    #
    # data_filter_ag_pass5 = data_filter_ag.loc[data_filter_ag.passage == 5]
    # data_filter_ag_pass5 = data_filter_ag_pass5[data_filter_ag_pass5["pval"] < 0.01]
    # data_filter_ag_pass5 = data_filter_ag_pass5.loc[data_filter_ag_pass5.Type == "Synonymous"]
    #
    # g_context = sns.catplot("ADAR_like", "frac_and_weight", data=data_filter_ag_pass5,
    #                         order=[True, False], palette=mutation_palette(2), kind="point",
    #                         join=False, estimator=weighted_varaint, orient="v", dodge=True)
    # # add_stat_annotation(g_context, data=data_filter_ag_pass5, x="ADAR_like", y="Frequency", order=[True, False],
    # #                     boxPairList=[(True, False)], test='Mann-Whitney', textFormat='star', loc='inside', verbose=2)
    # g_context.set_axis_labels("ADAR-like\nContext", "Variant Frequency")
    # g_context.set_xticklabels(rotation=45)
    # g_context.set(yscale='log')
    # g_context.set(ylim=(5 * 10 ** -5, 5 * 10 ** -3))
    # plt.tight_layout()
    # # plt.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/context_p5_point_plot.png", dpi=300)
    # plt.savefig(output_dir + "/context_p5_point_plot.png", dpi=300)
    # plt.close()
    #
    # ax = sns.boxplot("ADAR_like", "Frequency", data=data_filter_ag_pass5, palette=mutation_palette(2), order=(True, False))
    # ax = sns.stripplot("ADAR_like", "Frequency", data=data_filter_ag_pass5, color=".2", order=(True, False))
    # old_statannot.add_stat_annotation(ax, data=data_filter_ag_pass5, x="ADAR_like", y="Frequency",
    #                     boxPairList=[(True, False)], test='Mann-Whitney', textFormat='star', loc='outside', verbose=2)
    # ax.set_yscale('log')
    # ax.set_xlabel("ADAR-like\nContext")
    # ax.set_ylabel("Variant Frequency")
    # ax.set(ylim=(10 ** -5, 10 ** -2))
    # sns.despine()
    # plt.tight_layout()
    # # plt.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/context_p5_point_plot_v2.png", dpi=300)
    # plt.savefig(output_dir + "/context_p5_box_plot_RVB14.png", dpi=300)
    # plt.close()

    # data_codons = data_filter_ag[data_filter_ag["pval"] < 0.01]
    # data_codons = data_codons[data_codons["ADAR_like"] == True]
    # data_codons = data_codons[data_codons["Type"] == "Synonymous"]
    # g_codons = sns.catplot("label", "frac_and_weight", data=data_codons, hue="Consensus>Mutated_codon",
    #                        order=label_order, palette="tab20", kind="point",
    #                        join=False, estimator=weighted_varaint, orient="v", dodge=True)
    # g_codons.set_axis_labels("", "Variant Frequency")
    # g_codons.set_xticklabels(fontsize=5, rotation=45)
    # g_codons.set(yscale='log')
    # # g_codons.set(ylim=(10**-7, 10**-3))
    #
    # # plt.show()
    # g_codons.savefig(output_dir + "/codons_point_plot", dpi=300)
    # plt.close()
    #
    # g_codons_3c = sns.catplot("label", "frac_and_weight", data=data_codons[data_codons["Protein"] == "3C"],
    #                           hue="Consensus>Mutated_codon",
    #                           order=label_order, palette="tab20", kind="point",
    #                           join=False, estimator=weighted_varaint, orient="v", dodge=True)
    # g_codons_3c.set_axis_labels("", "Variant Frequency")
    # g_codons_3c.set_xticklabels(fontsize=5, rotation=45)
    # g_codons_3c.set(yscale='log')
    # # g_codons.set(ylim=(10**-7, 10**-3))
    #
    # # plt.show()
    # g_codons_3c.savefig(output_dir + "/codons_point_plot_3c", dpi=300)
    # plt.close()
    #
    # g6 = sns.relplot("passage", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", palette="Set2",
    #                  hue_order=[True, False], estimator=weighted_varaint, col="Type", kind="line",
    #                  col_order=type_order)
    #
    # g6.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -5)
    # g6.set(ylim=(0, 10 ** -2))
    # yaxis = plt.gca().yaxis
    # yaxis.set_minor_locator(fits_new_plotter.MinorSymLogLocator(1e-1))
    # g6.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g6.savefig(output_dir + "/Context_miseq_sample_time_line_plot", dpi=300)
    # plt.close()
    #
    # g7 = sns.relplot(x="passage", y="frac_and_weight", hue="Prev", data=data_filter_ag, palette="Paired",
    #                  kind="line",style="Type", style_order=type_order, hue_order=context_order,
    #                  estimator=weighted_varaint)
    # g7.set(yscale="log")
    # g7.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    # g7.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g7.savefig(output_dir + "/ADAR_like_AG_Mutation_Context_trajectories_Context_.png", dpi=300)
    # plt.close()
    #
    # g8 = sns.relplot(x="passage", y="frac_and_weight", hue="ADAR_like", data=data_filter_ag, palette="Paired",
    #                  kind="line", style="Type", style_order=type_order,
    #                  estimator=weighted_varaint)  # , hue_order=context_order)
    # g8.set(yscale="log")
    # g8.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    # g8.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g8.savefig(output_dir + "/ADAR_like_AG_Mutation_ADAR_trajectories_Context.png", dpi=300)
    # plt.close()
    #
    #
    # mutation_g9 = sns.relplot(x="passage", y="frac_and_weight", hue="ADAR_like", data=data_filter_ag, palette=flatui,
    #                           kind="line", estimator=weighted_varaint, hue_order=[True, False])
    # mutation_g9.set(yscale="log")
    # mutation_g9.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    # mutation_g9.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # mutation_g9.savefig(output_dir + "/ADAR_like_AG_Mutation_ADAR_trajectories_Context_all_mutation.png", dpi=300)
    # plt.close()
    # """U>C Context"""
    # g9 = sns.catplot("label", "frac_and_weight", data=data_filter_uc, hue="Next", order=label_order, palette="Set2",
    #                  hue_order=context_order_uc, estimator=weighted_varaint, orient="v", dodge=True, kind="point",
    #                  col="Type", join=False, col_order=type_order_ag)
    # g9.set_axis_labels("", "Variant Frequency")
    # g9.set(yscale='log')
    # g9.set(ylim=(10 ** -7, 10 ** -3))
    # g9.set_xticklabels(rotation=45)
    # # plt.show()
    # g9.savefig(output_dir + "/UC_Context_point_plot", dpi=300)
    # plt.close()
    #
    # g10 = sns.relplot("passage", "frac_and_weight", data=data_filter_uc, hue="Next", palette="Set2",
    #                   hue_order=context_order_uc, estimator=weighted_varaint, col="Type", kind="line",
    #                   col_order=type_order_ag)
    #
    # # g8.set(yscale="log")
    # g10.axes.flat[0].set_yscale('symlog', linthreshy=10 ** -4)
    # # g8.set(ylim=(0, 10 ** -6))
    # # yaxis = plt.gca().yaxis
    # yaxis.set_minor_locator(fits_new_plotter.MinorSymLogLocator(1e-1))
    # g10.set_axis_labels("Passage", "Variant Frequency")
    # # plt.show()
    # g10.savefig(output_dir + "/UC_Context_miseq_sample_time_line_plot", dpi=300)
    # plt.close()


if __name__ == "__main__":
    main()

    # data_filter_ag_grouped = data_filter_ag.groupby(["label", "Type", "passage", "replica",
    #                                                  "ADAR_grade_five", "5`_ADAR_Preference"])["frac_and_weight"].agg(
    #     lambda x: weighted_varaint(x))
    # data_filter_ag_grouped = data_filter_ag_grouped.reset_index()
    #
    # # print(data_filter_ag_grouped.to_string())
    #
    # data_filter_ag_grouped = data_filter_ag_grouped.rename(columns={"frac_and_weight": "Frequency"})
    # data_filter_ag_grouped["Frequency"] = data_filter_ag_grouped["Frequency"].astype(float)
    # data_filter_ag_grouped = data_filter_ag_grouped[data_filter_ag_grouped["label"] != "RNA Control\nPrimer ID"]
    # data_filter_ag_grouped = data_filter_ag_grouped[data_filter_ag_grouped["label"] != "RNA Control\nRND"]
    # data_filter_ag_grouped = data_filter_ag_grouped[data_filter_ag_grouped["replica"] == replica]
    #
    # data_reg_full_adar = data_filter_ag_grouped[data_filter_ag_grouped["ADAR_grade_five"] == 1]
    # data_reg_semi_adar = data_filter_ag_grouped[data_filter_ag_grouped["ADAR_grade_five"] == 0.5]
    # data_reg_nonadar = data_filter_ag_grouped[data_filter_ag_grouped["ADAR_grade_five"] == 0]
    #
    # data_reg_full_adar_syn = data_reg_full_adar[data_reg_full_adar["Type"] == "Synonymous"]
    # data_reg_semi_adar_syn = data_reg_semi_adar[data_reg_semi_adar["Type"] == "Synonymous"]
    # data_reg_nonadar_syn = data_reg_nonadar[data_reg_nonadar["Type"] == "Synonymous"]
    #
    # data_reg_full_adar_non_syn = data_reg_full_adar[data_reg_full_adar["Type"] == "Non-Synonymous"]
    # data_reg_semi_adar_non_syn = data_reg_semi_adar[data_reg_semi_adar["Type"] == "Non-Synonymous"]
    # data_reg_nonadar_non_syn = data_reg_nonadar[data_reg_nonadar["Type"] == "Non-Synonymous"]
    # fig, axes = plt.subplots(2, 3, sharey="all", sharex="all")
    #
    # from sklearn import linear_model
    # regr = linear_model.LinearRegression()
    # X = data_reg_full_adar_syn.passage.values.reshape(-1, 1)
    # y = data_reg_full_adar_syn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope1 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept1 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("100-ADAR SYN: y=%sx + %s" % (str(regr.coef_[0]).split("[")[1].split("]")[0], str(regr.intercept_).split("[")[1].split("]")[0]))
    #
    # g1 = sns.regplot(x="passage", y="Frequency", data=data_reg_full_adar_syn, ax=axes[0, 0],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope1, intercept1)})
    # g1.set(title="100-ADAR-like Synonymous")
    #
    # X = data_reg_semi_adar_syn.passage.values.reshape(-1, 1)
    # y = data_reg_semi_adar_syn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope2 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept2 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("50-ADAR SYN: y=%sx + %s" % (str(regr.coef_[0]).split("[")[1].split("]")[0], str(regr.intercept_).split("[")[1].split("]")[0]))
    #
    # g2 = sns.regplot(x="passage", y="Frequency", data=data_reg_semi_adar_syn, ax=axes[0, 1],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope2, intercept2)})
    # g2.set(title="50-ADAR-like Synonymous")
    #
    # X = data_reg_nonadar_syn.passage.values.reshape(-1, 1)
    # y = data_reg_nonadar_syn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope3 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept3 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("NON-ADAR SYN: y=%sx + %s" % (str(regr.coef_[0]).split("[")[1].split("]")[0], str(regr.intercept_).split("[")[1].split("]")[0]))
    #
    # g3 = sns.regplot(x="passage", y="Frequency", data=data_reg_nonadar_syn, ax=axes[0, 2],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope3, intercept3)})
    # g3.set(title="Non ADAR-like Synonymous")
    #
    # X = data_reg_full_adar_non_syn.passage.values.reshape(-1, 1)
    # y = data_reg_full_adar_non_syn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope4 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept4 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("100-ADAR NON-SYN: y=%sx + %s" % (str(regr.coef_[0]).split("[")[1].split("]")[0], str(regr.intercept_).split("[")[1].split("]")[0]))
    #
    #
    # g4 = sns.regplot(x="passage", y="Frequency", data=data_reg_full_adar_non_syn, ax=axes[1, 0],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope4, intercept4)})
    # g4.set(title="100-ADAR-like Non Synonymous")
    #
    # X = data_reg_semi_adar_non_syn.passage.values.reshape(-1, 1)
    # y = data_reg_semi_adar_non_syn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope5 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept5 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("50-ADAR NON-SYN: y=%sx + %s" % (str(regr.coef_[0]).split("[")[1].split("]")[0], str(regr.intercept_).split("[")[1].split("]")[0]))
    #
    # g5 = sns.regplot(x="passage", y="Frequency", data=data_reg_semi_adar_non_syn, ax=axes[1, 1],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope5, intercept5)})
    # g5.set(title="50-ADAR-like Non Synonymous")
    #
    # X = data_reg_nonadar_non_syn.passage.values.reshape(-1, 1)
    # y = data_reg_nonadar_non_syn.Frequency.values.reshape(-1, 1)
    # regr.fit(X, y)
    # slope6 = float(str(regr.coef_[0]).split("[")[1].split("]")[0])
    # intercept6 = float(str(regr.intercept_).split("[")[1].split("]")[0])
    # print("NON-ADAR NON-SYN: y=%sx + %s" % (str(regr.coef_[0]).split("[")[1].split("]")[0], str(regr.intercept_).split("[")[1].split("]")[0]))
    #
    # g6 = sns.regplot(x="passage", y="Frequency", data=data_reg_nonadar_non_syn, ax=axes[1, 2],
    #                  line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope6, intercept6)})
    # g6.set(title="Non ADAR-like Non Synonymous")
    #
    # axes[0, 0].legend()
    # axes[0, 0].legend(loc=2)
    # axes[0, 0].set_yscale("log")
    # axes[0, 0].set_xlabel('')
    # axes[0, 0].set_ylabel('')
    # axes[0, 0].set_xlim(0, 14, 2)
    # axes[0, 0].set_ylim(10 ** -6, 10 ** -2)
    # axes[0, 1].legend()
    # axes[0, 1].legend(loc=2)
    # axes[0, 1].set_yscale("log")
    # axes[0, 1].set_xlabel('')
    # axes[0, 1].set_ylabel('')
    # axes[0, 1].set_xlim(0, 14, 2)
    # axes[0, 1].set_ylim(10 ** -6, 10 ** -2)
    # axes[0, 2].legend()
    # axes[0, 2].legend(loc=2)
    # axes[0, 2].set_yscale("log")
    # axes[0, 2].set_xlabel('')
    # axes[0, 2].set_ylabel('')
    # axes[0, 2].set_xlim(0, 14, 2)
    # axes[0, 2].set_ylim(10 ** -6, 10 ** -2)
    # axes[1, 0].legend()
    # axes[1, 0].legend(loc=2)
    # axes[1, 0].set_yscale("log")
    # axes[1, 0].set_xlabel('')
    # axes[1, 0].set_ylabel('')
    # axes[1, 0].set_xlim(0, 14, 2)
    # axes[1, 0].set_ylim(10 ** -6, 10 ** -2)
    # axes[1, 1].legend()
    # axes[1, 1].legend(loc=2)
    # axes[1, 1].set_yscale("log")
    # axes[1, 1].set_xlabel('')
    # axes[1, 1].set_ylabel('')
    # axes[1, 1].set_xlim(0, 14, 2)
    # axes[1, 1].set_ylim(10 ** -6, 10 ** -2)
    # axes[1, 2].legend()
    # axes[1, 2].legend(loc=2)
    # axes[1, 2].set_yscale("log")
    # axes[1, 2].set_xlabel('')
    # axes[1, 2].set_ylabel('')
    # axes[1, 2].set_xlim(0, 14, 2)
    # axes[1, 2].set_ylim(10 ** -6, 10 ** -2)
    # fig.text(0.5, 0.01, 'Passage', ha='center')
    # fig.text(0.001, 0.5, 'Frequency', va='center', rotation='vertical')
    # plt.savefig(output_dir + "/regplot_AG_Mutation_Context_trajectories.png", dpi=300)
    # plt.close()

    # g2 = sns.catplot(x="label", y="frac_and_weight", data=data_filter, hue="Mutation", order=label_order,
    #                  palette=mutation_palette(4), kind="point", dodge=True, hue_order=Transitions_order, join=False,
    #                  estimator=weighted_varaint,
    #                  orient="v", legend=True)
    # g2.set_axis_labels("", "Variant Frequency")
    # g2.set(yscale='log', ylim=(10 ** -6, 10 ** -2), xlim=(0, 12, 2))
    # # g2.set_yticklabels(fontsize=12)
    # g2.set_xticklabels(fontsize=9, rotation=90)
    # plt.show()
    # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/Transitions_Mutations_point_plot_RV", dpi=300)
    # g2.savefig(output_dir + "/Transitions_Mutations_point_plot", dpi=300)
    # plt.close()