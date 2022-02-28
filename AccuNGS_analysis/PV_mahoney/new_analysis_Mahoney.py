
from datetime import datetime

import pandas as pd

from AccuNGS_analysis.new_analysis_fuctions import *


def main():
    # input_dir = "/Users/odedkushnir/PhD_Projects/After_review/CirSeq/PV/Mahoney/"
    input_dir = "C:/Users/odedku/PhD_Projects/After_review/Cirseq/PV/Mahoney/"
    date = datetime.today().strftime("%Y%m%d")
    prefix = "inosine_predict_context"
    output_dir = input_dir + prefix
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    q_file_name = "q30_data_mutation.csv"
    # data_adar = pd.read_csv("/Users/odedkushnir/PhD_Projects/fitness/CirSeq/PV/Mahoney/inosine_results/Output/PV1_adar1_trans.csv")
    data_adar = pd.read_csv("C:/Users/odedku/PhD_Projects/After_review/Inosine_Predict/Output/Mahoney_adar1_trans.csv")
    columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
               "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts",
               "Consensus>Mutated_codon", "fiveGrade", "threeGrade"]
    removed_mutation = "C>U"
    data_filter = analysis(input_dir, output_dir, q_file_name, data_adar, columns, removed_mutation)

    """Plots"""
    passage_order = ["p3", "p4", "p5", "p6", "p7", "p8"]
    transition_order = ["A>G", "U>C", "G>A"]
    pairs = [(("p3", "A>G"), ("p3", "G>A")), (("p4", "A>G"), ("p4", "G>A")),
             (("p5", "A>G"), ("p5", "G>A")), (("p6", "A>G"), ("p6", "G>A")),
             (("p7", "A>G"), ("p7", "G>A")), (("p8", "A>G"), ("p8", "G>A")),
             (("p3", "A>G"), ("p3", "U>C")), (("p4", "A>G"), ("p4", "U>C")),
             (("p5", "A>G"), ("p5", "U>C")), (("p6", "A>G"), ("p6", "U>C")),
             (("p7", "A>G"), ("p7", "U>C")), (("p8", "A>G"), ("p8", "U>C"))]

    pairs_adar = [(("p3", "High\nADAR-like\nA>G"), ("p3", "Low\nADAR-like\nA>G")),
                  (("p4", "High\nADAR-like\nA>G"), ("p4", "Low\nADAR-like\nA>G")),
                  (("p5", "High\nADAR-like\nA>G"), ("p5", "Low\nADAR-like\nA>G")),
                  (("p6", "High\nADAR-like\nA>G"), ("p6", "Low\nADAR-like\nA>G")),
                  (("p7", "High\nADAR-like\nA>G"), ("p7", "Low\nADAR-like\nA>G")),
                  (("p8", "High\nADAR-like\nA>G"), ("p8", "Low\nADAR-like\nA>G")),
                  (("p3", "High\nADAR-like\nU>C"), ("p3", "Low\nADAR-like\nU>C")),
                  (("p4", "High\nADAR-like\nU>C"), ("p4", "Low\nADAR-like\nU>C")),
                  (("p5", "High\nADAR-like\nU>C"), ("p5", "Low\nADAR-like\nU>C")),
                  (("p6", "High\nADAR-like\nU>C"), ("p6", "Low\nADAR-like\nU>C")),
                  (("p7", "High\nADAR-like\nU>C"), ("p7", "Low\nADAR-like\nU>C")),
                  (("p8", "High\nADAR-like\nU>C"), ("p8", "Low\nADAR-like\nU>C"))]
    label_order = ["PV-p3", "PV-p4", "PV-p5", "PV-p6", "PV-p7", "PV-p8"]
    x_order = range(0, 9, 1)
    plots(input_dir, date, data_filter, "PV1", passage_order, transition_order, pairs, label_order, pairs_adar,
          x_order, x_ticks=["", "", "", "3", "4", "5", "6", "7", "8"], filter_reads=True, ylim=(10**-5, 10**-2))


if __name__ == "__main__":
    main()

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



