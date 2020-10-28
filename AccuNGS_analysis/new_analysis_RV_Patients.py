
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

sns.set(font_scale=1.2)
sns.set_style("ticks")
sns.despine()

# print(plt.style.available)

def weighted_varaint(x, **kws):
    var, count = map(np.asarray, zip(*x))
    return var.sum() / count.sum()


def main():
    # input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/"
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/patients"
    output_dir = input_dir + "/20201027_plots"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)


    data_mutations = pd.read_csv(input_dir + "/q30_data_mutation.csv")

    columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
               "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts",
              "Consensus>Mutated_codon", "Protein"]
    data_filter = pd.DataFrame(data_mutations, columns=columns)
    data_filter["pval"] = data_filter["pval"].fillna(1)
    data_filter["no_variants"] = data_filter["Frequency"] * data_filter["Read_count"]

    # filter based on pval<0.01 and Prob>0.95
    data_filter["no_variants"] = np.where(data_filter["pval"] > 0.01, 0, data_filter["no_variants"])
    data_filter["no_variants"] = np.where(data_filter["Prob"] < 0.95, 0, data_filter["no_variants"])



    data_filter["frac_and_weight"] = list(zip(data_filter.no_variants, data_filter.Read_count))

    data_filter["Type"] = data_filter["Type"].fillna("NonCodingRegion")
    data_filter.to_csv(output_dir + "/data_filter.csv", sep=',', encoding='utf-8')
    data_filter_misense = data_filter[data_filter["Type"] == "Non-Synonymous"]

    label_order = ["RNA Control\nPrimer ID", "Cell Cultureֿ\nControl", "Patient-1", "Patient-4", "Patient-5",
                   "Patient-9", "Patient-16", "Patient-17", "Patient-20"]
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U"]
    transition_order = ["A>G", "U>C", "G>A", "C>U"]
    type_order1 = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]

    # A>G Prev Context
    data_filter_ag = data_filter[data_filter["Mutation"] == "A>G"]
    data_filter_ag = data_filter_ag.rename(columns={"Prev": "Context"})

    data_filter_ag['Context'].replace('AA', 'ApA', inplace=True)
    data_filter_ag['Context'].replace('UA', 'UpA', inplace=True)
    data_filter_ag['Context'].replace('CA', 'CpA', inplace=True)
    data_filter_ag['Context'].replace('GA', 'GpA', inplace=True)

    context_order = ["UpA", "ApA", "CpA", "GpA"]
    type_order2 = ["Synonymous", "Non-Synonymous"]

    data_filter_ag["ADAR_like"] = data_filter_ag.Context.str.contains('UpA') | data_filter_ag.Context.str.contains('ApA')
    print(data_filter_ag.to_string())
    data_filter_ag.to_csv(output_dir + "/data_mutation_AG_trajectories.csv", sep=',', encoding='utf-8')

    # U>C Prev Context
    data_filter_uc = data_filter[data_filter["Mutation"] == "U>C"]

    data_filter_uc['Next'].replace('UA', 'UpA', inplace=True)
    data_filter_uc['Next'].replace('UU', 'UpU', inplace=True)
    data_filter_uc['Next'].replace('UC', 'UpC', inplace=True)
    data_filter_uc['Next'].replace('UG', 'UpG', inplace=True)

    data_filter_uc.to_csv(output_dir + "/data_filter_uc.csv", sep=',', encoding='utf-8')
    context_order_uc = ["UpA", "UpU", "UpG",  "UpC"]

    #Plots
    g1 = sns.catplot(x="label", y="frac_and_weight", data=data_filter, hue="Mutation", order=label_order, palette="tab20",
                        kind="point", dodge=False, hue_order=mutation_order, join=True, estimator=weighted_varaint,
                     orient="v")
    g1.set_axis_labels("", "Variant Frequency")
    g1.set_xticklabels(fontsize=9, rotation=45)
    g1.set(yscale='log')
    # g1.set(ylim=(10**-7, 10**-3))

    # plt.show()
    g1.savefig(output_dir + "/All_Mutations_point_plot", dpi=300)
    plt.close()
    g2 = sns.catplot(x="label", y="frac_and_weight", data=data_filter, hue="Mutation", order=label_order, palette="tab20"
                        ,kind="point", dodge=True, hue_order=transition_order, join=False, estimator=weighted_varaint,
                     orient="v", legend=True)
    g2.set_axis_labels("", "Variant Frequency")
    g2.set(yscale='log')
    g2.set(ylim=(10 ** -6, 10 ** -2))
    # g2.set_yticklabels(fontsize=12)
    g2.set_xticklabels(fontsize=10, rotation=90)
    # plt.show()
    # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/Transition_Mutations_point_plot_RV", dpi=300)
    g2.savefig(output_dir + "/Transition_Mutations_point_plot", dpi=300)
    g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/Prgress reports/20200913 Final report/plots" +
                      "/Fig9a_Transition_Mutations_point_plot_Patients", dpi=300)
    plt.close()

    # g_rna = sns.catplot(x="RNA", y="frac_and_weight", data=data_filter, hue="Mutation", order=rna_order,
    #                  palette="tab20", kind="point", dodge=True, hue_order=transition_order, join=False, estimator=weighted_varaint,
    #                  orient="v", legend=True)
    # g_rna.set_axis_labels("", "Variant Frequency")
    # g_rna.set(yscale='log')
    # g_rna.set(ylim=(10 ** -6, 10 ** -2))
    # # g2.set_yticklabels(fontsize=12)
    # g_rna.set_xticklabels(fontsize=10, rotation=45)
    # plt.show()
    # g2.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/Transition_Mutations_point_plot_RV", dpi=300)
    # g_rna.savefig(output_dir + "/Transition_Mutations_point_RNA_plot", dpi=300)
    # plt.close()

    # A>G Prev Context
    flatui = ["#3498db", "#9b59b6"]
    g5 = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=label_order, palette=flatui,
                kind="point", dodge=True, hue_order=[True, False], estimator=weighted_varaint, orient="v",
                     col="Type", join=False, col_order=type_order2)
    g5.set_axis_labels("", "Variant Frequency")
    g5.set(yscale='log')
    g5.set(ylim=(7*10**-7, 4*10**-3))
    g5.set_xticklabels(rotation=90)
    # plt.show()
    g5.savefig(output_dir + "/Context_point_plot", dpi=300)
    g5.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/Prgress reports/20200913 Final report/plots" +
               "/Fig9b_Context_point_plot_Patients", dpi=300)
    plt.close()

    g6 = sns.catplot("label", "frac_and_weight", data=data_filter_ag, hue="ADAR_like", order=label_order, palette=flatui,
                kind="point", dodge=True, hue_order=[True, False], estimator=weighted_varaint, orient="v",
                     join=False)
    g6.set_axis_labels("", "Variant Frequency")
    g6.set(yscale='log')
    g6.set(ylim=(7*10**-7, 4*10**-3))
    g6.set_xticklabels(rotation=45)
    # plt.show()
    g6.savefig(output_dir + "/Context_point_all_mutations_type_plot", dpi=300)
    plt.close()

    g9 = sns.catplot("label", "frac_and_weight", data=data_filter_uc, hue="Next", order=label_order, palette="tab20",
                     hue_order=context_order_uc, estimator=weighted_varaint, orient="v", dodge=True, kind="point",
                     col="Type", join=False, col_order=type_order2)
    g9.set_axis_labels("", "Variant Frequency")
    g9.set(yscale='log')
    g9.set(ylim=(10 ** -7, 10 ** -3))
    g9.set_xticklabels(rotation=45)
    # plt.show()
    g9.savefig(output_dir + "/UC_Context_point_plot", dpi=300)
    plt.close()

    data_filter_ag_grouped = data_filter_ag.groupby(["ADAR_like", "label", "Type", "Pos", "Protein"])["frac_and_weight"].agg(lambda x: weighted_varaint(x))
    data_filter_ag_grouped = data_filter_ag_grouped.reset_index()
    data_filter_ag_grouped = data_filter_ag_grouped.rename(columns={"frac_and_weight": "Frequency"})
    data_filter_ag_grouped["Frequency"] = data_filter_ag_grouped["Frequency"].astype(float)
    print(data_filter_ag_grouped.to_string())

    data_filter_ag_grouped_silent = data_filter_ag_grouped[data_filter_ag_grouped["Type"] == "Synonymous"]
    data_filter_ag_grouped_silent = data_filter_ag_grouped_silent[data_filter_ag_grouped_silent["label"] == "Cell Cultureֿ\nControl"]
    position_g = sns.scatterplot("Pos", "Frequency", data=data_filter_ag_grouped_silent, hue="Protein",
                                 palette="tab10", style="ADAR_like", style_order=[True, False], legend="full")

    # position_g.set_axis_labels("", "Variant Frequency")
    position_g.set_yscale('log')
    position_g.set_ylim(10 ** -4, 10 ** -1)
    # position_g.set(xlim=(3500, 7500))
    position_g.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(output_dir + "/position_cells.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    main()

