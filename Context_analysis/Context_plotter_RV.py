
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker

import numpy as np
import seaborn as sns
from statannot import add_stat_annotation
from scipy import stats

# print(plt.style.available)

def main():
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/"
    output_dir = input_dir + "20190811_plots/"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    data_mutations = pd.read_csv(input_dir + "20190811_q38_data_mutation.csv")
    # print(data_mutations.to_string())
    mutation_order = ["A>G", "U>C", "C>U", "G>A", "A>U", "A>C", "U>A", "U>G", "C>A", "C>G", "G>U", "G>C"]
    transition_order = ["A>G", "U>C", "C>U", "G>A"]
    sample_order = ["RVB14-RNA Control", "RVB14-p2", "RVB14-p5", "RVB14-p7", "RVB14-p8", "RVB14-p10", "RVB14-p12"] #"RVB14-Next-RNA Control", "RVB14-p1",
    type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]
    data_control = data_mutations[data_mutations["label"] == "RVB14-RNA Control"]
    data_next = data_mutations[data_mutations["label"] == "RVB14-Next-RNA Control"]
    data_mutations = data_mutations[data_mutations["pval"] < 0.01]
    data_mutations = data_mutations[data_mutations["Prob"] > 0.95]
    data_mutations = pd.concat([data_mutations, data_control, data_next], sort=False)

    plt.style.use('classic')

    # All Mutations
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        all_plot = sns.catplot(x="Mutation", y="Frequency", hue="label", data=data_mutations, palette="tab20",
                               order=mutation_order, hue_order=sample_order, flierprops={"marker": "."},
                               kind="box", col="Type", col_order=type_order)
    all_plot.set(yscale='log')
    all_plot._legend.set_title('')
    all_plot.set(xlabel="")
    all_plot.fig.suptitle("All Mutations variant frequencies in RV", y=0.99)

    plt.savefig(output_dir + "All_Mutations_variant_frequencies_in_RV.png", dpi=300)
    plt.close()

    # mutation = "U>G"
    # data_mutations_stop = data_mutations[data_mutations["Type"] == "Premature Stop Codon"]
    # data_mutations_stop_cu = data_mutations_stop[data_mutations_stop["Mutation"]=="%s" % mutation]
    # # ata_mutations_stop_cu = data_mutations_stop_cu[data_mutations_stop_cu["Rank"] > 0]
    # data_mutations_stop_cu["passage"] = data_mutations_stop_cu["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    # data_mutations_stop_cu["passage"] = np.where(data_mutations_stop_cu["passage"] == "RNA Control", 0, data_mutations_stop_cu["passage"])
    # data_mutations_stop_cu["passage"] = data_mutations_stop_cu["passage"].astype(int)
    # data_mutations_stop_cu = data_mutations_stop_cu[data_mutations_stop_cu["label"] != "RVB14-p7"]
    # # data_mutations_stop_cu["Frequency"] = np.log(data_mutations_stop_cu["Frequency"])
    #
    # slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(data_mutations_stop_cu['passage'], data_mutations_stop_cu['Frequency'])
    # print("pval: %s" % p_value1)
    #
    # g1 = sns.regplot(x="passage", y="Frequency", data=data_mutations_stop_cu, line_kws={'label':"y={0:.3g}x+{1:.7g}".format(slope1, intercept1)})
    # g1.legend()
    # g1.legend(loc=2)
    # g1.set(title="Mutation rate")
    #
    # g1.set_yscale("log")
    # plt.show()
    # plt.savefig(output_dir + "Mutation_rate_%s.png" % mutation, dpi=300)
    # plt.close()

    #stat
    passage = "5"
    data_p = data_mutations[data_mutations["label"] == "RVB14-p%s" % passage]

    transiotion_plot = sns.boxplot(x="Mutation", y="Frequency", hue="label", data=data_p, palette="tab20",
                                   order=transition_order, flierprops={"marker": "."})#, hue_order=sample_order)
    transiotion_plot.set_yscale('log')
    add_stat_annotation(transiotion_plot, data=data_p, x="Mutation", y="Frequency", order=transition_order,
                        boxPairList=[("A>G", "U>C"), ("A>G", "C>U"), ("A>G", "G>A")],
                        test='Mann-Whitney', textFormat='star', loc='inside', verbose=2)
    transiotion_plot.legend(bbox_to_anchor=(1.05, 0.5), loc=3, borderaxespad=0., fontsize='small')
    plt.tight_layout()
    plt.savefig(output_dir + "p%s_Stat_Transitions_variant_frequency_in_RV.png" % passage, dpi=300)
    plt.close()

# Transitions mutations plot
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        transiotion_plot = sns.catplot(x="Mutation", y="Frequency", hue="label", data=data_mutations, palette="tab20",
                                       order=transition_order, hue_order=sample_order, flierprops={"marker": "."},
                                       kind="box")
    transiotion_plot.set(yscale='log')
    transiotion_plot._legend.set_title('')
    transiotion_plot.set(xlabel="")
    transiotion_plot.fig.suptitle("Transitions variant frequency in RV", y=0.99)

    plt.savefig(output_dir + "Transitions_variant_frequency_in_RV.png", dpi=300)
    plt.close()

# facet plot for all mutations
    data_mutations = data_mutations.rename(columns={"label": "Sample"})
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        cat_plot = sns.catplot(x="Mutation", y="Frequency", hue="Sample", data=data_mutations, palette="tab20",
                               order=transition_order, hue_order=sample_order, col="Type", kind="box",
                               flierprops={"marker": "."}, col_order=type_order)

    cat_plot.set(yscale='log')
    cat_plot.set_xlabels("")
    cat_plot._legend.set_title('')
    cat_plot.fig.suptitle("Transitions variant frequency in RV", y=0.99)
    plt.savefig(output_dir + "Facet_Transitions_variant_frequency_in_RV.png", dpi=300)
    plt.close()

# Just silent mutations
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        silent_plot = sns.catplot(x="Mutation", y="Frequency", hue="Sample",
                                  data=data_mutations[data_mutations["Type"] =="Synonymous"], palette="tab20",
                                  order=transition_order, hue_order=sample_order, flierprops={"marker": "."}, kind="box")
    silent_plot.set(yscale='log')
    silent_plot._legend.set_title('')
    silent_plot.set(xlabel="")
    silent_plot.fig.suptitle("Transitions variant frequency in RV, only silent mutations", y=0.99)

    plt.savefig(output_dir + "Silent_Transition_variant_frequency_in_RV.png", dpi=300)
    plt.close()

    data_context = pd.read_csv(input_dir + "20190811_q38_data_XpA_by_mutation.csv")

    data_context['Prev'].replace('AA', 'ApA', inplace=True)
    data_context['Prev'].replace('UA', 'UpA', inplace=True)
    data_context['Prev'].replace('CA', 'CpA', inplace=True)
    data_context['Prev'].replace('GA', 'GpA', inplace=True)

    data_context['Type'].replace("Synonymous", "Silent", inplace=True)
    data_context['Type'].replace("Non-Synonymous", "Missense", inplace=True)
    data_context = data_context[data_context["pval"] < 0.01]
    data_context = data_context[data_context["Prob"] > 0.95]

    context_order = ["UpA", "ApA", "CpA", "GpA"]

    # 5’ neighbors preferences

    # 2 plots in one fig - silent and missense
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        context_plot = sns.catplot(x="Prev", y="Frequency", hue="label", data=data_context, palette="tab20",
                                   order=context_order, hue_order=sample_order, col="Type", kind="box",
                                   flierprops={"marker": "."})

    context_plot.set(yscale='log')
    context_plot._legend.set_title('')
    context_plot.set(xlabel='')

    context_plot.fig.suptitle("A>G 5’ neighbors preferences in cell culture infected with RV", y=0.99)

    plt.savefig(output_dir + "2_plots_ADAR_Context.png", dpi=300)
    plt.close()

    with sns.plotting_context(rc={"legend.fontsize": 6}):
        context_plot = sns.catplot(x="Prev", y="Frequency", hue="label", data=data_context[data_context["Type"] == "Silent"],
                                   palette="tab20", order=context_order, hue_order=sample_order, kind="box",
                                   flierprops={"marker": "."})

    context_plot.set(yscale='log')
    context_plot.set_xlabels("")
    context_plot._legend.set_title('')
    context_plot.fig.suptitle("A>G 5’ neighbors preferences in RV, only silent mutations", y=0.99)

    context_plot.savefig(output_dir + "ADAR_Context_silent.png", dpi=300)
    plt.close()

    # data_context = data_context[data_context["Protein"] != "2A"]
    # sns.set_context("talk")
    # ax = sns.scatterplot(x="Pos", y="Frequency", hue="Protein", data=data_context[data_context["Type"] == "Silent"]
    #                      , palette="tab20",style="Prev", style_order=context_order) #
    # ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize='xx-small')
    #
    # ax.set_yscale('log')
    # ax.set(ylim=(10 ** -6, 1))
    # ax.set(xlabel='')
    # plt.tight_layout()
    # plt.savefig(output_dir + "Pos_silent.png", dpi=300)
    # plt.close()

# anti-sense neighbors preferences

    data_context_UC = pd.read_csv(input_dir + "20190811_q38_data_UpX_by_mutation.csv")
    # data_context_UC['label'].replace('RV-p71', 'p7', inplace=True)
    # data_context_UC['label'].replace('RV-IVT', 'RNA Control', inplace=True)
    # data_context_UC['label'].replace('RV-p11', 'p1', inplace=True)

    data_context_UC['Next'].replace('UU', 'UpU', inplace=True)
    data_context_UC['Next'].replace('UA', 'UpA', inplace=True)
    data_context_UC['Next'].replace('UC', 'UpC', inplace=True)
    data_context_UC['Next'].replace('UG', 'UpG', inplace=True)

    data_context_UC['Type'].replace("Synonymous", "Silent", inplace=True)
    data_context_UC['Type'].replace("Non-Synonymous", "Missense", inplace=True)
    data_context_UC = data_context_UC[data_context_UC["pval"] < 0.01]
    data_context_UC = data_context_UC[data_context_UC["Prob"] > 0.95]

    context_order_UC = ["UpU", "UpA", "UpC", "UpG"]
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        context_UC_plot = sns.catplot(x="Next", y="Frequency", hue="label",
                                      data=data_context_UC[data_context_UC["Type"] == "Silent"],
                                      palette="tab20", order=context_order_UC, hue_order=sample_order,
                                      flierprops={"marker": "."}, kind="box")
    context_UC_plot.set(yscale='log')
    context_UC_plot._legend.set_title('')
    context_UC_plot.set_xlabels("")

    context_UC_plot.fig.suptitle("U>C 3' (Anti-sense) neighbors preferences in RV only silent mutations", y=0.99)
    plt.savefig(output_dir + "3'_ADAR_Context.png", dpi=300)
    plt.close()

    sns.set(font_scale=0.8)


if __name__ == "__main__":
    main()
