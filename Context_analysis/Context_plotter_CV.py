
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker

import numpy as np
import seaborn as sns
from statannot import add_stat_annotation

# print(plt.style.available)

def main():
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/CVB3/"
    output_dir = input_dir + "plots_q38/"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    data_mutations = pd.read_csv(input_dir + "q38_data_mutation.csv")
    # print(data_mutations.to_string())
    mutation_order = ["A>G", "U>C", "C>U", "G>A", "A>U", "A>C", "U>A", "U>G", "C>A", "C>G", "G>U", "G>C"]
    transition_order = ["A>G", "U>C", "C>U", "G>A"]
    sample_order = ["CVB3-RNA Control", "CVB3-p2", "CVB3-p5", "CVB3-p8", "CVB3-p10", "CVB3-p12"]
    type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]

    plt.style.use('classic')

    #stat
    data_p = data_mutations[data_mutations["label"] =="CVB3-p12"]

    transiotion_plot = sns.boxplot(x="Mutation", y="Frequency", hue="label", data=data_mutations, palette="Set2",
                                   order=transition_order, flierprops={"marker": "."}, hue_order=sample_order)
    transiotion_plot.set_yscale('log')
    add_stat_annotation(transiotion_plot, data=data_mutations, x="Mutation", y="Frequency", order=transition_order,
                        boxPairList=[("A>G", "U>C"), ("A>G", "C>U"), ("A>G", "G>A")],
                        test='Mann-Whitney', textFormat='star', loc='inside', verbose=2)
    transiotion_plot.legend(bbox_to_anchor=(1.05, 0.5), loc=3, borderaxespad=0., fontsize='small')
    plt.tight_layout()
    # plt.show()
    plt.savefig(output_dir + "Stat_Transitions_variant_frequency_in_CV.png", dpi=300)
    plt.close()

    # Transitions mutations plot
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        transiotion_plot = sns.catplot(x="Mutation", y="Frequency", hue="label", data=data_mutations, palette="Set2",
                                       order=mutation_order, hue_order=sample_order, flierprops={"marker": "."},
                                       kind="box")
    transiotion_plot.set(yscale='log')
    transiotion_plot._legend.set_title('')
    transiotion_plot.set(xlabel="")
    transiotion_plot.fig.suptitle("Transitions variant frequency in CV", y=0.99)

    plt.savefig(output_dir + "Transitions_variant_frequency_in_CV.png", dpi=300)
    plt.close()

    # facet plot for all mutations
    data_mutations = data_mutations.rename(columns={"label": "Sample"})
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        cat_plot = sns.catplot(x="Mutation", y="Frequency", hue="Sample", data=data_mutations, palette="Set2",
                               order=transition_order, hue_order=sample_order, col="Type", kind="box",
                               flierprops={"marker": "."}, col_order=type_order)

    cat_plot.set(yscale='log')
    cat_plot.set_xlabels("")
    cat_plot._legend.set_title('')
    cat_plot.fig.suptitle("Transitions variant frequency in CV", y=0.99)
    plt.savefig(output_dir + "Facet_Transitions_variant_frequency_in_CV.png", dpi=300)
    plt.close()

    # Just silent mutations
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        silent_plot = sns.catplot(x="Mutation", y="Frequency", hue="Sample",
                                  data=data_mutations[data_mutations["Type"] == "Synonymous"], palette="Set2",
                                  order=transition_order, hue_order=sample_order, flierprops={"marker": "."}, kind="box")
    silent_plot.set(yscale='log')
    silent_plot._legend.set_title('')
    silent_plot.set(xlabel="")
    silent_plot.fig.suptitle("Transitions variant frequency in CV, only silent mutations", y=0.99)

    plt.savefig(output_dir + "Silent_Transition_variant_frequency_in_CV.png", dpi=300)
    plt.close()

    data_context = pd.read_csv(input_dir + "q38_data_XpA_by_mutation.csv")

    data_context['Prev'].replace('AA', 'ApA', inplace=True)
    data_context['Prev'].replace('UA', 'UpA', inplace=True)
    data_context['Prev'].replace('CA', 'CpA', inplace=True)
    data_context['Prev'].replace('GA', 'GpA', inplace=True)

    data_context['Type'].replace("Synonymous", "Silent", inplace=True)
    data_context['Type'].replace("Non-Synonymous", "Missense", inplace=True)

    context_order = ["UpA", "ApA", "CpA", "GpA"]

    # 5’ neighbors preferences

    # 2 plots in one fig - silent and missense
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        context_plot = sns.catplot(x="Prev", y="Frequency", hue="label", data=data_context, palette="Set2",
                                   order=context_order,
                                   hue_order=sample_order, col="Type", kind="box", flierprops={"marker": "."})

    context_plot.set(yscale='log')
    context_plot._legend.set_title('')
    context_plot.set(xlabel='')

    context_plot.fig.suptitle("A>G 5’ neighbors preferences in cell culture infected with CV", y=0.99)

    plt.savefig(output_dir + "2_plots_ADAR_Context.png", dpi=300)
    plt.close()

    with sns.plotting_context(rc={"legend.fontsize": 6}):
        context_plot = sns.catplot(x="Prev", y="Frequency", hue="label",
                                   data=data_context[data_context["Type"] == "Silent"],
                                   palette="Set2", order=context_order, hue_order=sample_order, kind="box",
                                   flierprops={"marker": "."})

    context_plot.set(yscale='log')
    context_plot.set_xlabels("")
    context_plot._legend.set_title('')
    context_plot.fig.suptitle("A>G 5’ neighbors preferences in CV, only silent mutations", y=0.99)

    context_plot.savefig(output_dir + "ADAR_Context_silent.png", dpi=300)
    plt.close()

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

    data_context_UC = pd.read_csv(input_dir + "q38_data_UpX_by_mutation.csv")
    data_context_UC['Next'].replace('UU', 'UpU', inplace=True)
    data_context_UC['Next'].replace('UA', 'UpA', inplace=True)
    data_context_UC['Next'].replace('UC', 'UpC', inplace=True)
    data_context_UC['Next'].replace('UG', 'UpG', inplace=True)

    data_context_UC['Type'].replace("Synonymous", "Silent", inplace=True)
    data_context_UC['Type'].replace("Non-Synonymous", "Missense", inplace=True)

    context_order_UC = ["UpU", "UpA", "UpC", "UpG"]
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        context_UC_plot = sns.catplot(x="Next", y="Frequency", hue="label",
                                      data=data_context_UC[data_context_UC["Type"] == "Silent"],
                                      palette="Set2", order=context_order_UC, hue_order=sample_order,
                                      flierprops={"marker": "."}, kind="box")
    context_UC_plot.set(yscale='log')
    context_UC_plot._legend.set_title('')
    context_UC_plot.set_xlabels("")

    context_UC_plot.fig.suptitle("U>C 3' (Anti-sense) neighbors preferences in CV only silent mutations", y=0.99)
    plt.savefig(output_dir + "3'_ADAR_Context.png", dpi=300)
    plt.close()

    data_mutations = data_mutations.rename(columns={"label": "Sample"})
    data_mutations_PMSC = data_mutations[data_mutations["Mutation"] == "C>U"]

    with sns.plotting_context(rc={"legend.fontsize": 6}):
        cat_plot_mut = sns.catplot(x="Mutated_codon", y="Frequency", hue="Sample",
                                   data=data_mutations_PMSC[data_mutations_PMSC["Type"] == "Premature Stop Codon"],
                                   palette="Set2",
                                   hue_order=sample_order, col="Type", kind="box", flierprops={"marker": "."})

    cat_plot_mut.set(yscale='log')
    cat_plot_mut.set_xlabels("")
    cat_plot_mut._legend.set_title('')
    cat_plot_mut.fig.suptitle("Premature Stop Codon C>T mutated codons in CV", y=0.99)

    plt.savefig(output_dir + "Facet_Transitions_variant_frequency_in_CV_stop.png", dpi=300)
    plt.close()

    sns.set(font_scale=0.8)

    # All Mutations
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        all_plot = sns.catplot(x="Mutation", y="Frequency", hue="Sample", data=data_mutations, palette="Set2",
                               order=mutation_order, hue_order=sample_order, flierprops={"marker": "."},
                               kind="box", col="Type", col_order=type_order)
    all_plot.set(yscale='log')
    all_plot._legend.set_title('')
    all_plot.set(xlabel="")
    all_plot.fig.suptitle("All Mutations variant frequencies in CV", y=0.99)

    plt.savefig(output_dir + "All_Mutations_variant_frequencies_in_CV.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    main()