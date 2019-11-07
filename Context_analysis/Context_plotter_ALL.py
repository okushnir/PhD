
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker

import numpy as np
import seaborn as sns; sns.set(font_scale=1.2)

# print(plt.style.available)

def main():
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/"
    output_dir = input_dir + "plots/"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    data_mutations = pd.read_csv(input_dir + "data_mutation.csv")

    data_mutations = data_mutations.rename(columns={"label": "Sample"})
    # print(data_mutations.to_string())
    mutation_order = ["A>G", "U>C", "C>U", "G>A"]

    sample_order = ["RVB14-RNA Control", "RVB14-p2", "RVB14-p5", "RVB14-p7", "RVB14-p8", "RVB14-p10", "RVB14-p12",
                    "CVB3-RNA Control", "CVB3-p2", "CVB3-p5", "CVB3-p8", "CVB3-p10", "CVB3-p12","ATCG1-RV-minus",
                    "ATCG1-RV-plus", "FLNA-RV-plus"]
    type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]

    plt.style.use('classic')

# Transitions mutations plot
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        transiotion_plot = sns.catplot(x="Mutation", y="Frequency", hue="Sample", data=data_mutations, palette="tab20",
                                       order=mutation_order, hue_order=sample_order, flierprops={"marker": "."},
                                       kind="box")
    transiotion_plot.set(yscale='log')
    transiotion_plot.set_xlabels("")
    transiotion_plot._legend.set_title('')
    transiotion_plot.fig.suptitle("Transitions variant frequency in all samples")

    plt.savefig(output_dir + "Transitions_variant_frequency.png", dpi=600)
    plt.close()

# facet plot for all mutations

    with sns.plotting_context(rc={"legend.fontsize": 6}):
        cat_plot = sns.catplot(x="Mutation", y="Frequency", hue="Sample", data=data_mutations, palette="tab20",
                               order=mutation_order, hue_order=sample_order, col="Type", kind="box",
                               flierprops={"marker": "."}, col_order=type_order)

    cat_plot.set(yscale='log')
    cat_plot.set_xlabels("")
    cat_plot._legend.set_title('')
    cat_plot.fig.suptitle("Transitions variant frequency in all samples", y=0.99)
    # cat_plot.legend(loc="upper right")

    plt.savefig(output_dir + "Facet_Transitions_variant_frequency.png", dpi=600)
    plt.close()

# Just silent mutations
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        silent_plot = sns.catplot(x="Mutation", y="Frequency", hue="Sample",
                                  data=data_mutations[data_mutations["Type"] =="Synonymous"], palette="tab20",
                                  order=mutation_order, hue_order=sample_order, flierprops={"marker": "."}, kind="box")
    silent_plot.set(yscale='log')
    silent_plot._legend.set_title('')
    silent_plot.set(xlabel="")
    silent_plot.fig.suptitle("Transitions variant frequency, only silent mutations")

    plt.savefig(output_dir + "Silent_Transition_variant_frequency.png", dpi=300)
    plt.close()

    data_context = pd.read_csv(input_dir + "data_XpA_by_mutation.csv")

    data_context['Prev'].replace('AA', 'ApA', inplace=True)
    data_context['Prev'].replace('UA', 'UpA', inplace=True)
    data_context['Prev'].replace('CA', 'CpA', inplace=True)
    data_context['Prev'].replace('GA', 'GpA', inplace=True)

    data_context['Type'].replace("Synonymous", "Silent", inplace=True)
    data_context['Type'].replace("Non-Synonymous", "Missense", inplace=True)
    # data_context.rename(columns={"Freq": "Frequency"}, inplace=True)


    context_order = ["UpA", "ApA", "CpA", "GpA"]



# 5’ neighbors preferences

    # 2 plots in one fig - silent and missense
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        context_plot = sns.catplot(x="Prev", y="Frequency", hue="label", data=data_context, palette="tab20", order=context_order,
                                   hue_order=sample_order, col="Type", kind="box", flierprops={"marker": "."})

    context_plot.set(yscale='log')
    context_plot._legend.set_title('')
    context_plot.set(xlabel='')

    context_plot.fig.suptitle("A>G 5’ neighbors preferences in all samples", y=0.99)

    plt.savefig(output_dir + "2_plots_ADAR_Context.png", dpi=300)
    plt.close()

    with sns.plotting_context(rc={"legend.fontsize": 6}):
        context_plot = sns.catplot(x="Prev", y="Frequency", hue="label", data=data_context,#[data_context["Type"] == "Silent"],
                                   palette="tab20", order=context_order, hue_order=sample_order, kind="box",
                                   flierprops={"marker": "."})

    context_plot.set(yscale='log')
    context_plot.set_xlabels("")
    context_plot._legend.set_title('')
    context_plot.fig.suptitle("A>G 5’ neighbors preferences in all samples", y=0.99)

    context_plot.savefig(output_dir + "ADAR_Context_both_types.png", dpi=300)
    plt.close()

    with sns.plotting_context(rc={"legend.fontsize": 6}):
        context_plot = sns.catplot(x="Prev", y="Frequency", hue="label", data=data_context[data_context["Type"] == "Silent"],
                                   palette="tab20", order=context_order, hue_order=sample_order, kind="box",
                                   flierprops={"marker": "."})

    context_plot.set(yscale='log')
    context_plot.set_xlabels("")
    context_plot._legend.set_title('')
    context_plot.fig.suptitle("A>G 5’ neighbors preferences in all samples, only silent mutations", y=0.99)

    context_plot.savefig(output_dir + "ADAR_Context_silent.png", dpi=300)
    plt.close()



# anti-sense neighbors preferences
    data_context_UC = pd.read_csv(input_dir + "data_UpX_by_mutation.csv")


    data_context_UC['Next'].replace('UU', 'UpU', inplace=True)
    data_context_UC['Next'].replace('UA', 'UpA', inplace=True)
    data_context_UC['Next'].replace('UC', 'UpC', inplace=True)
    data_context_UC['Next'].replace('UG', 'UpG', inplace=True)

    data_context_UC['Type'].replace("Synonymous", "Silent", inplace=True)
    data_context_UC['Type'].replace("Non-Synonymous", "Missense", inplace=True)
    # data_context_UC = data_context_UC.rename(columns={"Freq": "Frequency"})
    data_context_UC = data_context_UC.rename(columns={"label": "Sample"})

    context_order_UC = ["UpU", "UpA", "UpC", "UpG"]
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        context_UC_plot = sns.catplot(x="Next", y="Frequency", hue="Sample",
                                      data=data_context_UC[data_context_UC["Type"] == "Silent"],
                                      palette="tab20", order=context_order_UC, hue_order=sample_order,
                                      flierprops={"marker": "."}, kind="box")
    context_UC_plot.set(yscale='log')

    context_UC_plot.set_xlabels("")
    context_UC_plot._legend.set_title('')

    context_UC_plot.fig.suptitle("U>C 3'(anti-sense) neighbors preferences in all samples\nonly silent mutations", y=0.99)

    plt.savefig(output_dir + "3'_ADAR_Context.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()