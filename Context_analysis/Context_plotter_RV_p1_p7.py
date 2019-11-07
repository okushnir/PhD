
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker

import numpy as np
import seaborn as sns; sns.set(font_scale=1.2)

# print(plt.style.available)

def main():
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/"
    output_dir = input_dir + "plots/"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    data_mutations = pd.read_csv(input_dir + "data_mutation.csv")
    # print(data_mutations.to_string())
    mutation_order = ["A>G", "U>C", "C>U", "G>A"]

    # data_mutations['label'].replace('RV-p71', 'p7', inplace=True)
    # data_mutations['label'].replace('RV-IVT', 'RNA Control', inplace=True)
    # data_mutations['label'].replace('RV-p11', 'p1', inplace=True)
    sample_order = ["RVB14-RNA Control", "RVB14-p2", "RVB14-p5", "RVB14-p7", "RVB14-p8", "RVB14-p10", "RVB14-p12"]

    plt.style.use('classic')

    # data_transitions = data_mutations
    # data_transitions["Transitions"] = data_transitions.Mutation.str.contains("A>G") | data_transitions.Mutation.str.contains("U>C") | data_transitions.Mutation.str.contains("C>U") | data_transitions.Mutation.str.contains("G>A")
    # data_transitions = data_transitions[data_transitions["Transitions"] == True]
    #
    # data_transitions["Samples"] = data_transitions.label.str.contains('RNA Control') | data_transitions.label.str.contains('p7') | data_transitions.label.str.contains("p1")
    # data_transitions = data_transitions[data_transitions["Samples"] == True]

# Transitions mutations plot
    transiotion_plot = sns.boxplot(x="Mutation", y="Frequency", hue="label", data=data_mutations, palette="Set2",
                                   order=mutation_order, hue_order=sample_order, flierprops={"marker": "."})
    transiotion_plot.set_yscale('log')
    transiotion_plot.legend().set_title('')
    transiotion_plot.set(xlabel="")
    transiotion_plot.set_title("Transitions variant frequency in RV")
    transiotion_plot.legend(loc="upper right")

    plt.savefig(output_dir + "Transitions_variant_frequency_in_RV.png", dpi=300)
    plt.close()

# facet plot for all mutations
    data_mutations = data_mutations.rename(columns={"label": "Sample"})
    # data_mutations = data_mutations.rename(columns={"Freq": "Frequency"})
    cat_plot = sns.catplot(x="Mutation", y="Frequency", hue="Sample", data=data_mutations, palette="Set2", order=mutation_order,
                               hue_order=sample_order, col="Type", kind="box", flierprops={"marker": "."})

    cat_plot.set(yscale='log')
    cat_plot.set_xlabels("")

    cat_plot.fig.suptitle("Transitions variant frequency in RV", y=0.99)
    # cat_plot.legend(loc="upper right")

    plt.savefig(output_dir + "Facet_Transitions_variant_frequency_in_RV.png", dpi=300)
    plt.close()

# Just silent mutations

    silent_plot = sns.boxplot(x="Mutation", y="Frequency", hue="Sample", data=data_mutations[data_mutations["Type"] ==
                                                                                            "Synonymous"], palette="Set2",
                                   order=mutation_order, hue_order=sample_order, flierprops={"marker": "."})
    silent_plot.set_yscale('log')
    silent_plot.legend().set_title('')
    silent_plot.set(xlabel="")
    silent_plot.set_title("Transitions variant frequency in RV, only silent mutations")
    silent_plot.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(output_dir + "Silent_Transition_variant_frequency_in_RV.png", dpi=300)
    plt.close()

    data_context = pd.read_csv(input_dir + "data_XpA_by_mutation.csv")
    # data_context['label'].replace('RV-p71', 'p7', inplace=True)
    # data_context['label'].replace('RV-IVT', 'RNA Control', inplace=True)
    # data_context['label'].replace('RV-p11', 'p1', inplace=True)

    data_context['Prev'].replace('AA', 'ApA', inplace=True)
    data_context['Prev'].replace('UA', 'UpA', inplace=True)
    data_context['Prev'].replace('CA', 'CpA', inplace=True)
    data_context['Prev'].replace('GA', 'GpA', inplace=True)

    data_context['Type'].replace("Synonymous", "Silent", inplace=True)
    data_context['Type'].replace("Non-Synonymous", "Missense", inplace=True)
    data_context.rename(columns={"Freq": "Frequency"}, inplace=True)


    context_order = ["UpA", "ApA", "GpA", "CpA"]
    type_order = ["Silent", "Missense", "Premature Stop Codon"]

    # data_context_silent = data_context
    # data_context_silent["silent"] = data_context_silent.Type.str.contains("Silent")
    # data_context_silent = data_context_silent[data_context_silent["silent"] == True]


# 5’ neighbors preferences

    # 2 plots in one fig - silent and missense
    # print(data_context.to_string())
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        context_plot = sns.catplot(x="Prev", y="Frequency", hue="label", data=data_context, palette="set2", order=context_order,
                                   hue_order=sample_order, col="Type", flierprops={"marker": "."}, kind="box")


    context_plot.set(yscale='log')
    context_plot._legend.set_title('')
    context_plot.set(xlabel='')


    context_plot.fig.suptitle("A>G 5’ neighbors preferences in cell culture infected with RV", y=0.99)

    plt.savefig(output_dir + "2_plots_ADAR_Context.png", dpi=300)
    plt.close()


    context_plot = sns.boxplot(x="Prev", y="Frequency", hue="label", data=data_context[data_context["Type"] == "Silent"],
                               palette="Set2", order=context_order, hue_order=sample_order, flierprops={"marker": "."})
    context_plot.set_yscale('log')
    context_plot.legend().set_title('')
    context_plot.set_xlabel('')
    # context_plot.set_ylabel('Frequency')
    context_plot.set_title("A>G 5’ neighbors preferences in cell culture infected with RV")
    context_plot.legend(loc="upper right")
    plt.savefig(output_dir + "ADAR_Context.png", dpi=300)
    plt.close()



# anti-sense neighbors preferences
    data_context_UC = pd.read_csv(input_dir + "data_UpX_by_mutation.csv")
    # data_context_UC['label'].replace('RV-p71', 'p7', inplace=True)
    # data_context_UC['label'].replace('RV-IVT', 'RNA Control', inplace=True)
    # data_context_UC['label'].replace('RV-p11', 'p1', inplace=True)

    data_context_UC['Next'].replace('UU', 'UpU', inplace=True)
    data_context_UC['Next'].replace('UA', 'UpA', inplace=True)
    data_context_UC['Next'].replace('UC', 'UpC', inplace=True)
    data_context_UC['Next'].replace('UG', 'UpG', inplace=True)

    data_context_UC['Type'].replace("Synonymous", "Silent", inplace=True)
    data_context_UC['Type'].replace("Non-Synonymous", "Missense", inplace=True)

    context_order_UC = ["UpU", "UpA", "UpC", "UpG"]

    context_UC_plot = sns.boxplot(x="Next", y="Frequency", hue="label",
                                  data=data_context_UC[data_context_UC["Type"] == "Silent"],
                                  palette="Set2", order=context_order_UC, hue_order=sample_order, flierprops={"marker": "."})
    context_UC_plot.set_yscale('log')
    context_UC_plot.legend().set_title('')
    context_UC_plot.set_xlabel('')
    # context_UC_plot.set_ylabel('Frequency')
    context_UC_plot.set_title("RV in cell culture anti-sense neighbors preferences")
    context_UC_plot.legend(loc="upper right")
    plt.savefig(output_dir + "3'_ADAR_Context.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()