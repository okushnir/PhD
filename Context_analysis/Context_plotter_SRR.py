
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation

def plots_for_srr(input_dir, output_dir, virus_path, virus):
    data_mutations = pd.read_csv(input_dir + virus_path + "data_mutation.csv")
    data_mutations["Mutation"] = data_mutations["Mutation"].apply(lambda x: x.split("->")[0] + ">" + x.split("->")[1])
    mutation_order = ["A>G", "U>C", "C>U", "G>A", "A>U", "A>C", "U>A", "U>G", "C>A", "C>G", "G>U", "G>C"]
    transition_order = ["A>G", "U>C", "C>U", "G>A"]
    type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]

    # data_mutations = data_mutations[data_mutations["pval"] < 0.01]
    data_mutations = data_mutations[data_mutations["Prob"] > 0.95]

    # plt.style.use('classic')
    sns.set_style("ticks")
    mypalette = ["#9b59b6", "#3498db"]

    # All Mutations
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        all_plot = sns.catplot(x="Mutation", y="Frequency", hue="label", data=data_mutations, palette="tab20",
                               order=mutation_order, flierprops={"marker": "."},
                               kind="box", col="Type", col_order=type_order)
    all_plot.set(yscale='log')
    all_plot._legend.set_title('')
    all_plot.set(xlabel="")
    all_plot.fig.suptitle("All Mutations variant frequencies in %s" % virus, y=0.99)

    plt.savefig(output_dir + "All_Mutations_variant_frequencies_in_%s.png" % virus, dpi=300)
    plt.close()

    # Transitions mutations plot
    with sns.plotting_context(rc={"legend.fontsize": 6}):
        transiotion_plot = sns.catplot(x="Mutation", y="Frequency", hue="label", data=data_mutations, palette="tab20",
                                       order=transition_order, flierprops={"marker": "."},
                                       kind="box")
    transiotion_plot.set(yscale='log')
    transiotion_plot._legend.set_title('')
    transiotion_plot.set(xlabel="")
    transiotion_plot.fig.suptitle("Transitions variant frequency in %s" % virus, y=0.99)

    plt.savefig(output_dir + "Transitions_variant_frequency_in_%s.png" % virus, dpi=300)
    plt.close()

##Context
    data_context = pd.read_csv(input_dir + virus_path + "data_XpA_by_mutation.csv")
    if virus == "EnteroA":
        data_context = data_context.loc[data_context.Organism != "Coxsackievirus A16"]
    data_context["Mutation"] = data_context["Mutation"].apply(lambda x: x.split("->")[0] + ">" + x.split("->")[1])
    data_context = data_context.rename(columns={"Freq": "Frequency"})
    data_context['Prev'].replace('AA', 'ApA', inplace=True)
    data_context['Prev'].replace('UA', 'UpA', inplace=True)
    data_context['Prev'].replace('CA', 'CpA', inplace=True)
    data_context['Prev'].replace('GA', 'GpA', inplace=True)
    data_context['Type'].replace('Synonymous', 'Silent', inplace=True)
    data_context['Type'].replace('Non-Synonymous', 'Missense', inplace=True)
    # data_context_rv = data_context_rv[data_context_rv["pval"] < 0.01]
    data_context = data_context[data_context["Prob"] > 0.95]

    context_order = ["UpA", "ApA", "CpA", "GpA"]

    data_adar = data_context.loc[data_context.Type == "Silent"]
    # print(type(data_adar["Context"]))
    data_adar["ADAR_like"] = data_adar.Prev.str.contains('UpA') | data_adar.Prev.str.contains('ApA')
    no_organism = data_adar.Organism.value_counts()



    # 5’ neighbors preferences
    #stat
    if virus == "EnteroA":
        context_stat_plot = sns.catplot("ADAR_like", "Frequency", data=data_adar, palette=mypalette, col="Organism",
                                        flierprops={"marker": "."}, kind="box", col_wrap=int(len(no_organism)/2))
    else:
        context_stat_plot = sns.catplot("ADAR_like", "Frequency", data=data_adar, palette=mypalette, col="Organism",
                                        flierprops={"marker": "."}, kind="box", col_wrap=len(no_organism))

    context_stat_plot.set(yscale='log')
    context_stat_plot.set(ylim=(10 ** -5, 0))
    context_stat_plot.set_axis_labels("", "Variant Frequency")
    organism_list = data_adar["Organism"].unique()
    for ax in context_stat_plot.axes.flat:
        print(ax.get_title())
        for organism in organism_list:
            if ax.get_title().split(" = ")[-1] == organism:
                ax.set_title(organism, pad=20)
                add_stat_annotation(ax, data=data_adar[data_adar["Organism"] == organism], x="ADAR_like", y="Frequency",
                                boxPairList=[(True, False)], test='Mann-Whitney', textFormat='star', loc='inside',
                                    verbose=0, lineOffsetToBox=2, lineHeight=0, stack=False, useFixedOffset=True)

    # context_stat_plot.legend(bbox_to_anchor=(1.05, 0.5), loc=3, borderaxespad=0., fontsize='small')
    plt.tight_layout()
    plt.savefig(output_dir + "Stat_synonymous_variant_frequency_in_%s.png" % virus, dpi=300)
    plt.close()

#

def main():
    input_dir = "/Users/odedkushnir/Projects/signatures/ADAR/SRA/"
    virus_path_dict = {"RV": "SRP006391_RV/", "PV": "SRP064468_PV/", "EnteroA": "ERP014415_Entero_A/"}
    output_dir = input_dir + "20190811_plots/"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)
    for virus in virus_path_dict:
        plots_for_srr(input_dir, output_dir, virus_path_dict[virus], virus)


if __name__ == "__main__":
        main()



