import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
from statannotations.Annotator import Annotator

def main():
    freq = pd.read_pickle("D:/My Drive/Studies/PhD/Projects/RV/RVB14/p12_1/data_filter.pkl")
    replica = 3
    table = freq.loc[freq["replica"] != replica]
    # table["Filter"] = np.where(table["passage"] == 2, True, np.where(table["passage"] == 5, True, False))
    # table = table.loc[table["Filter"] == True]
    transition_order = ["A>G", "U>C", "G>A", "C>U"]
    table["Mutation_Type"] = np.where(table["Mutation"] == "A>G", "transition", np.where(table["Mutation"] == "U>C", "transition",
                                                                                         np.where(table["Mutation"] == "G>A", "transition",
                                                                                                  np.where(table["Mutation"] == "C>U", "transition", "transversion"))))
    table = table.loc[table["Mutation_Type"] == "transition"]
    first_quantile = table["Frequency"].quantile(0.25)
    last_quantile = table["Frequency"].quantile(0.75)
    table["Filter"] = np.where(table["Frequency"] >= first_quantile, True, np.where(table["Frequency"] <= last_quantile, True, False))
    table = table.loc[table["Filter"] == True]

    fig1 = sns.catplot(y="Frequency", x="Mutation", data=table, hue="Type", kind="box", col="Mutation_Type",
                       hue_order=["Synonymous", "Premature Stop Codon", "Non-Synonymous"], order=transition_order)
    fig1.set(yscale='log', ylim=(10 ** -5, 10 ** -2))
    plt.axhline(y=5*10**-4, color='r', linestyle='--')
    # plt.show()
    plt.savefig("D:/My Drive/Studies/PhD/Projects/RV/RVB14/p12_1/transition_freq_Accu_quant.png", dpi=300)
    plt.close()

    fig3 = sns.catplot(y="Frequency", x="passage", data=table, hue="Type", kind="box",
                      hue_order=["Synonymous", "Premature Stop Codon"], order=range(0, 13, 1))
    fig3.set(yscale='log', ylim=(10 ** -5, 10 ** -2))
    fig3.set(xticklabels=["RNA\nControl", "", "2", "", "", "5", "", "", "8", "", "10", "", "12"])
    plt.axhline(y=5*10**-4, color='r', linestyle='--')
    plt.savefig("D:/My Drive/Studies/PhD/Projects/RV/RVB14/p12_1/GA_freq_Accu_boxplot_quant.png", dpi=300)
    plt.close()

    table_ga = table.loc[table["Mutation"] == "G>A"]
    table_ga = table.loc[table["passage"] != 0]
    table_ga_syn = table.loc[table["Type"] == "Synonymous"]
    table_ga_pmsc = table.loc[table["Type"] == "Premature Stop Codon"]
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(table_ga_syn['passage'],
                                                                        table_ga_syn[
                                                                            'Frequency'])
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(table_ga_pmsc['passage'],
                                                                        table_ga_pmsc[
                                                                            'Frequency'])
    fig2 = sns.lmplot(y="Frequency", x="passage", data=table_ga, hue="Type",
                       hue_order=["Synonymous", "Premature Stop Codon"])#, order=range(0, 13, 1))
    fig2.set(yscale='log', ylim=(10 ** -5, 10 ** -2))

    fig2.fig.subplots_adjust(wspace=.02)
    ax = fig2.axes[0, 0]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    label_line_1 = "log(y)=1e{0:.3g} X+{1:.3g} pval={2:.3g}".format(np.log10(slope1), intercept1, p_value1)
    label_line_2 = "log(y)=1e{0:.3g} X+{1:.3g} pval={2:.3g}".format(np.log10(slope2), intercept2, p_value2)
    L_labels[0].set_text(label_line_1)
    L_labels[1].set_text(label_line_2)
    plt.savefig("D:/My Drive/Studies/PhD/Projects/RV/RVB14/p12_1/GA_freq_Accu_lmplot_quant.png", dpi=300)
    plt.close()

    table = freq.loc[freq["replica"] != replica]
    # table["Filter"] = np.where(table["passage"] == 2, True, np.where(table["passage"] == 5, True, False))
    # table = table.loc[table["Filter"] == True]
    transition_order = ["A>G", "U>C", "G>A", "C>U"]
    table["Mutation_Type"] = np.where(table["Mutation"] == "A>G", "transition", np.where(table["Mutation"] == "U>C", "transition",
                                                                                         np.where(table["Mutation"] == "G>A", "transition",
                                                                                                  np.where(table["Mutation"] == "C>U", "transition", "transversion"))))
    table = table.loc[table["Mutation_Type"] == "transition"]
    table = table.loc[table["Mutation"] == "G>A"]
    table = table.loc[table["passage"] != 0]
    table_syn = table.loc[table["Type"] == "Synonymous"]
    table_pmsc = table.loc[table["Type"] == "Premature Stop Codon"]
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(table_syn['passage'],
                                                                        table_syn[
                                                                            'Frequency'])
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(table_pmsc['passage'],
                                                                        table_pmsc[
                                                                            'Frequency'])
    fig2 = sns.lmplot(y="Frequency", x="passage", data=table, hue="Type",
                       hue_order=["Synonymous", "Premature Stop Codon"])#, order=range(0, 13, 1))
    fig2.set(yscale='log', ylim=(10 ** -5, 10 ** -2))

    fig2.fig.subplots_adjust(wspace=.02)
    ax = fig2.axes[0, 0]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    label_line_1 = "log(y)=1e{0:.3g} X+{1:.3g} pval={2:.3g}".format(slope1, intercept1, p_value1)
    label_line_2 = "log(y)=1e{0:.3g} X+{1:.3g} pval={2:.3g}".format(slope2, intercept2, p_value2)
    L_labels[0].set_text(label_line_1)
    L_labels[1].set_text(label_line_2)
    plt.savefig("D:/My Drive/Studies/PhD/Projects/RV/RVB14/p12_1/GA_freq_Accu_lmplot.png", dpi=300)
    plt.close()
    #OdedKushnir


if __name__ == "__main__":
    main()