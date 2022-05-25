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
    first_quantile = table["Frequency"].quantile(0.25)
    last_quantile = table["Frequency"].quantile(0.75)
    table = table[((table["Frequency"] >= first_quantile) | (table["Frequency"] <= last_quantile))]
    # table["Filter"] = np.where(table["passage"] == 2, True, np.where(table["passage"] == 5, True, False))
    # table = table.loc[table["Filter"] == True]
    transition_order = ["A>G", "U>C", "G>A", "C>U"]
    table["Mutation_Type"] = np.where(table["Mutation"] == "A>G", "transition", np.where(table["Mutation"] == "U>C", "transition",
                                                                                         np.where(table["Mutation"] == "G>A", "transition",
                                                                                                  np.where(table["Mutation"] == "C>U", "transition", "transversion"))))


    table = table.loc[table["Mutation_Type"] == "transition"]
    fig1 = sns.catplot(y="Frequency", x="Mutation", data=table, hue="Type", kind="box", col="Mutation_Type",
                       hue_order=["Synonymous", "Premature Stop Codon", "Non-Synonymous"], order=transition_order)
    fig1.set(yscale='log', ylim=(10 ** -5, 10 ** -2))
    plt.axhline(y=5*10**-4, color='r', linestyle='--')
    # plt.show()
    plt.savefig("D:/My Drive/Studies/PhD/Projects/RV/RVB14/p12_1/transition_freq_Accu_quant.png", dpi=300)
    plt.close()

    fig2 = sns.catplot(y="Frequency", x="passage", data=table, hue="Type", kind="box",
                      hue_order=["Synonymous", "Premature Stop Codon"], order=range(0, 13, 1))
    fig2.set(yscale='log', ylim=(10 ** -5, 10 ** -2))
    fig2.set(xticklabels=["RNA\nControl", "", "2", "", "", "5", "", "", "8", "", "10", "", "12"])
    plt.axhline(y=5*10**-4, color='r', linestyle='--')
    plt.savefig("D:/My Drive/Studies/PhD/Projects/RV/RVB14/p12_1/transition_freq_Accu_boxplot_quant.png", dpi=300)
    plt.close()

    table_ga = table.loc[table["Mutation"] == "G>A"]
    fig3 = sns.catplot(y="Frequency", x="passage", data=table_ga, hue="Type", kind="box",
                      hue_order=["Synonymous", "Premature Stop Codon"], order=range(0, 13, 1))
    fig3.set(yscale='log', ylim=(10 ** -5, 10 ** -2))
    fig3.set(xticklabels=["RNA\nControl", "", "2", "", "", "5", "", "", "8", "", "10", "", "12"])
    plt.axhline(y=5*10**-4, color='r', linestyle='--')
    plt.savefig("D:/My Drive/Studies/PhD/Projects/RV/RVB14/p12_1/GA_freq_Accu_boxplot_quant.png", dpi=300)
    plt.close()

    table_ga = table.loc[table["passage"] != 0]
    table_ga["log10(Frequency)"] = np.log10(table_ga["Frequency"])
    table_ga_syn = table_ga.loc[table_ga["Type"] == "Synonymous"]
    table_ga_pmsc = table_ga.loc[table_ga["Type"] == "Premature Stop Codon"]
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(table_ga_syn['passage'], table_ga_syn["log10(Frequency)"])
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(table_ga_pmsc['passage'], table_ga_pmsc["log10(Frequency)"])
    fig4 = sns.lmplot(y="log10(Frequency)", x="passage", data=table_ga, hue="Type",
                      hue_order=["Synonymous", "Premature Stop Codon"])
    # fig4.set(yscale='log', ylim=(10 ** -5, 10 ** -2))

    fig4.fig.subplots_adjust(wspace=.02)
    ax = fig4.axes[0, 0]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    label_line_1 = "log(y)={0:.3g}x +({1:.3g}) pval={2:.3g}".format(slope1, intercept1, p_value1)
    label_line_2 = "log(y)={0:.3g}x +({1:.3g}) pval={2:.3g}".format(slope2, intercept2, p_value2)
    L_labels[0].set_text(label_line_1)
    L_labels[1].set_text(label_line_2)
    plt.savefig("D:/My Drive/Studies/PhD/Projects/RV/RVB14/p12_1/GA_freq_Accu_lmplot_quant.png", dpi=300)
    plt.close()

    #OdedKushnir


if __name__ == "__main__":
    main()