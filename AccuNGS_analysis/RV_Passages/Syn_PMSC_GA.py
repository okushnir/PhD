import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from statannotations.Annotator import Annotator

def main():
    freq = pd.read_pickle("D:/My Drive/Studies/PhD/Projects/RV/RVB14/p12_1/data_filter.pkl")
    replica = 3
    table = freq.loc[freq["passage"] == 12]
    table = table.loc[table["replica"] != replica]
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
    plt.savefig("D:/My Drive/Studies/PhD/Projects/RV/RVB14/p12_1/transition_freq_p12_Accu.png", dpi=300)
    plt.close()
    #OdedKushnir


if __name__ == "__main__":
    main()