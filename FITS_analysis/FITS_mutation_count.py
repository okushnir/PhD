import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def main():
    data = pd.read_csv("/Users/odedkushnir/Library/CloudStorage/GoogleDrive-odedkushnir@gmail.com/My Drive/Studies/PhD/Thesis/FinalVer/AfterReview/mutation_rate_FITS_count.csv")
    plt.style.use('classic')
    sns.set_palette("Set2")
    virus_order = ["RVB14 #1", "RVB14 #2", "RVB14 #3", "CVB3", "OPV2", "PV1"]
    mutation_order = ["C>U", "G>A", "U>C", "A>G"]
    g2 = sns.catplot(x="Mutation", y="Count", data=data, hue="Virus", order=mutation_order, hue_order=virus_order, kind="bar", legend=False)
    g2.set(yscale="log", ylim=(10**0, 10**5))
    # plt.gcf().canvas.get_renderer()
    plt.savefig("/Users/odedkushnir/Library/CloudStorage/GoogleDrive-odedkushnir@gmail.com/My Drive/Studies/PhD/Thesis/FinalVer/AfterReview/mutation_rate_FITS_count.png", dpi=300)


if __name__ == "__main__":
    main()

