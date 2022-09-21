import matplotlib.pyplot as plt

import seaborn as sns
import pandas as pd


def main():
    data = pd.read_csv("/Users/odedkushnir/Library/CloudStorage/GoogleDrive-odedkushnir@gmail.com/My Drive/Studies/PhD/Thesis/FinalVer/AfterReview/mutation_rate_FITS_count.csv")
    plt.style.use('classic')
    sns.set_palette("Set2")
    virus_order = ["RVB14 #1", "RVB14 #2", "RVB14 #3", "CVB3", "OPV2", "PV1"]
    g2 = sns.barplot(x="Mutation", y="Count", data=data, hue="Virus", hue_order=virus_order)
    # plt.gcf().canvas.get_renderer()
    plt.show()
    #.savefig("/Users/odedkushnir/Library/CloudStorage/GoogleDrive-odedkushnir@gmail.com/My Drive/Studies/PhD/Thesis/FinalVer/AfterReview/mutation_rate_FITS_count.png", dpi=300)


if __name__ == "__main__":
    main()

