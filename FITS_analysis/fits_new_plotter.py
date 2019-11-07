
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.ticker import Locator
import math
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker
import numpy as np
from scipy import stats
import seaborn as sns

sns.set_style("ticks")


class MinorSymLogLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.
    """
    def __init__(self, linthresh):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically
        """
        self.linthresh = linthresh

    def __call__(self):
        'Return the locations of the ticks'
        majorlocs = self.axis.get_majorticklocs()

        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in range(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i-1]
            if abs(majorlocs[i-1] + majorstep/2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i-1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minorlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '%s type.' % type(self))


def post_data_mutation(input_dir):
    files = glob.glob(input_dir + "/posterior_mutation_*")
    columns = ["distance", "allele0_1", "Mutation", "label"]
    df = pd.DataFrame()
    for file in files:
        data = pd.read_csv(file, sep="\t")
        data["Mutation"] = file.split("_")[-1].split(".")[0]
        data["label"] = file.split("/")[-6]
        df = df.append(data)
    df = pd.DataFrame(df, columns=columns)
    return df

def post_data_fitness(input_dir):
    mutation_lst = glob.glob(input_dir + "/*")
    columns = ["pos", "inferred_w", "category", "levenes_p", "filename", "Mutation", "label"]
    df = pd.DataFrame()
    for mutation in mutation_lst:
        mutation = mutation.split("/")[-1]
        file = input_dir + "/" + str(mutation) + "/all.txt"
        data = pd.read_csv(file, sep="\t")
        data["Mutation"] = file.split("/")[-2]
        data["label"] = file.split("/")[-7]
        df = df.append(data)
    df = pd.DataFrame(df, columns=columns)
    return df

def syn_fitness(input_dir):
    files = glob.glob(input_dir + "/posterior_fitness_syn_*")
    columns = ["distance", "allele1", "Mutation", "label"]
    df = pd.DataFrame()
    for file in files:
        data = pd.read_csv(file, sep="\t")
        data["Mutation"] = file.split("_")[-1].split(".")[0]
        data["label"] = file.split("/")[-7]
        df = df.append(data)
    df = pd.DataFrame(df, columns=columns)
    return df

def x_round(x):
    return math.ceil(x*10)/10

def qqplot(x, y, **kwargs):
    _, xr = stats.probplot(x, fit=False)
    _, yr = stats.probplot(y, fit=False)
    plt.scatter(xr, yr, **kwargs)

def main():
    passages = "p0-p12"
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV"
    rv_mutation_data = post_data_mutation(input_dir + "/RVB14/fits/output/mutation/%s" % passages)
    cv_mutation_data = post_data_mutation(input_dir + "/CVB3/fits/output/mutation/%s" % passages)
    output_dir = input_dir + "/RVB14/fits/output/fitness/plots"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    all_data = pd.concat([rv_mutation_data, cv_mutation_data], sort=False)
    # print(all_data.to_string())
    all_data = all_data.rename(columns={"allele0_1": "Transition rate"})
    # all_data["Transition rate"] = np.log10(all_data["Transition rate"])
    # all_data["allele0_1"] = all_data["allele0_1"]*10**5

    plt.style.use('classic')
    g1 = sns.violinplot(x="Mutation", hue="label", y="Transition rate", data=all_data, split=True, cut=0,
                        order=["AG", "adar", "nonadar", "UC", "GA", "CU"])

    # g1.set_ylabel("Transition rate")
    g1.set_yscale("log")
    g1.set_yscale('symlog', linthreshy=1*10**-6)
    yaxis = plt.gca().yaxis
    yaxis.set_minor_locator(MinorSymLogLocator(1e-1))
    g1.set_title("Transition rate distribution in RVB14 and CVB3")
    # g1.set_yticks(ticks=[10**-5, 10**-6, 0], minor=True)
    # g1.set_ylim(10 ** -8, 10 ** -3)
    # g1.legend(bbox_to_anchor=(1.05, 0.5), loc="upper right", borderaxespad=0., fontsize='small')
    plt.tight_layout()
    # plt.show()
    plt.savefig(output_dir + "/20190901_posterior_dist.png", dpi=300)
    plt.close()

    mutation_lst = ["AG", "AG_adar", "AG_nonadar", "GA", "UC", "CU"]

    rv_fitness_data = post_data_fitness(input_dir + "/RVB14/fits/output/fitness/%s" % passages)
    rv_fitness_data = rv_fitness_data.rename(columns={"pos": "Pos"})
    rv_fitness_data = rv_fitness_data.rename(columns={"inferred_w": "Fitness"})
    rv_fitness_data["Protein"] = np.where(rv_fitness_data["Pos"] <= 629, "5'UTR",
                                             np.where(rv_fitness_data["Pos"] <= 835, "VP4",
                                                      np.where(rv_fitness_data["Pos"] <= 1621, "VP2",
                                                               np.where(rv_fitness_data["Pos"] <= 2329, "VP3",
                                                                        np.where(rv_fitness_data["Pos"] <= 3196,
                                                                                 "VP1",
                                                                                 np.where(
                                                                                     rv_fitness_data["Pos"] <= 3634,
                                                                                     "2A",
                                                                                     np.where(rv_fitness_data[
                                                                                                  "Pos"] <= 3925, "2B",
                                                                                              np.where(
                                                                                                  rv_fitness_data[
                                                                                                      "Pos"] <= 4915,
                                                                                                  "2C",
                                                                                                  np.where(
                                                                                                      rv_fitness_data[
                                                                                                          "Pos"] <= 5170,
                                                                                                      "3A",
                                                                                                      np.where(
                                                                                                          rv_fitness_data[
                                                                                                              "Pos"] <= 5239,
                                                                                                          "3B",
                                                                                                          np.where(
                                                                                                              rv_fitness_data[
                                                                                                                  "Pos"] <= 5785,
                                                                                                              "3C",
                                                                                                              np.where(
                                                                                                                  rv_fitness_data[
                                                                                                                      "Pos"] <= 7168,
                                                                                                                  "3D",
                                                                                                                  "3'UTR"))))))))))))
    rv_fitness_data = rv_fitness_data[rv_fitness_data["Pos"] > 3635] # just 2B
    rv_fitness_data["Filter"] = rv_fitness_data["Protein"] != "3'UTR"
    rv_fitness_data = rv_fitness_data[rv_fitness_data["Filter"] == True]

    cv_fitness_data = post_data_fitness(input_dir + "/CVB3/fits/output/fitness/%s" % passages)
    cv_fitness_data = cv_fitness_data.rename(columns={"pos": "Pos"})
    cv_fitness_data = cv_fitness_data.rename(columns={"inferred_w": "Fitness"})
    cv_fitness_data["Protein"] = np.where(cv_fitness_data["Pos"] <= 739, "5'UTR",
                                             np.where(cv_fitness_data["Pos"] <= 948, "VP4",
                                                      np.where(cv_fitness_data["Pos"] <= 1737, "VP2",
                                                               np.where(cv_fitness_data["Pos"] <= 2451, "VP3",
                                                                        np.where(cv_fitness_data["Pos"] <= 3303,
                                                                                 "VP1",
                                                                                 np.where(
                                                                                     cv_fitness_data["Pos"] <= 3744,
                                                                                     "2A",
                                                                                     np.where(cv_fitness_data[
                                                                                                  "Pos"] <= 4041, "2B",
                                                                                              np.where(
                                                                                                  cv_fitness_data[
                                                                                                      "Pos"] <= 5028,
                                                                                                  "2C",
                                                                                                  np.where(
                                                                                                      cv_fitness_data[
                                                                                                          "Pos"] <= 5295,
                                                                                                      "3A",
                                                                                                      np.where(
                                                                                                          cv_fitness_data[
                                                                                                              "Pos"] <= 5361,
                                                                                                          "3B",
                                                                                                          np.where(
                                                                                                              cv_fitness_data[
                                                                                                                  "Pos"] <= 5910,
                                                                                                              "3C",
                                                                                                              np.where(
                                                                                                                  cv_fitness_data[
                                                                                                                      "Pos"] <= 7299,
                                                                                                                  "3D",
                                                                                                                  "3'UTR"))))))))))))
    cv_fitness_data = cv_fitness_data[cv_fitness_data["Pos"] > 3745] # just 2B
    cv_fitness_data["Pos"] = cv_fitness_data["Pos"]-110 #align CVB3 and RVB14

    fitness_data = pd.concat([rv_fitness_data, cv_fitness_data], sort=False)
    # fitness_data = cv_fitness_data.rename(columns={"inferred_w": "Fitness"})
    # print(all_data.to_string())
    # rv_fitness_data = rv_fitness_data.rename(columns={"inferred_w": "Fitness"})
    # rv_fitness_data["inferred_w"] = round(rv_fitness_data["inferred_w"], 2)
    # rv_fitness_data["Fitness"] = rv_fitness_data["inferred_w"].apply(lambda x: x_round(x))


    # print(rv_fitness_data.to_string())
    mutation_order = ["AG", "AG_adar", "AG_nonadar", "UC", "GA", "CU"]


    #Plots
    bins = np.arange(0, 2, 0.1)
    g2 = sns.FacetGrid(data=fitness_data, col="Mutation", col_order=mutation_order
                       , col_wrap=3,  legend_out=True, hue="label", palette="Set1")
    kws = dict(linewidth=.5, edgecolor="w", alpha=0.5)
    g2 = (g2.map(plt.hist, "Fitness", **kws, bins=bins).set_titles("{col_name}").add_legend())

    g2.fig.suptitle("DFE of RVB14 and CVB3 (%s)" % passages, fontsize=16, weight='bold', ha='center', va='bottom')
    g2.savefig(output_dir + "/20190908_DFE_RVB14_CVB3.png", dpi=300)
    plt.close()

    # pal = dict(CVB3="seagreen", RVB14="gray")
    g3 = sns.FacetGrid(fitness_data, row="Mutation", hue="label", row_order=mutation_order, col="label")#, col_wrap=3)
    g3 = (g3.map(plt.scatter, "Pos", "Fitness", **kws).add_legend())
    # g3.fig.suptitle(" RVB14 and CVB3 (%s)" % passages, fontsize=16, weight='bold', ha='center', va='bottom')
    g3.savefig(output_dir + "/20190909_Mutation_Pos_RVB14_CVB3.png", dpi=300)
    plt.close()

    g4 = sns.catplot(x="Mutation", y="Fitness", data=fitness_data, col="Protein", hue="label",
                     col_wrap=3, kind="box", order=mutation_order)
    # g4 = (g4.map(plt.boxplot, ).add_legend())
    # g4.fig.suptitle(" RVB14 and CVB3 (%s)" % passages, fontsize=16, weight='bold', ha='center', va='bottom')
    g4.savefig(output_dir + "/20190909_Mutation_Protein_RVB14_CVB3_box.png", dpi=300)
    plt.close()

    g5 = sns.FacetGrid(fitness_data, col="Mutation", hue="label", col_order=mutation_order, col_wrap=3)
    g5 = (g5.map(qqplot, "Fitness", "Pos", **kws).add_legend())
    g5.fig.suptitle("qqplot of RVB14 and CVB3 (%s)" % passages, fontsize=16, weight='bold', ha='center', va='bottom')
    g5.savefig(output_dir + "/20190909_qqplot_RVB14_CVB3.png", dpi=300)
    plt.close()

    rvb14 = fitness_data.loc[fitness_data.label == "RVB14"]
    cvb3 = fitness_data.loc[fitness_data.label == "CVB3"]
    ax = sns.kdeplot(rvb14.Fitness, label="RVB14")
    ax = sns.kdeplot(cvb3.Fitness, label="CVB3")
    plt.legend(title='Virus', loc='upper right')
    plt.suptitle("Fitness KDE in RVB14 and CVB3 (%s)" % passages, fontsize=16, weight='bold', ha='center')
    # plt.tight_layout()
    plt.savefig(output_dir + "/20190909_KDE_RVB14_CVB3.png", dpi=300)
    plt.close()

    rvb14 = fitness_data.loc[fitness_data.label == "RVB14"]
    rvb14_AG = rvb14.loc[rvb14.Mutation == "AG"]
    rvb14_AG_adar = rvb14.loc[rvb14.Mutation == "AG_adar"]
    rvb14_AG_nonadar = rvb14.loc[rvb14.Mutation == "AG_nonadar"]
    rvb14_UC = rvb14.loc[rvb14.Mutation == "UC"]
    rvb14_GA = rvb14.loc[rvb14.Mutation == "GA"]
    rvb14_CU = rvb14.loc[rvb14.Mutation == "CU"]

    cvb3 = fitness_data.loc[fitness_data.label == "CVB3"]
    cvb3_AG = cvb3.loc[cvb3.Mutation == "AG"]
    cvb3_AG_adar = cvb3.loc[cvb3.Mutation == "AG_adar"]
    cvb3_AG_nonadar = cvb3.loc[cvb3.Mutation == "AG_nonadar"]
    cvb3_UC = cvb3.loc[cvb3.Mutation == "UC"]
    cvb3_GA = cvb3.loc[cvb3.Mutation == "GA"]
    cvb3_CU = cvb3.loc[cvb3.Mutation == "CU"]

    #["AG", "AG_adar", "AG_nonadar", "GA", "UC", "CU"]
    fig, axes = plt.subplots(2, 3, sharey=True, sharex=True)
    sns.kdeplot(rvb14_AG.Fitness, label="RVB14", ax=axes[0, 0], legend=False)
    sns.kdeplot(cvb3_AG.Fitness, label="CVB3", ax=axes[0, 0], legend=False)

    sns.kdeplot(rvb14_AG_adar.Fitness, label="RVB14", ax=axes[0, 1], legend=False)
    sns.kdeplot(cvb3_AG_adar.Fitness, label="CVB3", ax=axes[0, 1], legend=False)

    sns.kdeplot(rvb14_AG_nonadar.Fitness, label="RVB14", ax=axes[0, 2], legend=False)
    sns.kdeplot(cvb3_AG_nonadar.Fitness, label="CVB3", ax=axes[0, 2], legend=False)

    sns.kdeplot(rvb14_GA.Fitness, label="RVB14", ax=axes[1, 0], legend=False)
    sns.kdeplot(cvb3_GA.Fitness, label="CVB3", ax=axes[1, 0], legend=False)

    sns.kdeplot(rvb14_UC.Fitness, label="RVB14", ax=axes[1, 1], legend=False)
    sns.kdeplot(cvb3_UC.Fitness, label="CVB3", ax=axes[1, 1], legend=False)

    sns.kdeplot(rvb14_CU.Fitness, label="RVB14", ax=axes[1, 2], legend=False)
    sns.kdeplot(cvb3_CU.Fitness, label="CVB3", ax=axes[1, 2], legend=False)

    plt.legend(title='Virus', loc="upper right", bbox_to_anchor=(1.05, 1.05), borderaxespad=0.)
    # plt.suptitle("Fitness KDE in RVB14 and CVB3\n(%s)" % passages, fontsize=16, weight='bold', ha='center')
    plt.tight_layout()
    plt.savefig(output_dir + "/20190909_KDE_Mutations.png", dpi=300)
    plt.close()


    # rv_fitness_syn_data = syn_fitness(input_dir + "/RVB14/fits/output/fitness/%s/synonymous" % passages)
    # cv_fitness_syn_data = syn_fitness(input_dir + "/CVB3/fits/output/fitness/%s/synonymous" % passages)
    #
    # fitness_data = pd.concat([rv_fitness_syn_data, cv_fitness_syn_data], sort=False)
    # print(all_data.to_string())
    # fitness_data = fitness_data.rename(columns={"allele1": "Fitness"})
    #
    # g3 = sns.violinplot(x="Mutation", hue="label", y="Fitness", data=fitness_data, split=True, cut=0,
    #                     order=["AG", "adar", "nonadar", "UC", "GA", "CU"])
    # g3.set_title("Fitness distribution of synonymous mutations in RVB14 and CVB3\n(%s)" % passages)
    # plt.tight_layout()
    # # plt.show()
    # plt.savefig(output_dir + "/20190901_posterior_dist_fitness_synon.png", dpi=300)
    # plt.close()


if __name__ == "__main__":
    main()
