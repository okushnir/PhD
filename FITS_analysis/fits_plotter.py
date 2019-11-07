
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.ticker import Locator
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
    files = glob.glob(input_dir + "/posterior_fitness_all_*")
    columns = ["distance", "allele1", "Mutation", "label"]
    df = pd.DataFrame()
    for file in files:
        data = pd.read_csv(file, sep="\t")
        data["Mutation"] = file.split("_")[-1].split(".")[0]
        data["label"] = file.split("/")[-6]
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

def main():
    passages = "p0-p12"
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV"
    rv_mutation_data = post_data_mutation(input_dir + "/RVB14/fits/output/mutation/%s" % passages)
    cv_mutation_data = post_data_mutation(input_dir + "/CVB3/fits/output/mutation/%s" % passages)
    output_dir = input_dir + "/RVB14/fits/output/fitness/%s/plots" % passages
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    all_data = pd.concat([rv_mutation_data, cv_mutation_data], sort=False)
    print(all_data.to_string())
    all_data = all_data.rename(columns={"allele0_1": "Transition rate"})
    # all_data["Transition rate"] = np.log10(all_data["Transition rate"])
    # all_data["allele0_1"] = all_data["allele0_1"]*10**5

    g1 = sns.violinplot(x="Mutation", hue="label", y="Transition rate", data=all_data, split=True, cut=0,
                        order=["AG", "adar", "nonadar", "UC", "GA", "CU"])

    # g1.set_ylabel("Transition rate")
    g1.set_yscale("log")
    # g1.set_yscale('symlog', linthreshy=5*10**-9)
    # yaxis = plt.gca().yaxis
    # yaxis.set_minor_locator(MinorSymLogLocator(1e-1))
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
    cv_fitness_data = post_data_fitness(input_dir + "/CVB3/fits/output/fitness/%s" % passages)

    fitness_data = pd.concat([rv_fitness_data, cv_fitness_data], sort=False)
    print(all_data.to_string())
    fitness_data = fitness_data.rename(columns={"allele1": "Fitness"})

    g2 = sns.violinplot(x="Mutation", hue="label", y="Fitness", data=fitness_data, split=True, cut=0,
                        order=["AG", "adar", "nonadar", "UC", "GA", "CU"])
    g2.set_title("Fitness distribution in RVB14 and CVB3 (%s)" % passages)
    plt.tight_layout()
    # plt.show()
    plt.savefig(output_dir + "/20190901_posterior_dist_fitness.png", dpi=300)
    plt.close()

    rv_fitness_syn_data = syn_fitness(input_dir + "/RVB14/fits/output/fitness/%s/synonymous" % passages)
    cv_fitness_syn_data = syn_fitness(input_dir + "/CVB3/fits/output/fitness/%s/synonymous" % passages)

    fitness_data = pd.concat([rv_fitness_syn_data, cv_fitness_syn_data], sort=False)
    print(all_data.to_string())
    fitness_data = fitness_data.rename(columns={"allele1": "Fitness"})

    g3 = sns.violinplot(x="Mutation", hue="label", y="Fitness", data=fitness_data, split=True, cut=0,
                        order=["AG", "adar", "nonadar", "UC", "GA", "CU"])
    g3.set_title("Fitness distribution of synonymous mutations in RVB14 and CVB3\n(%s)" % passages)
    plt.tight_layout()
    # plt.show()
    plt.savefig(output_dir + "/20190901_posterior_dist_fitness_synon.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
