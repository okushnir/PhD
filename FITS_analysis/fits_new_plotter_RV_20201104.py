
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import Locator
import math
from scipy import stats
import seaborn as sns
from statannot import add_stat_annotation
from AccuNGS_analysis.add_Protein_to_pd_df import *
from statannotations.Annotator import Annotator
import datetime

sns.set(font_scale=1.2)
sns.set_style("ticks")
sns.despine()


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


def suplabel(axis,label,label_prop=None,
             labelpad=5,
             ha='center',va='center'):
    ''' Add super ylabel or xlabel to the figure
    Similar to matplotlib.suptitle
    axis       - string: "x" or "y"
    label      - string
    label_prop - keyword dictionary for Text
    labelpad   - padding from the axis (default: 5)
    ha         - horizontal alignment (default: "center")
    va         - vertical alignment (default: "center")
    '''
    fig = pylab.gcf()
    xmin = []
    ymin = []
    for ax in fig.axes:
        xmin.append(ax.get_position().xmin)
        ymin.append(ax.get_position().ymin)
    xmin,ymin = min(xmin),min(ymin)
    dpi = fig.dpi
    if axis.lower() == "y":
        rotation=90.
        x = xmin-float(labelpad)/dpi
        y = 0.5
    elif axis.lower() == 'x':
        rotation = 0.
        x = 0.5
        y = ymin - float(labelpad)/dpi
    else:
        raise Exception("Unexpected axis: x or y")
    if label_prop is None:
        label_prop = dict()
    pylab.text(x,y,label,rotation=rotation,
               transform=fig.transFigure,
               ha=ha,va=va,
               **label_prop)

def post_data_mutation(input_dir):
    mutation_lst = glob.glob(input_dir + "/*")
    columns = ["pos", "inferred_mu", "levenes_p", "filename", "Mutation"]
    df = pd.DataFrame()
    for mutation in mutation_lst:
        mutation = mutation.split("/")[-1]
        file = input_dir + "/" + str(mutation) + "/all.txt"
        data = pd.read_csv(file, sep="\t")
        data["Mutation"] = file.split("/")[-2]
        # data["label"] = file.split("/")[-7]
        df = df.append(data)
    df = pd.DataFrame(df, columns=columns)
    return df

# def post_data_fitness(input_dir):
#     mutation_lst = glob.glob(input_dir + "/*")
#     columns = ["pos", "inferred_w", "category", "levenes_p", "filename", "Mutation"]
#     df = pd.DataFrame()
#     for mutation in mutation_lst:
#         mutation = mutation.split("/")[-1]
#         file = input_dir + "/" + str(mutation) + "/all.txt"
#         data = pd.read_csv(file, sep="\t")
#         data["Mutation"] = file.split("/")[-2]
#         # data["label"] = file.split("/")[-7]
#         df = df.append(data)
#     df = pd.DataFrame(df, columns=columns)
#     return df

"""All At once"""
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

def syn_fitness(input_dir):
    files = glob.glob(input_dir + "/posterior_fitness_syn_*")
    columns = ["distance", "allele1", "Mutation"]
    df = pd.DataFrame()
    for file in files:
        data = pd.read_csv(file, sep="\t")
        data["Mutation"] = file.split("_")[-1].split(".")[0]
        # data["label"] = file.split("/")[-7]
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
    date = datetime.date.today().strftime("%Y%m%d")
    passages = "p2-p12"
    opv_passages = "p1-p7"
    pv_passages = "p3-p8"
    input_dir = "/Users/odedkushnir/Projects/fitness"
    rv_replica1_mutation_data = post_data_mutation("/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/"
                                                   "20201008RV-202329127/merged/passages/fits_all_pos_at_once_sampling/"
                                                   "replica1_syn/output/mutation/%s" % passages)
    rv_replica1_mutation_data["Virus"] = "RVB14 #1"
    rv_replica2_mutation_data = post_data_mutation("/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/"
                                                   "20201008RV-202329127/merged/passages/fits_all_pos_at_once_sampling/"
                                                   "replica2_syn/output/mutation/%s" % passages)
    rv_replica2_mutation_data["Virus"] = "RVB14 #2"
    rv_replica3_mutation_data = post_data_mutation("/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/"
                                                   "20201008RV-202329127/merged/passages/fits_all_pos_at_once_sampling/"
                                                   "replica3_syn/output/mutation/%s" % passages)
    rv_replica3_mutation_data["Virus"] = "RVB14 #3"
    cv_mutation_data = post_data_mutation("/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/190627_RV_CV"
                                                      "/merged/CVB3/Rank0_data_mutation/fits/output/mutation/%s" % passages)
    cv_mutation_data["Virus"] = "CVB3"
    opv_mutataion_data = post_data_mutation("/Volumes/STERNADILABHOME$/volume3/okushnir/CirSeq/OPV/fits/output/mutation/"
                                            "all_positions_p1-p7")
    opv_mutataion_data["Virus"] = "OPV2"

    pv_mutataion_data = post_data_mutation("/Volumes/STERNADILABHOME$/volume3/okushnir/CirSeq/Mahoney/fits/output/"
                                           "mutation/p3-p8")
    pv_mutataion_data["Virus"] = "PV1"



    output_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/" \
                             "%s_fits_syn_plots" % date
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    all_data = pd.concat([rv_replica1_mutation_data, rv_replica2_mutation_data, rv_replica3_mutation_data,
                          cv_mutation_data, opv_mutataion_data, pv_mutataion_data], sort=False)
    all_data = all_data.rename(columns={"allele0_1": "Transition rate"})
    all_data["Transition rate"] = all_data["Transition rate"].astype(float)
    # print(all_data.to_string())
    # all_data = all_data.rename(columns={"inferred_mu": "Mutation rate"})
    # # print(all_data["Mutation rate"].dtype)
    # all_data["Mutation rate"] = all_data["Mutation rate"].map(lambda x: str(x).lstrip('*'))
    # all_data["Mutation rate"] = pd.to_numeric(all_data["Mutation rate"], errors='coerce')#.astype(float)
    # # print(all_data["Mutation rate"].dtype)
    # all_data["Mutation"] = all_data["Mutation"].apply(lambda x: x[0]+">"+x[1:]if len(x)<=2 else x)# if len(x)==2 else x[0]+">"+x[1:])
    # all_data["Mutation"] = all_data["Mutation"].apply(lambda x: x.split("_")[0] + "\n" + x.split("_")[-1] + "-like" if len(x)>3 else x)
    all_data["Mutation"] = np.where(all_data["Mutation"] == "nonadar", "A>G\nNon-ADAR-like", all_data["Mutation"])
    all_data["Mutation"] = np.where(all_data["Mutation"] == "adar", "A>G\nADAR-like", all_data["Mutation"])
    all_data["Mutation"] = np.where(all_data["Mutation"] == "AG", "A>G", all_data["Mutation"])
    all_data["Mutation"] = np.where(all_data["Mutation"] == "UC", "U>C", all_data["Mutation"])
    all_data["Mutation"] = np.where(all_data["Mutation"] == "GA", "G>A", all_data["Mutation"])
    all_data["Mutation"] = np.where(all_data["Mutation"] == "CU", "C>U", all_data["Mutation"])
    # all_data = all_data[(all_data["pos"] >= 5785) & (all_data["pos"] <= 7212)]



    # q1 = all_data["Transition rate"].quantile(0.25)
    # q3 = all_data["Transition rate"].quantile(0.75)
    # all_data = all_data[all_data["Transition rate"] > q1]
    # all_data = all_data[all_data["Transition rate"] < q3]

    all_data = all_data[all_data["Mutation"] != "A>G\nADAR-like"]
    all_data = all_data[all_data["Mutation"] != "A>G\nNon-ADAR-like"]
    print(all_data.shape[0])
    all_data.to_csv("/Users/odedkushnir/PhD_Projects/fitness/all_data.csv")

    #Plots - local
    all_data = pd.read_csv("/Users/odedkushnir/PhD_Projects/fitness/all_data.csv")
    output_dir = "/Users/odedkushnir/PhD_Projects/fitness/{0}_fits_syn_plots".format(date)
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    plt.style.use('classic')

    sns.set_palette("Set2")
    mutation_order = ["C>U", "G>A", "U>C", "A>G"]
    virus_order = ["RVB14 #1", "RVB14 #2", "RVB14 #3", "CVB3", "OPV2", "PV1"]
    g1 = sns.boxenplot(x="Mutation", y="Transition rate", data=all_data, order=mutation_order, hue="Virus",
                       hue_order=virus_order)
    g1.set_yscale("log")
    """[((cat1, hue1), (cat2, hue2)), ((cat3, hue3), (cat4, hue4))]"""
    pairs = [(("A>G", "RVB14 #1"), ("C>U", "RVB14 #1")), (("A>G", "RVB14 #1"), ("G>A","RVB14 #1")), (("A>G", "RVB14 #1"), ("U>C","RVB14 #1")),
             (("A>G", "RVB14 #2"), ("C>U", "RVB14 #2")), (("A>G", "RVB14 #2"), ("G>A","RVB14 #2")), (("A>G", "RVB14 #2"), ("U>C","RVB14 #2")),
             (("C>U", "RVB14 #3"), ("A>G", "RVB14 #3")), (("A>G", "RVB14 #3"), ("G>A","RVB14 #3")), (("A>G", "RVB14 #3"), ("U>C","RVB14 #3")),
             (("A>G", "CVB3"), ("C>U", "CVB3")), (("A>G", "CVB3"), ("G>A", "CVB3")), (("A>G", "CVB3"), ("U>C", "CVB3")),
             (("A>G", "OPV2"), ("C>U", "OPV2")), (("A>G", "OPV2"), ("G>A", "OPV2")), (("A>G", "OPV2"), ("U>C", "OPV2")),
             (("A>G", "PV1"), ("C>U", "PV1")), (("A>G", "PV1"), ("G>A", "PV1")), (("A>G", "PV1"), ("U>C", "PV1"))]

    # [(("A>G", "RVB14 #1"), ("A>G", "CVB3")), (("A>G", "RVB14"), ("A>G","OPV")),
    #  (("U>C", "RVB14"), ("U>C", "CVB3")), (("U>C", "RVB14"),("U>C", "OPV")),
    #  (("G>A", "RVB14"), ("G>A", "CVB3")), (("G>A", "RVB14"),("G>A", "OPV")),
    #  (("C>U", "RVB14"), ("C>U", "CVB3")), (("C>U", "RVB14"),("C>U", "OPV"))]

    annotator = Annotator(g1, pairs, x="Mutation", y="Transition rate", data=all_data, order=mutation_order,
                          hue="Virus", hue_order=virus_order)
    annotator.configure(test='Mann-Whitney', text_format='star', loc='outside', comparisons_correction="Bonferroni")
    annotator.apply_and_annotate()
    # g1.set_xticklabels(labels=mutation_order, fontsize=8)
    # add_stat_annotation(g1, data=all_data, x="Mutation", y="Mutation rate", hue="Virus", order=mutation_order,
    #                     boxPairList=[(("A>G\nall", "RVB14"), ("A>G\nall", "CVB3")), (("A>G\nall", "RVB14"), ("A>G\nall",
    #                                                                                                          "OPV")),
    #                                  (("A>G\nadar-like", "RVB14"), ("A>G\nadar-like", "CVB3")),
    #                                  (("A>G\nadar-like", "RVB14"),("A>G\nadar-like", "OPV")),
    #                                  (("A>G\nnonadar-like", "RVB14"), ("A>G\nnonadar-like", "CVB3")),
    #                                  (("A>G\nnonadar-like", "RVB14"),("A>G\nnonadar-like", "OPV")),
    #                                  (("U>C", "RVB14"), ("U>C", "CVB3")), (("U>C", "RVB14"),("U>C", "OPV")),
    #                                  (("G>A", "RVB14"), ("G>A", "CVB3")), (("G>A", "RVB14"),("G>A", "OPV")),
    #                                  (("C>U", "RVB14"), ("C>U", "CVB3")), (("C>U", "RVB14"),("C>U", "OPV"))],
    #                     test='Mann-Whitney', textFormat='star', loc='inside', verbose=2)
    # add_stat_annotation(g1, data=all_data, x="Mutation", y="Mutation rate", hue="Virus", order=mutation_order,
    #                     boxPairList=[(("A>G\nall", "RVB14"), ("A>G\nadar-like", "RVB14")), (("A>G\nall", "RVB14"), ("A>G\nnonadar-like", "RVB14")),
    #                                  (("A>G\nall", "RVB14"), ("U>C", "RVB14")),
    #                                  (("A>G\nall", "RVB14"),("G>A", "RVB14")), (("A>G\nall", "RVB14"), ("C>U", "RVB14"))],
    #                     test='Mann-Whitney', textFormat='star', loc='inside', verbose=2)

    g1.set(xlabel="Type of mutation")
    g1.set(ylabel="Mutation rate inferred")

    # g1.set_yscale('symlog', linthreshy=1*10**-6)
    # yaxis = plt.gca().yaxis
    # yaxis.set_minor_locator(MinorSymLogLocator(1e-1))
    # g1.set_title("Mutation rate distribution")
    # g1.set_yticks(ticks=[10**-5, 10**-6, 0], minor=True)
    g1.set_ylim(10 ** -10, 10 ** -1)
    g1.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), borderaxespad=0., fontsize=7)
    # _, xlabels = plt.xticks()
    # g1.set_xticklabels(xlabels, size=7)
    # sns.set(font_scale=0.6)
    plt.tight_layout()
    # plt.show()
    plt.savefig(output_dir + "/%s_mutation_rate.png" % date, dpi=600)
    # plt.savefig(
    #     "/Users/odedkushnir/Google Drive/Studies/PhD/MyArticle-Signatures/plots/20191006_Mutation_rate.png",
    #     dpi=300)

    plt.close()

if __name__ == "__main__":
    main()
