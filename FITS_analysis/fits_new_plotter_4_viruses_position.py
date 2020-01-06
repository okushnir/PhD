
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import Locator
import math
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker
import numpy as np
from scipy import stats
import seaborn as sns
from statannot import add_stat_annotation
from AccuNGS_analysis.add_Protein_to_pd_df import *

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
    columns = ["pos", "inferred_mu", "levenes_p", "filename", "Mutation", "label"]
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
    date = "20191009"
    passages = "p0-p12"
    opv_passages = "p1-p7"
    pv_passages = "p3-p8"
    p5_p12 = "p5-p12"
    input_dir = "/Users/odedkushnir/Projects/fitness"
    rv_mutation_data = post_data_mutation(input_dir + "/AccuNGS/190627_RV_CV/RVB14/fits/output/mutation/%s_all_pass" % passages)
    rv_p0_p5_p12_mutation_data = post_data_mutation(input_dir + "/AccuNGS/190627_RV_CV/RVB14/fits/output/mutation/p0-p12_without_p2")
    rv_p0_p5_p12_mutation_data["label"] = "RVB14_p0-p12\nwithout_p2"
    rv_p5_p12_mutation_data = post_data_mutation(input_dir + "/AccuNGS/190627_RV_CV/RVB14/fits/output/mutation/%s" % p5_p12)
    rv_p5_p12_mutation_data["label"] = "RVB14_p5-p12"
    cv_mutation_data = post_data_mutation(input_dir + "/AccuNGS/190627_RV_CV/CVB3/fits/output/mutation/%s_all_pass" % passages)
    cv_p0_p5_p12_mutation_data = post_data_mutation(input_dir + "/AccuNGS/190627_RV_CV/CVB3/fits/output/mutation/p0-p12_without_p2")
    cv_p0_p5_p12_mutation_data["label"] = "CVB3_p0-p12\nwithout_p2"
    cv_p5_p12_mutation_data = post_data_mutation(input_dir + "/AccuNGS/190627_RV_CV/CVB3/fits/output/mutation/%s" % p5_p12)
    cv_p5_p12_mutation_data["label"] = "CVB3_p5-p12"
    opv_mutataion_data = post_data_mutation(input_dir + "/CirSeq/PV/OPV/fits/output/mutation/%s" % opv_passages)
    pv_mutataion_data = post_data_mutation(input_dir + "/CirSeq/PV/Mahoney/fits/output/mutation/%s" % pv_passages)
    output_dir = input_dir + "/CirSeq/PV/Mahoney/fits/output/%s_syn_plots" % date
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    all_data = pd.concat([rv_mutation_data, rv_p0_p5_p12_mutation_data, rv_p5_p12_mutation_data, cv_mutation_data,
                          cv_p0_p5_p12_mutation_data, cv_p5_p12_mutation_data, opv_mutataion_data, pv_mutataion_data],
                         sort=False)
    # print(all_data.to_string())
    all_data = all_data.rename(columns={"inferred_mu": "Mutation rate"})
    print(all_data["Mutation rate"].dtype)
    all_data["Mutation rate"] = all_data["Mutation rate"].map(lambda x: str(x).lstrip('*'))
    all_data["Mutation rate"] = pd.to_numeric(all_data["Mutation rate"], errors='coerce')#.astype(float)
    print(all_data["Mutation rate"].dtype)
    all_data = all_data.rename(columns={"label": "Virus"})
    all_data["Mutation"] = all_data["Mutation"].apply(lambda x: x[0]+">"+x[1:]if len(x)<=2 else x)# if len(x)==2 else x[0]+">"+x[1:])
    # all_data["Mutation"] = all_data["Mutation"].apply(lambda x: x.split("_")[0] + "\n" + x.split("_")[-1] + "-like" if len(x)>3 else x)
    all_data["Mutation"] = np.where(all_data["Mutation"] == "AG_nonadar", "A>G\nnonadar-like", all_data["Mutation"])
    all_data["Mutation"] = np.where(all_data["Mutation"] == "AG_adar", "A>G\nadar-like", all_data["Mutation"])
    all_data["Mutation"] = np.where(all_data["Mutation"] == "A>G", "A>G\nall", all_data["Mutation"])
    # all_data["Mutation rate"] = np.log10(all_data["Mutation rate"])

    mutation_order = ["A>G\nall", "A>G\nadar-like", "A>G\nnonadar-like", "U>C", "G>A", "C>U"]

    #Plots
    plt.style.use('classic')

    sns.set_palette("Set2")

    g1 = sns.boxenplot(x="Mutation", y=all_data["Mutation rate"], hue="Virus", data=all_data, order=mutation_order)
    g1.set_yscale("log")
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

    # g1.set(ylabel="Muataion rate (log10)")

    # g1.set_yscale('symlog', linthreshy=1*10**-6)
    # yaxis = plt.gca().yaxis
    # yaxis.set_minor_locator(MinorSymLogLocator(1e-1))
    # g1.set_title("Mutation rate distribution")
    # g1.set_yticks(ticks=[10**-5, 10**-6, 0], minor=True)
    # g1.set_ylim(10 ** -8, 10 ** -3)
    g1.legend(bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0.)
    sns.set(font_scale=0.8)
    plt.tight_layout()
    # plt.show()
    plt.savefig(output_dir + "/%s_mutation_rate.png" % date, dpi=300)
    # plt.savefig(
    #     "/Users/odedkushnir/Google Drive/Studies/PhD/MyArticle-Signatures/plots/20191006_Mutation_rate.png",
    #     dpi=300)

    plt.close()

    # mutation_lst = ["AG", "AG_adar", "AG_nonadar", "GA", "UC", "CU"]
    #
    # rv_fitness_data = post_data_fitness(input_dir + "/AccuNGS/190627_RV_CV/RVB14/fits/output/fitness/syn_mutations_%s" % passages)
    # rv_fitness_data = rv_fitness_data.rename(columns={"pos": "Pos", "inferred_w": "Fitness"})
    # rv_protein_region = [629, 835, 1621, 2329, 3196, 3634, 3925, 4915, 5170, 5239, 5785, 7168]
    # rv_fitness_data = add_Protein_to_pd_df(rv_fitness_data, rv_protein_region)
    # rv_fitness_data = rv_fitness_data[rv_fitness_data["Pos"] > 3635] # just 2B
    # rv_fitness_data = rv_fitness_data.loc[rv_fitness_data.Protein != "3'UTR"]
    # # rv_fitness_data = rv_fitness_data[rv_fitness_data["Filter"] == True]
    #
    # rv_p5_p12_fitness_data = post_data_fitness(input_dir + "/AccuNGS/190627_RV_CV/RVB14/fits/output/fitness/syn_mutations_%s" % p5_p12)
    # rv_p5_p12_fitness_data = rv_p5_p12_fitness_data.rename(columns={"pos": "Pos", "inferred_w": "Fitness"})
    # rv_p5_p12_fitness_data["label"] = "RVB14_p5-p12"
    # rv_p5_p12_fitness_data = add_Protein_to_pd_df(rv_p5_p12_fitness_data, rv_protein_region)
    #
    # rv_p5_p12_fitness_data = rv_p5_p12_fitness_data[rv_p5_p12_fitness_data["Pos"] > 3635] # just 2B
    # rv_p5_p12_fitness_data = rv_p5_p12_fitness_data.loc[rv_p5_p12_fitness_data.Protein != "3'UTR"]
    #
    # cv_fitness_data = post_data_fitness(input_dir + "/AccuNGS/190627_RV_CV/CVB3/fits/output/fitness/syn_mutations_%s" % passages)
    # cv_fitness_data = cv_fitness_data.rename(columns={"pos": "Pos", "inferred_w": "Fitness"})
    # cv_protein_region = [739, 948, 1737, 2451, 3303, 3744, 4041, 5028, 5295, 5361, 5910, 7299]
    # cv_fitness_data = add_Protein_to_pd_df(cv_fitness_data, cv_protein_region)
    # cv_fitness_data = cv_fitness_data[cv_fitness_data["Pos"] > 3745] # just 2B
    # cv_fitness_data["Pos"] = cv_fitness_data["Pos"]-110 #align CVB3 and RVB14
    #
    # cv_p5_p12_fitness_data = post_data_fitness(input_dir + "/AccuNGS/190627_RV_CV/CVB3/fits/output/fitness/syn_mutations_%s" % p5_p12)
    # cv_p5_p12_fitness_data = cv_p5_p12_fitness_data.rename(columns={"pos": "Pos", "inferred_w": "Fitness"})
    # cv_p5_p12_fitness_data["label"] = "CVB3_p5-p12"
    # cv_p5_p12_fitness_data = add_Protein_to_pd_df(cv_p5_p12_fitness_data, cv_protein_region)
    # cv_p5_p12_fitness_data = cv_p5_p12_fitness_data[cv_p5_p12_fitness_data["Pos"] > 3745] # just 2B
    # cv_p5_p12_fitness_data["Pos"] = cv_p5_p12_fitness_data["Pos"]-110 #align CVB3 and RVB14
    #
    # opv_fitness_data = post_data_fitness(input_dir + "/CirSeq/PV/OPV/fits/output/fitness/syn_mutations_%s" % opv_passages)
    # opv_fitness_data = opv_fitness_data.rename(columns={"pos": "Pos", "inferred_w": "Fitness"})
    # opv_protein_region = [748, 954, 1767, 2481, 3384, 3831, 4122, 5109, 5370, 5436, 5985, 7368]
    #
    # opv_fitness_data = add_Protein_to_pd_df(opv_fitness_data, opv_protein_region)
    #
    # opv_fitness_data = opv_fitness_data.loc[opv_fitness_data.Pos >= 3832]  # just 2B
    # opv_fitness_data["Pos"] = opv_fitness_data["Pos"] - 197
    # opv_fitness_data = opv_fitness_data.loc[opv_fitness_data.Protein != "3'UTR"]
    #
    # fitness_data = pd.concat([rv_fitness_data, rv_p5_p12_fitness_data, cv_fitness_data, cv_p5_p12_fitness_data, opv_fitness_data], sort=False)
    # fitness_data = fitness_data.rename(columns={"label": "Virus"})
    # fitness_data["Mutation"] = fitness_data["Mutation"].apply(lambda x: x[0]+">"+x[1:])# if len(x)==2 else x[0]+">"+x[1:])
    # fitness_data["Mutation"] = fitness_data["Mutation"].apply(lambda x: x.split("_")[0] + "\n" + x.split("_")[-1] + "-like" if len(x)>3 else x)
    # fitness_data["Mutation"] = np.where(fitness_data["Mutation"] == "A>G", "A>G\nall", fitness_data["Mutation"])
    # # fitness_data = cv_fitness_data.rename(columns={"inferred_w": "Fitness"})
    # # print(all_data.to_string())
    # # rv_fitness_data = rv_fitness_data.rename(columns={"inferred_w": "Fitness"})
    # # rv_fitness_data["inferred_w"] = round(rv_fitness_data["inferred_w"], 2)
    # # rv_fitness_data["Fitness"] = rv_fitness_data["inferred_w"].apply(lambda x: x_round(x))
    #
    #
    # # print(rv_fitness_data.to_string())
    #
    #
    #
    # #Plots
    # #TODO: Normalize the number of pos per virus
    #
    # bins = np.arange(0, 2, 0.1)
    # g2 = sns.FacetGrid(data=fitness_data, col="Mutation", col_order=mutation_order
    #                    , col_wrap=3,  legend_out=True, hue="Virus", palette="Set2")
    # kws = dict(linewidth=.5, edgecolor="w", alpha=0.5)
    # g2 = (g2.map(plt.hist, "Fitness", **kws, bins=bins).set_titles("{col_name}").add_legend())
    #
    # # g2.fig.suptitle("Distribution of Fitness Effect", fontsize=16, weight='bold', ha='center', va='bottom')
    # g2.savefig(output_dir + "/20191006_DFE_RVB14_CVB3_OPV.png", dpi=300)
    # plt.close()
    #
    # # # pal = dict(CVB3="seagreen", RVB14="gray")
    # g3 = sns.FacetGrid(fitness_data, row="Mutation", hue="Virus", row_order=mutation_order, col="Virus", palette="Set2")#, col_wrap=3)
    # g3 = (g3.map(plt.scatter, "Pos", "Fitness", **kws).add_legend())
    # # g3.fig.suptitle(" RVB14 and CVB3 (%s)" % passages, fontsize=16, weight='bold', ha='center', va='bottom')
    # g3.savefig(output_dir + "/20191006_Mutation_Pos_RVB14_CVB3_OPV.png", dpi=300)
    # plt.close()
    #
    # g3 = sns.FacetGrid(fitness_data, col="Mutation", hue="Virus", col_order=mutation_order, palette="Set2", col_wrap=3)
    # g3 = (g3.map(plt.scatter, "Pos", "Fitness", **kws).add_legend())
    # # g3.fig.suptitle(" RVB14 and CVB3 (%s)" % passages, fontsize=16, weight='bold', ha='center', va='bottom')
    # g3.savefig(output_dir + "/20191006_Mutation_Pos_RVB14_CVB3_OPV_in_one_plot.png", dpi=300)
    # plt.close()
    #
    # # with sns.plotting_context(rc={"xVirus.fontsize": 8}):
    # g4 = sns.catplot(x="Mutation", y="Fitness", data=fitness_data, col="Protein", hue="Virus",
    #                  col_wrap=3, kind="box", order=mutation_order, palette="Set2")
    # g4.set_xticklabels(fontsize=10, rotation=45)
    #     # g4 = (g4.map(plt.boxplot, ).add_legend())
    # # g4.fig.suptitle("Fitness vs Protein", fontsize=16, weight='bold', ha='center', va='bottom')
    # g4.savefig(output_dir + "/20191006_Mutation_Protein_RVB14_CVB3_OPV_box.png", dpi=300)
    # plt.close()
    #
    # stat_df = fitness_data.loc[fitness_data.Protein == "3C"]
    # stat_df.to_csv(output_dir + "/3B_df.csv", sep=',', encoding='utf-8')
    # # stat_df = stat_df.loc[stat_df.Mutation == "AG_adar"]
    # # g4_stat = sns.boxplot(x="label", y="Fitness", data=stat_df, palette="Set2", hue="Mutation", hue_order=mutation_order)
    # #     # add_stat_annotation(g4_stat, data=stat_df, x="label", y="Fitness", hue="Mutation", hue_order=mutation_order,
    # #     #                     boxPairList=[(("RVB14", "AG_adar"), ("RVB14", "AG")), (("RVB14", "AG_adar"),("RVB14", "AG_nonadar")),
    # #     #                                  (("RVB14", "AG_adar"),("RVB14", "UC")), (("RVB14", "AG_adar"),("RVB14", "GA")),
    # #     #                                  (("RVB14", "AG_adar"),("RVB14", "CU"))],
    # #     #                     test='Mann-Whitney', textFormat='star', loc='inside', verbose=2)
    #
    # g4_stat = sns.boxenplot(x="Mutation", y="Fitness", data=stat_df, palette="Set2", hue="Virus",
    #                       order=mutation_order)
    # add_stat_annotation(g4_stat, data=stat_df, x="Mutation", y="Fitness", hue="Virus", order=mutation_order,
    #                     boxPairList=[(("A>G\nadar-like", "RVB14"), ("A>G\nadar-like", "CVB3")), (("A>G\nadar-like", "RVB14"),("A>G\nadar-like", "OPV"))],
    #                     test='Mann-Whitney', textFormat='star', loc='inside', verbose=2)
    # g4_stat.set_xticklabels(g4_stat.get_xticklabels(), rotation=45)
    # g4_stat.set(xlabel="")
    # plt.subplots_adjust(right=0.8)
    # plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))
    # plt.tight_layout()
    #
    # # g4_stat.set_xticklabels(fontsize=10)
    #     # g4 = (g4.map(plt.boxplot, ).add_legend())
    # # g4_stat.fig.suptitle("Fitness in 3B Protein", fontsize=16, weight='bold', ha='center', va='bottom')
    # plt.savefig(output_dir + "/20191006_Mutation_Protein_RVB14_CVB3_OPV_box_stat_3C.png", dpi=300)
    # # plt.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/20191006_Mutation_Protein_RVB14_CVB3_OPV_box_stat_3C.png", dpi=300)
    # plt.close()
    #
    # #
    # g5 = sns.FacetGrid(fitness_data, col="Mutation", hue="Virus", col_order=mutation_order, col_wrap=3)
    # g5 = (g5.map(qqplot, "Fitness", "Pos", **kws).add_legend())
    # g5.fig.suptitle("qqplot" , fontsize=16, weight='bold', ha='center', va='bottom')
    # g5.savefig(output_dir + "/20191006_qqplot_RVB14_CVB3_OPV.png", dpi=300)
    # plt.close()
    #
    # rvb14 = fitness_data.loc[fitness_data.Virus == "RVB14"]
    # cvb3 = fitness_data.loc[fitness_data.Virus == "CVB3"]
    # opv = fitness_data.loc[fitness_data.Virus == "OPV"]
    # ax = sns.kdeplot(rvb14.Fitness, label="RVB14")
    # ax = sns.kdeplot(cvb3.Fitness, label="CVB3")
    # ax = sns.kdeplot(opv.Fitness, label="OPV")
    # plt.legend(title='Virus', loc='upper right')
    # plt.suptitle("Kernel Density Estimation of the Fitness", fontsize=16, weight='bold', ha='center')
    # plt.xlabel("Fitness")
    # plt.ylabel("Density")
    # plt.axvline(1, color='red', ls='--')
    # # plt.tight_layout()
    # plt.savefig(output_dir + "/20191006_KDE_RVB14_CVB3_OPV.png", dpi=300)
    # plt.close()
    #
    # rvb14 = fitness_data.loc[fitness_data.Virus == "RVB14"]
    # rvb14_AG = rvb14.loc[rvb14.Mutation == "A>G\nall"]
    # rvb14_AG_adar = rvb14.loc[rvb14.Mutation == "A>G\nadar-like"]
    # rvb14_AG_nonadar = rvb14.loc[rvb14.Mutation == "A>G\nnonadar-like"]
    # rvb14_UC = rvb14.loc[rvb14.Mutation == "U>C"]
    # rvb14_GA = rvb14.loc[rvb14.Mutation == "G>A"]
    # rvb14_CU = rvb14.loc[rvb14.Mutation == "C>U"]
    #
    # cvb3 = fitness_data.loc[fitness_data.Virus == "CVB3"]
    # cvb3_AG = cvb3.loc[cvb3.Mutation == "A>G\nall"]
    # cvb3_AG_adar = cvb3.loc[cvb3.Mutation == "A>G\nadar-like"]
    # cvb3_AG_nonadar = cvb3.loc[cvb3.Mutation == "A>G\nnonadar-like"]
    # cvb3_UC = cvb3.loc[cvb3.Mutation == "U>C"]
    # cvb3_GA = cvb3.loc[cvb3.Mutation == "G>A"]
    # cvb3_CU = cvb3.loc[cvb3.Mutation == "C>U"]
    #
    # opv = fitness_data.loc[fitness_data.Virus == "OPV"]
    # opv_AG = opv.loc[opv.Mutation == "A>G\nall"]
    # opv_AG_adar = opv.loc[opv.Mutation == "A>G\nadar-like"]
    # opv_AG_nonadar = opv.loc[opv.Mutation == "A>G\nnonadar-like"]
    # opv_UC = opv.loc[opv.Mutation == "U>C"]
    # opv_GA = opv.loc[opv.Mutation == "G>A"]
    # opv_CU = opv.loc[opv.Mutation == "C>U"]
    #
    # #["AG", "AG_adar", "AG_nonadar", "UC", "GA", "CU"]
    # fig, axes = plt.subplots(2, 3, sharey='all', sharex='all')
    # sns.kdeplot(rvb14_AG.Fitness, label="RVB14", ax=axes[0, 0], legend=False)
    # sns.kdeplot(cvb3_AG.Fitness, label="CVB3", ax=axes[0, 0], legend=False)
    # sns.kdeplot(opv_AG.Fitness, label="OPV", ax=axes[0, 0], legend=False)
    #
    # sns.kdeplot(rvb14_AG_adar.Fitness, label="RVB14", ax=axes[0, 1], legend=False)
    # sns.kdeplot(cvb3_AG_adar.Fitness, label="CVB3", ax=axes[0, 1], legend=False)
    # sns.kdeplot(opv_AG_adar.Fitness, label="OPV", ax=axes[0, 1], legend=False)
    #
    # sns.kdeplot(rvb14_AG_nonadar.Fitness, label="RVB14", ax=axes[0, 2], legend=False)
    # sns.kdeplot(cvb3_AG_nonadar.Fitness, label="CVB3", ax=axes[0, 2], legend=False)
    # sns.kdeplot(opv_AG_nonadar.Fitness, label="OPV", ax=axes[0, 2], legend=False)
    #
    # sns.kdeplot(rvb14_UC.Fitness, label="RVB14", ax=axes[1, 0], legend=False)
    # sns.kdeplot(cvb3_UC.Fitness, label="CVB3", ax=axes[1, 0], legend=False)
    # sns.kdeplot(opv_UC.Fitness, label="OPV", ax=axes[1, 0], legend=False)
    #
    # sns.kdeplot(rvb14_GA.Fitness, label="RVB14", ax=axes[1, 1], legend=False)
    # sns.kdeplot(cvb3_GA.Fitness, label="CVB3", ax=axes[1, 1], legend=False)
    # sns.kdeplot(opv_GA.Fitness, label="OPV", ax=axes[1, 1], legend=False)
    #
    # sns.kdeplot(rvb14_CU.Fitness, label="RVB14", ax=axes[1, 2], legend=False)
    # sns.kdeplot(cvb3_CU.Fitness, label="CVB3", ax=axes[1, 2], legend=False)
    # sns.kdeplot(opv_CU.Fitness, label="OPV", ax=axes[1, 2], legend=False)
    #
    #
    # axes[0, 0].title.set_text('A>G\nall')
    # axes[0, 1].title.set_text('A>G\nadar-like')
    # axes[0, 2].title.set_text('A>G\nnonadar-like')
    # axes[1, 0].title.set_text('U>C')
    # axes[1, 1].title.set_text('G>A')
    # axes[1, 2].title.set_text('C>U')
    #
    #
    # # axes[1, 0].tick_params(axis="x", rotation=45)#, labelsize=9
    # # axes[1, 1].tick_params(axis="x", rotation=45)#, labelsize=9
    # # axes[1, 2].tick_params(axis="x", rotation=45)#, labelsize=9
    #
    # plt.xticks(np.arange(0, 2, step=0.5))
    #
    # axes[0, 0].axvline(1, color='red', ls='--')
    # axes[0, 1].axvline(1, color='red', ls='--')
    # axes[0, 2].axvline(1, color='red', ls='--')
    # axes[1, 0].axvline(1, color='red', ls='--')
    # axes[1, 1].axvline(1, color='red', ls='--')
    # axes[1, 2].axvline(1, color='red', ls='--')
    #
    # # plt.legend(title='Virus', loc="upper right", bbox_to_anchor=(1.05, 1.05), borderaxespad=0.)
    # # plt.suptitle("Kernel Density Estimation of the Fitness per mutation", fontsize=16, weight='bold', ha='center')
    # suplabel("x", "Fitness")
    # suplabel("y", "Density")
    #
    # # Put a legend to the right of the current axis
    # plt.legend(title='Virus', bbox_to_anchor=(1.04, 1), loc="upper left")
    # fig.subplots_adjust(right=0.8)
    # plt.tight_layout()
    # plt.savefig(output_dir + "/20191006_KDE_Mutations.png", dpi=300)
    # # plt.savefig("/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/20190912_KDE_Mutations.png", dpi=300)
    # plt.close()
    # #
    # #
    # # # rv_fitness_syn_data = syn_fitness(input_dir + "/RVB14/fits/output/fitness/%s/synonymous" % passages)
    # # # cv_fitness_syn_data = syn_fitness(input_dir + "/CVB3/fits/output/fitness/%s/synonymous" % passages)
    # # #
    # # # fitness_data = pd.concat([rv_fitness_syn_data, cv_fitness_syn_data], sort=False)
    # # # print(all_data.to_string())
    # # # fitness_data = fitness_data.rename(columns={"allele1": "Fitness"})
    # # #
    # # # g3 = sns.violinplot(x="Mutation", hue="Virus", y="Fitness", data=fitness_data, split=True, cut=0,
    # # #                     order=["AG", "adar", "nonadar", "UC", "GA", "CU"])
    # # # g3.set_title("Fitness distribution of synonymous mutations in RVB14 and CVB3\n(%s)" % passages)
    # # # plt.tight_layout()
    # # # # plt.show()
    # # # plt.savefig(output_dir + "/20190901_posterior_dist_fitness_synon.png", dpi=300)
    # # # plt.close()


if __name__ == "__main__":
    main()
