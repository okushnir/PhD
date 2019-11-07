#! /usr/local/python/anaconda_python-3.6.1

"""
@Author: odedkushnir

"""

import sys
import pandas as pd
import os
import numpy as np
import subprocess
from Context_analysis_RV import checkKey
import rpy2.robjects as robjects
import glob
from Context_analysis_RV import checkKey
from fits_plotter import *
from fits_parameters import *


def fits_data_construction(input_dir, output_dir, from_passage, to_passage):
    all = pd.read_csv(input_dir + "q38_data_mutation.csv")
    all = all[all["label"] != "RVB14-Next-RNA Control"]
    all = all[all["label"] != "RVB14-p1"]

    # all = all[all["label"] != "%s-RNA Control" % virus]
    # all = all[all["label"] != "%s-%s" % (virus, from_passage)]
    # all = all[all["label"] != "%s-%s" % (virus, from_passage)]
    all["passage"] = all["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    all["passage"] = np.where(all["passage"] == "RNA Control", 0, all["passage"])
    all["passage"] = all["passage"].astype(int)
    all = all[all["passage"] >= from_passage]
    all = all[all["passage"] <= to_passage]
    all["Type"] = all["Type"].fillna("NonCodingRegion")
    all["pval"] = all["pval"].fillna(1)
    all["Frequency"] = np.where(all["pval"] > 0.01, 0, all["Frequency"])
    all["Frequency"] = np.where(all["Prob"] < 0.95, 0, all["Frequency"])

    syn = all[all['Type'] == "Synonymous"]

    df_dict = {"syn": syn, "all": all}

    transitions = ["A>G", "G>A", "U>C", "C>U"]
    for key in df_dict:
        for mutation in transitions:
            df_mutation = checkKey(df_dict, key)
            df_mutation = df_mutation[df_mutation["Mutation"] == mutation]
            quarte_allelic_mapping = {'A': 0, 'C': 1, 'G': 2, 'U': 3}
            if mutation == "A>G":
                df_mutation["ADAR_like"] = df_mutation.Prev.str.contains('UA') | df_mutation.Prev.str.contains('AA')
                adar = df_mutation[df_mutation["ADAR_like"] == True]
                nonadar = df_mutation[df_mutation["ADAR_like"] == False]
                adar_dict = {"adar": adar, "nonadar": nonadar}
                for adar_key in adar_dict:
                    like = checkKey(adar_dict, adar_key)
                    like['Base'] = like['Base'].apply(lambda x: quarte_allelic_mapping[x])
                    like['Ref'] = like['Ref'].apply(lambda x: quarte_allelic_mapping[x])
                    like = like.rename(columns={'passage': 'gen', 'Base': 'base', 'Frequency': 'freq', 'Ref': 'ref',
                                                'Read_count': 'read_count', 'Rank': 'rank', 'Pos': 'pos'})
                    like = like[['gen', 'base', 'freq', 'pos']]
                    like = like.sort_values(['pos', 'gen', 'freq'])

                    like.to_csv(output_dir + "%s_mutations_%s_%s.csv" % (key, mutation, adar_key), index=False)

            df_mutation['Base'] = df_mutation['Base'].apply(lambda x: quarte_allelic_mapping[x])
            df_mutation['Ref'] = df_mutation['Ref'].apply(lambda x: quarte_allelic_mapping[x])
            df_mutation = df_mutation.rename(columns={'passage': 'gen', 'Base': 'base', 'Frequency': 'freq', 'Ref': 'ref',
                                    'Read_count': 'read_count', 'Rank': 'rank', 'Pos': 'pos'})
            df_mutation = df_mutation[['gen', 'base', 'freq', 'pos']]
            df_mutation = df_mutation.sort_values(['pos', 'gen', 'freq'])
            df_mutation.to_csv(output_dir + "%s_mutations_%s.csv" % (key, mutation), index=False)



def main():
    virus = "RVB14"
    passages = "p0-p12"
    from_passage = int(passages.split("p")[1].split("-")[0])
    to_passage = int(passages.split("-p")[1])
    input_dir = "/sternadi/home/volume3/okushnir/AccuNGS/190627_RV_CV/merged/%s/" % virus #"/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/%s/"
    output_dir = input_dir + "fits/input/%s/" % passages
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    fits_data_construction(input_dir,output_dir, from_passage, to_passage)

    # R-script modify 'virus' and 'passage' as previous
    test = robjects.r(''' library (plyr)
                        library(stringr)
                        virus = "RVB14" ### MODIFY
                        passage = "p0-p12" ### MODIFY
                        setwd <- setwd(sprintf("/sternadi/home/volume3/okushnir/AccuNGS/190627_RV_CV/merged/%s/fits/input/%s", virus, passage)) ### MODIFY
                        setwd <- setwd(sprintf("/sternadi/home/volume3/okushnir/AccuNGS/190627_RV_CV/merged/%s/fits/input/%s", virus, passage)) ### MODIFY

                        mutation_list <- list("A>G", "C>U", "G>A", "U>C", "A>G_nonadar", "A>G_adar")
                        type_list  <- list("all", "syn")
                        for (type in type_list){
                          for (mutation in mutation_list){
                            all.pos.df <- read.csv(sprintf("%s_mutations_%s.csv", type, mutation))
                            mutation <- str_replace_all(mutation, ">", "")
                        final.pos.df <- ddply( all.pos.df, c("gen", "pos"), function(my.subset) {
                          my.newline <- my.subset
                          my.newline$base <- 0
                          my.newline$freq <- 1 - my.subset$freq[1]
                          my.subset$base <- 1

                          rbind(my.newline, my.subset)	
                        })

                        final.pos.sorted.df <-  final.pos.df[ order( final.pos.df$pos, final.pos.df$gen, final.pos.df$base), ]
                        write.table(final.pos.sorted.df, file = sprintf("final_%sMutations_sorted_%s.txt", type, mutation), quote = F, row.names = F, sep = "\t")
                        }}''')
    mutation_lst = ["AG", "AG_adar", "AG_nonadar", "GA", "UC", "CU"]
    for mutation in mutation_lst:
        df_mutation = pd.read_table(output_dir + "final_allMutations_sorted_%s.txt" % mutation)
        for position in df_mutation["pos"].iloc[0:-1]:
            df_pos = df_mutation.groupby(["pos"]).get_group(position)
            # print(df_pos.to_string())
            os.mkdir(output_dir + "%s" % mutation)
            df_pos.to_csv(output_dir + "%s/final_allMutations_sorted_%s.txt" % (mutation, position), index=False, sep="\t")
    # Run FITS_jobarray_mutation.cmd

    # Run fits_parameters.py
    fitness_parameters(input_dir=output_dir, output_dir=input_dir)

    # Run FITS_jobarray_fitness.cmd

    # Run fits_plotter.py


    # passages = "p0-p12_with_err"
    # input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV"
    # rv_mutation_data = post_data_mutation(input_dir + "/RVB14/fits/output/mutation/%s" % passages)
    # cv_mutation_data = post_data_mutation(input_dir + "/CVB3/fits/output/mutation/%s" % passages)
    # output_dir = input_dir + "/RVB14/fits/output/fitness/%s/plots" % passages
    # try:
    #     os.mkdir(output_dir)
    # except OSError:
    #     print("Creation of the directory %s failed" % output_dir)
    # else:
    #     print("Successfully created the directory %s " % output_dir)
    # all_data = pd.concat([rv_mutation_data, cv_mutation_data], sort=False)
    # print(all_data.to_string())
    # all_data = all_data.rename(columns={"allele0_1": "Transition rate"})
    # # all_data["Transition rate"] = np.log10(all_data["Transition rate"])
    # # all_data["allele0_1"] = all_data["allele0_1"]*10**5
    #
    # g1 = sns.violinplot(x="Mutation", hue="label", y="Transition rate", data=all_data, split=True, cut=0,
    #                     order=["AG", "adar", "nonadar", "UC", "GA", "CU"])
    #
    # # g1.set_ylabel("Transition rate")
    # g1.set_yscale("log")
    # # g1.set_yscale('symlog', linthreshy=5*10**-9)
    # # yaxis = plt.gca().yaxis
    # # yaxis.set_minor_locator(MinorSymLogLocator(1e-1))
    # g1.set_title("Transition rate distribution in RVB14 and CVB3")
    # # g1.set_yticks(ticks=[10**-5, 10**-6, 0], minor=True)
    # # g1.set_ylim(10 ** -8, 10 ** -3)
    # # g1.legend(bbox_to_anchor=(1.05, 0.5), loc="upper right", borderaxespad=0., fontsize='small')
    # plt.tight_layout()
    # # plt.show()
    # plt.savefig(output_dir + "/20190901_posterior_dist.png", dpi=300)
    # plt.close()
    #
    # rv_fitness_data = post_data_fitness(input_dir + "/RVB14/fits/output/fitness/%s" % passages)
    # cv_fitness_data = post_data_fitness(input_dir + "/CVB3/fits/output/fitness/%s" % passages)
    #
    # fitness_data = pd.concat([rv_fitness_data, cv_fitness_data], sort=False)
    # print(all_data.to_string())
    # fitness_data = fitness_data.rename(columns={"allele1": "Fitness"})
    #
    # g2 = sns.violinplot(x="Mutation", hue="label", y="Fitness", data=fitness_data, split=True, cut=0,
    #                     order=["AG", "adar", "nonadar", "UC", "GA", "CU"])
    # g2.set_title("Fitness distribution in RVB14 and CVB3 (%s)" % passages)
    # plt.tight_layout()
    # # plt.show()
    # plt.savefig(output_dir + "/20190901_posterior_dist_fitness.png", dpi=300)
    # plt.close()
    #
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