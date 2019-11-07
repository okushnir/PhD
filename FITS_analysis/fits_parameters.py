#! /usr/local/python/anaconda_python-3.6.1

"""
@Author: odedkushnir

"""


import glob
import subprocess
import os
from Context_analysis_RV import checkKey


def fitness_parameters(input_dir, output_dir):
    """
    Creates the fitness parameters file for each mutation
    :param input_dir: the output directory for fitness -mutation
    :param output_dir: the input directory for FITS
    :return: Dictionary for all mutation and their mutation rates
    """
    files_lst = glob.glob(input_dir + "/summary_mutation_syn*")
    rate_dict = {}
    for file in files_lst:
        with open(file, "r") as file_handler:
            lines = file_handler.readlines()
            # print(lines)
            mutation = file.split("/")[-1].split(".txt")[0].split("_")[-1]
            rate = lines[15].split("     ")[2].split("    ")[0]
            adar_rev = lines[16].split("     ")[2].split("   ")[0]
            if adar_rev.find("*") == 0:
                adar_rev = adar_rev.split("*")[-1]
            rate_dict[mutation] = rate
            if (mutation == "adar") | (mutation == "nonadar"):
                mutation = "AG_" + file.split("/")[-1].split(".txt")[0].split("_")[-1]
                rate_dict[mutation] = rate_dict[mutation.split("AG_")[-1]]
                del rate_dict[mutation.split("AG_")[-1]]
                # print(mutation)
                rev_mutation = mutation[::-1]
                rate_dict[rev_mutation] = adar_rev
    print(rate_dict)

    for mutation in rate_dict:
        if (mutation == "AG_adar") | (mutation == "AG_nonadar"):
            with open(output_dir + "/parameters_fitness_%s.txt" % mutation, "w+") as out_handler:
                out_handler.writelines("verbose 1\n"
                                       "N 756000\n"
                                       "fitness_prior smoothed_composite\n"
                                       "mutation_rate0_1 %s\n"
                                       "mutation_rate1_0 %s\n"
                                       "min_fitness_allele0 1.0\n"
                                       "max_fitness_allele0 1.0\n"
                                       "min_fitness_allele1 0.0\n"
                                       "max_fitness_allele1 2.0\n"
                                       "num_samples_from_prior 100000\n"
                                       "acceptance_rate 0.01" % (rate_dict["AG_nonadar"], rate_dict["radanon_GA"]))
        elif (mutation != "rada_GA") & (mutation != "radanon_GA"):
            rev_mutation = mutation[::-1]
            with open(output_dir + "/parameters_fitness_%s.txt" % mutation, "w+") as out_handler:
                out_handler.writelines("verbose 1\n"
                                        "N 756000\n"
                                        "fitness_prior smoothed_composite\n"
                                        "mutation_rate0_1 %s\n"
                                        "mutation_rate1_0 %s\n"
                                        "min_fitness_allele0 1.0\n"
                                        "max_fitness_allele0 1.0\n"
                                        "min_fitness_allele1 0.0\n"
                                        "max_fitness_allele1 2.0\n"
                                        "num_samples_from_prior 100000\n"
                                        "acceptance_rate 0.01" % (rate_dict[mutation], rate_dict[rev_mutation]))


    return rate_dict

def main():
    virus = "RVB14"
    passages = "p0-p12"
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/190627_RV_CV/merged/%s/fits/output/mutation/%s_all_positions" % (
    virus, passages)
    output_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/190627_RV_CV/merged/%s/fits/input/%s" % (
    virus, passages)

    fitness_parameters(input_dir, output_dir)

if __name__ == "__main__":
    main()
