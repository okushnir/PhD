#! /usr/local/python/anaconda_python-3.6.1

"""
@Author: odedkushnir

"""


import glob
import pandas as pd
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
    list = []
    columns = ["median", "MAD", "min", "max", "pval", "mutation"]
    all_df = pd.DataFrame(columns=columns)
    for file in files_lst:
        with open(file, "r") as file_handler:
            lines = file_handler.readlines()
            # print(lines)
            mutation = file.split("/")[-1].split(".txt")[0].split("_")[-1]
            median = lines[-2].split("     ")[2].split("    ")[0]
            MAD = lines[-2].split("     ")[2].split("    ")[1]
            min = lines[-2].split("     ")[2].split("    ")[2]
            max = lines[-2].split("     ")[2].split("    ")[3]
            pval = lines[-2].split("     ")[2].split("    ")[4]
            list.append([median, MAD, min, max, pval, mutation])
            df = pd.DataFrame(list, columns=columns, index=range(1))
            all_df = pd.concat([all_df, df])
            list = []
    all_df.reset_index(inplace=True)
    all_df.drop(["index"], axis=1, inplace=True)
    return all_df


def main():

    input_dir = "/Users/odedkushnir/Projects/fitness"
    rvb14 = input_dir + "/AccuNGS/190627_RV_CV/RVB14/fits/output/mutation/p0-p12"
    cvb3 = input_dir + "/AccuNGS/190627_RV_CV/CVB3/fits/output/mutation/p0-p12"
    opv = input_dir + "/CirSeq/PV/OPV/fits/output/mutation/p1-p7"
    output_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/fits/all_viruses_mutation rates"
    rvb14_rates = fitness_parameters(rvb14, output_dir)
    rvb14_rates["label"] = rvb14.split("/")[-5]
    rvb14_rates.sort_values(by="mutation", ascending=True, inplace=True)
    cvb3_rates = fitness_parameters(cvb3, output_dir)
    cvb3_rates["label"] = cvb3.split("/")[-5]
    cvb3_rates.sort_values(by="mutation", ascending=True, inplace=True)
    opv_rates = fitness_parameters(opv, output_dir)
    opv_rates["label"] = opv.split("/")[-5]
    opv_rates.sort_values(by="mutation", ascending=True, inplace=True)
    columns = ["label", "median", "MAD", "min", "max", "pval", "mutation"]
    df = pd.concat([rvb14_rates, cvb3_rates, opv_rates])
    df = df[columns]
    print(df.to_string())
    df.to_csv(output_dir + "/viruses_df.csv", sep=',', encoding='utf-8')


if __name__ == "__main__":
    main()
