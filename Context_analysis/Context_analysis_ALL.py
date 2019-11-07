#! /usr/local/python/anaconda_python-3.6.1

"""
@Author: odedkushnir

"""


import os.path
from sequnce_utilities import *
import glob
import urllib



def checkKey(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        raise Exception()

def main():
    # for Cluster
    # parser = OptionParser("usage: %prog [options]")
    # parser.add_option("-s", "--SRR_dir", dest="SRR_dir", help="the directory of all the SRR's")
    # parser.add_option("-n", "--SRR_no", dest="SRR_no", help="SRR Number")
    # (options, args) = parser.parse_args()
    #
    # SRR_dir = options.SRR_dir
    # SRR_No = options.SRR_no
    # freqs_file = check_dirname(freqs_file)


    #for Local
    suffix = "data_mutation.csv"
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV"

    data_RV = pd.read_csv(input_dir +"/RVB14/" + suffix)
    data_CV = pd.read_csv(input_dir +"/CVB3/" + suffix)
    # data_atcg1 = pd.read_csv(input_dir +"/ATCG1/" + suffix)
    # data_flna = pd.read_csv(input_dir +"/FLNA/" + suffix)
    data = pd.concat([data_RV, data_CV], sort=False) #, data_atcg1, data_flna

    data.to_csv(input_dir + "/data_mutation.csv", sep=',', encoding='utf-8')


# Context data


    mutation_for_rna = ["AG"]
    dataForPlotting_AG = data[(data["mutation_type"].isin(mutation_for_rna))]#& (data["Type"] == "Synonymous")]
    dataForPlotting_AG.to_csv(input_dir + "/data_XpA_by_mutation.csv", sep=',', encoding='utf-8')

    mutation_for_rna = ["UC"]
    dataForPlotting_UC = data[(data["mutation_type"].isin(mutation_for_rna))]  # & (data["Type"] == "Synonymous")]
    dataForPlotting_UC.to_csv(input_dir + "/data_UpX_by_mutation.csv", sep=',', encoding='utf-8')


if __name__ == "__main__":
    main()

