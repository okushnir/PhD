#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import sys, argparse
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re
import glob
import Bio.Seq as Seq
import matplotlib.gridspec as gridspec
import time
from Bio import Entrez
from Bio import SeqIO

def find_a_gene(ncbi_id):
    """
    :param ncbi_id: NCBI_ID number
    :return:
    """
    try:
        Entrez.email = "A.N.Other@example.com"
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=ncbi_id)
        ncbi_gb = SeqIO.read(handle, "gb")
        handle.close()
        start_pos = 0
        end_pos = 0
        for feature in ncbi_gb.features:
            if feature.type == 'CDS':
                print(feature.qualifiers['gene'])
                if (feature.qualifiers['gene'] == ['N']) | (feature.qualifiers['gene'] == "N"):
                    start_pos = feature.location.start + 1
                    end_pos = feature.location.end + 0
                    protein_seq = str(feature.qualifiers['translation']).strip('[]')
                    return start_pos, end_pos, protein_seq
    except Exception as e:
        print(e)

def main():
    # sars = "NC_004718"
    # find_a_gene(sars)
    # "con229e": "NC_002645", "mers": "NC_019843",
    cov_dic = {"sars": "NC_004718"} #"oc43": "NC_006213","sars_cov2": "LC534419"

    for key in cov_dic:
        start_pos, end_pos, protein_seq = find_a_gene(cov_dic[key])
        print("%s: start_pos=%s; end_pos=%s; protein_seq=%s" % (cov_dic[key], start_pos, end_pos, protein_seq))







if __name__ == "__main__":
    main()

# def main(args):
    # parser = argparse.ArgumentParser()
    # parser.add_argument("virus", type=str, help="name of the virus, RVB14")
    # parser.add_argument("passages", type=str, help="from which passages, p0-p12")
    # parser.add_argument("without", type=int, help="Exclude passage no.")
    # parser.add_argument("input_dir", type=str, help="the path to the directory that contains data_mutation.csv")
    # parser.add_argument("quality", type=str, help="what is the prefix for the data_mutation.csv file; quality of the pipline ; for example: q38")
    # parser.add_argument("mutation_type", type=str, help="all/syn")
    # args = parser.parse_args(sys.argv[1:])
    # main(args)