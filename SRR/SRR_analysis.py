#! /usr/local/python/anaconda_python-3.6.1

"""
@Author: odedkushnir

"""

import sys
import os
import numpy as np
import re
import glob
import pandas as pd
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO
import pbs_runners
import subprocess


def find_srr_ref(srr_id, out_dir):
    """
    :param ncbi_id: NCBI_ID number
    :return:the start and the end positions of the Coding region
    """
    try:
        Entrez.email = "A.N.Other@example.com"
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=srr_id)
        ncbi_gb = SeqIO.read(handle, "gb")
        handle.close()
        file_out = out_dir + "%s.fasta" % (srr_id)
        with open(file_out, "w") as handle:
            SeqIO.write(ncbi_gb, handle, "fasta")
            return ncbi_gb
    except Exception as e:
        print("type error: " + str(e))
        print('Failed to fetch record.')
        return Exception


def checkKey(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        raise Exception()

def main():
    SRR_dir = "/Users/odedkushnir/Projects/fitness/Pra_SRA/"
    dirs = glob.glob(SRR_dir + "*RP*")
    lst_srr = []
    for SRR in dirs:
        file_path = glob.glob(SRR)
        """ Read the sra table"""
        sra_table = pd.read_table(file_path[0] + "/SraRunTable.txt")

        """ Filter the tanle by organism and MBases"""
        pattern = re.compile("[E][c](.){1,}|[E][n](.){1,}|[C](.){1,}|[R][h](.){1,}|[H][u](.){1,}[r][h](.){1,}|[H][u](.){1,}[p][o](.){1,}", re.DOTALL)
        sra_table['Interest'] = sra_table['Organism'].str.contains(pattern)
        sra_table = sra_table[(sra_table['Interest'] == True)]
        sra_table = sra_table.sort_values(by=['MBases'], ascending=False)
        sub_sra_table = sra_table.head(10)
        sub_sra_table.to_csv(file_path[0] + "/Top10SRR.csv", sep=',', encoding='utf-8')
        """ Manualy downlaod refs"""

        srr_dic = {"SRR2558452": "KJ170553", "SRR2558453": "KJ170542", "SRR2558454": "KJ170565",
                   "SRR2558455": "KJ170571", "SRR2558456": "KJ170550", "SRR2558457": "KJ170560",
                   "SRR2558518": "KJ170670", "SRR2558529": "KJ170579", "SRR2558449": "KJ170541",
                   "SRR2558492": "KJ170587"}
        org_dic = {"CVB3": "M16572", "RVB14": "NC_001490", "Echovirus E7": "MH732737", "Echovirus E6": "JX976771",
                    "Coxsackievirus A16": "NC_001612", "Enterovirus A": "NC_001612", "Echovirus E3": "KX808644",
                   "Coxsackievirus B2": "AF081485", "Echovirus E25": "KJ957190 ", "Human poliovirus 1": "V01149",
                   "Human poliovirus 2": "V01149", "Human poliovirus 3": "V01149", "Enterovirus C": "V01149",
                   "Enterovirus D68": "NC_001430", "Enterovirus D": "NC_001430", "Rhinovirus B": "NC_001490",
                   "Coxsackievirus B3 (strain Nancy)": "JN048468"}
        try:
            sub_sra_table["Run"].apply(lambda x: os.mkdir(SRR + "/" + x))
            sub_sra_table["Run"].apply(lambda x: os.mkdir(SRR + "/" + x + "/fastq"))
        except Exception as e:
            print("type error: " + str(e))

        fastq_dump = "/Users/odedkushnir/Tools/sratoolkit.2.9.4-mac64/bin/fastq-dump --split-files "
        print("Downloading fastq files...")
        sub_sra_table["Run"].apply(lambda x: os.popen(fastq_dump + x + " -O " + file_path[0] + "/" + x + "/fastq").close())
        print("Done")
        try:
            sub_sra_table["ncbi_id"] = sub_sra_table.apply(lambda x: checkKey(srr_dic, x["Run"]), axis=1)
        except Exception as e:
            print("type error: " + str(e))
            print("The SRR is not the dictionary")
            sub_sra_table["ncbi_id"] = sub_sra_table.apply(lambda x: checkKey(org_dic, x["Organism"]), axis=1)


        try:
            sub_sra_table["Run"].apply(lambda x: os.mkdir(SRR + "/" + x + "/ref/"))
        except Exception as e:
            print("type error: " + str(e))

        sub_sra_table["ref_dir"] = SRR + "/" + sub_sra_table["Run"] + "/ref/"
        # try:
        #     sub_sra_table.apply(lambda x: find_srr_ref(x["Sample_Name"], x["ref_dir"]), axis=1)
        # except Exception as e:
        #     print("type error: " + str(e) + "HIIII")

        print("Fetching References...")
        sub_sra_table.apply(lambda x: find_srr_ref(x["ncbi_id"], x["ref_dir"]), axis=1)
        print("DONE Fetching.")

        sub_sra_table.to_csv(file_path[0] + "/Top10SRR.csv", sep=',', encoding='utf-8')
        print(sub_sra_table)

        """Run pipline using qsub """
        # 4th run pipeline:
        # folders = glob.glob("/sternadi/home/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/*")
        #
        # for d in folders:
        # output_dir = SRR + "/q30_consensus_1e-03"
        # # cmd = "python /sternadi/home/volume1/shared/SternLab/pipeline_runner.py -i %s -o %s -r %s " \
        # #       "-s 1 -e 4 -g Y -NGS_or_Cirseq 1 -q 30 -rep 1 -ev 1e-03" % (SRR, output_dir, ref_file)
        # cmd = "python /sternadi/home/volume1/shared/SternLab/pipeline_runner.py"
        # sub_sra_table.apply(lambda x: pbs_runners.script_runner(cmd + "-i %s -o %s -r %s " \
        #       "-s 1 -e 4 -g Y -NGS_or_Cirseq 1 -q 30 -rep 1 -ev 1e-03" % (SRR + "/" + x["Run"] + "/fastq/", output_dir,
        #                                                                   x["ref_dir"]), alias="pipeline_d"))
        #
        # pbs_runners.script_runner(cmd, alias="pipeline_d")

if __name__ == "__main__":
    main()