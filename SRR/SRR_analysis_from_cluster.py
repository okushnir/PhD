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
import subprocess
import pbs_jobs
from time import sleep
#import make_reference_from_consensus
from os import listdir
from os.path import isfile, join


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

def script_runner(cmds, alias = "script", load_python=False):
    """
    run script on cluster
    :param cmds: script running line
    :param alias: job name (default: script)
    :return: job id
    """
    cmdfile = pbs_jobs.get_cmdfile_dir("script", alias); tnum=1; gmem=2
    print(cmdfile, alias, tnum, gmem, cmds)
    pbs_jobs.create_pbs_cmd(cmdfile, alias=alias, gmem=gmem, cmds=cmds, load_python=load_python)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def check_pbs(job_id):
    """
    :param job_id: The PBS job id
    :return: "Done!", when the job is done
    """
    status = "Running..."
    try:
        process = subprocess.check_output("qstat | grep " + str(job_id), shell=True)
        while process != "":
            process = subprocess.check_output("qstat | grep " + str(job_id), shell=True)
            sleep(0.30)
        print("")
    except (subprocess.CalledProcessError):
        process = ""
    if process == "":
        status = "Done"
    return status


def subtract_file(path):
    with open(path, "r") as handle:
        file = handle.readlines()
        print("NO. of lines in the fastq file are: %i" % (len(file)))
        no_line = len(file) % 4
        no_line = -(no_line)
        print("NO. of lines to cut are: %i" % (no_line))
        if no_line != 0:
            file = file[0:no_line]
            # print(file)
            with open(path, "w") as output_handle:
                output_handle.writelines(file)
        else:
            return path

def does_file_exist_in_dir(path):
    return any(isfile(join(path, i)) for i in listdir(path))


def download_SRR(SRR_dir):
    """
    :param SRR_dir: path directory of the SRA's to analyze. should contain at least 1 dir of SRA whcih contains SraRunTable.txt
    :return: pd.Series of job id in the qsub
    """
    print("check1")
    dirs = glob.glob(SRR_dir + "/SRP*")
    job_status = pd.Series()
    for SRR in dirs:
        print(SRR)
        file_path = glob.glob(SRR)
        """ Read the sra table"""
        sra_table = pd.read_table(file_path[0] + "/SraRunTable.txt")

        """ Filter the table by organism and MBases"""
        pattern = re.compile("[E][c](.){1,}|[E][n](.){1,}|[C](.){1,}|[R][h](.){1,}|[H][u](.){1,}[r][h](.){1,}|[H][u](.){1,}[p][o](.){1,}", re.DOTALL)
        sra_table['Interest'] = sra_table['Organism'].str.contains(pattern)
        sra_table = sra_table[(sra_table['Interest'] == True)]
        sra_table = sra_table.sort_values(by=['MBases'], ascending=False)
        sub_sra_table = sra_table.head(10)
        sub_sra_table.to_csv(file_path[0] + "/Top10SRR1.csv", sep=',', encoding='utf-8')

        """Download fastqs"""
        try:
            sub_sra_table["Run"].apply(lambda x: os.mkdir(SRR + "/" + x))
            sub_sra_table["Run"].apply(lambda x: os.mkdir(SRR + "/" + x + "/fastq"))
        except Exception as e:
            print("type error: " + str(e))
        sub_sra_table["fastq_dir"] = SRR + "/" + sub_sra_table["Run"] + "/fastq/"
    #
    #
    #     fastq_dump = "HOME=/sternadi/home/volume3/okushnir/\n" \
    #                  "cd /sternadi/home/volume3/okushnir\n" \
    #                  "/sternadi/home/volume3/okushnir/Tools/sratoolkit.2.9.2-centos_linux64/bin/fastq-dump --split-files "
    #     print("Downloading fastq files...")
    #     sub_sra_table["fastq_exist"] = sub_sra_table.apply(lambda x: does_file_exist_in_dir(x["fastq_dir"]), axis=1)
    #
    #     if sub_sra_table.iloc[0]["fastq_exist"] == False:
    #         sub_sra_table["job_id"] = sub_sra_table["Run"].apply(
    #             lambda x: script_runner(fastq_dump + x + " -O " + file_path[0] + "/" + x + "/fastq"))
    #         # sub_sra_table["job_id"].apply(lambda x: (check_pbs(x) != "Done"))
    #         # print(sub_sra_table)
    #         print("Done!")
    #         # return sub_sra_table, sub_sra_table["job_id"]
    #     else:
    #         # sub_sra_table["fastq_file"] = sub_sra_table.apply(lambda x: glob.glob(x["fastq_dir"] + "*.fastq")[0], axis=1)
    #         # sub_sra_table["fastq_file"].apply(lambda x: subtract_file(x))
    #         # sub_sra_table["job_id"] = 0
    #         print("Done!")
    #     sub_sra_table.to_csv(file_path[0] + "/Top10SRR2.csv", sep=',', encoding='utf-8')
    #     job_status.append(sub_sra_table["job_id"])
    #     # download_refs_and_run_pipeline(sub_sra_table, SRR_dir)
    # return job_status


def download_refs_and_run_pipeline(sub_sra_table, SRR_dir):
    """
    :param sub_sra_table:
    :param SRR_dir:
    :return:
    """
    srr_dic = {"SRR2558452": "KJ170553", "SRR2558453": "KJ170542", "SRR2558454": "KJ170565",
               "SRR2558455": "KJ170571", "SRR2558456": "KJ170550", "SRR2558457": "KJ170560",
               "SRR2558518": "KJ170670", "SRR2558529": "KJ170579", "SRR2558449": "KJ170541",
               "SRR2558492": "KJ170587"}
    org_dic = {"CVB3": "M16572", "RVB14": "NC_001490", "Echovirus E7": "MH732737", "Echovirus E6": "JX976771",
               "Coxsackievirus A16": "NC_001612", "Enterovirus A": "NC_001612", "Echovirus E3": "KX808644",
               "Coxsackievirus B2": "AF081485", "Echovirus E25": "KJ957190", "Human poliovirus 1": "V01149",
               "Human poliovirus 2": "V01149", "Human poliovirus 3": "V01149", "Enterovirus C": "V01149",
               "Enterovirus D68": "NC_001430", "Enterovirus D": "NC_001430", "Rhinovirus B": "NC_001490",
               "Coxsackievirus B3 (strain Nancy)": "JN048468", "Rhinovirus C": "LC428177"}

    print("check2")
    print(SRR_dir)
    dirs = glob.glob(SRR_dir + "*RP*")
    for SRR in dirs:
        print(SRR)
        file_path = glob.glob(SRR)
        """Reference download"""
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

        print("Fetching References...")
        # TO ADD - check if file exists....
        sub_sra_table.apply(lambda x: find_srr_ref(x["ncbi_id"], x["ref_dir"]), axis=1)
        print("DONE Fetching.")

        sub_sra_table.to_csv(file_path[0] + "/Top10SRR3.csv", sep=',', encoding='utf-8')
        print(sub_sra_table)

        """Run pipline using qsub """
        try:
            # Makes the dir for the pipeline output
            sub_sra_table["Run"].apply(lambda x: os.mkdir(SRR + "/" + x + "/q30_consensus_1e-03/"))

        except Exception as e:
            print(e)

        sub_sra_table["output_dir"] = SRR + "/" + sub_sra_table["Run"] + "/q30_consensus_1e-03/"
        sub_sra_table["ref"] = sub_sra_table.apply(lambda x: glob.glob(x["ref_dir"] + "*")[0], axis=1)
        sub_sra_table.to_csv(file_path[0] + "/Top10SRR4.csv", sep=',', encoding='utf-8')
        cmd = "python /sternadi/home/volume1/shared/SternLab/pipeline_runner.py"

        sub_sra_table["ref_exist"] = sub_sra_table.apply(lambda x: does_file_exist_in_dir(x["ref_dir"]), axis=1)
        if sub_sra_table.iloc[0]["ref_exist"] == False:
            sub_sra_table["job_id"] = sub_sra_table.apply(lambda x: script_runner(cmd + " -i " + x["fastq_dir"] + " -o " + x["output_dir"] + " -r " +
                                                    x["ref"] + " -s 1 -e 4 -t f -g Y -NGS_or_Cirseq 1 -q 30 -rep 1 -ev "
                                                               "1e-03"), axis=1)
            print("Done!")
        else:
            print("Done!")
    return sub_sra_table


def ref_from_con(sub_sra_table, SRR_dir):
    """

    :param sub_sra_table: path of the TOP10SRR3.csv
    :param SRR_dir:
    :return:
    """
    dirs = glob.glob(SRR_dir + "*RP*")
    for SRR in dirs:
        print(SRR)
        file_path = glob.glob(SRR)
        sub_sra_table = pd.read_csv(sub_sra_table, sep=',', encoding='utf-8')
        sub_sra_table["freqs_path"] = sub_sra_table.apply(lambda x: glob.glob(x["output_dir"] + "*.freqs")[0], axis=1)

        try:
            sub_sra_table["Run"].apply(lambda x: os.mkdir(SRR + "/" + x + "/ref_from_conc_dir/"))
        except Exception as e:
            print("type error: " + str(e))

        sub_sra_table["ref_from_conc_dir"] = SRR + "/" + sub_sra_table["Run"] + "/ref_from_conc/"
        cmd = "python /sternadi/home/volume3/okushnir/Python_script/make_reference_from_consensus.py"
        sub_sra_table.apply(lambda x: script_runner(cmd + " -f " + x["ref"] + "-p" + x["freqs_path"] + "-o" +
                                                        x["ref_from_conc_dir"] + "-c 1000", alias="Reference_Con"), axis=1)

    return sub_sra_table


def main():
    # SRR_dir = "/sternadi/home/volume3/okushnir/Pra_SRA2/"

    #Local
    SRR_dir = "/Users/odedkushnir/Projects/fitness/SRA"


    print(SRR_dir)

    job_status = download_SRR(SRR_dir)
    print(job_status)
    # job_status.apply(lambda x: (check_pbs(x) != "Done"))
    # print(dataframe)
    # print(job_lst)

    # job_status = pd.Series()

    # print(job_status)

    # dataframe["check_id"] = dataframe["job_id"].apply(lambda x: check_pbs(x), axis=1)


    # dataframe = pd.read_csv("/sternadi/home/volume3/okushnir/Pra_SRA2/SRP155896_Entero_D/Top10SRR2.csv", sep=',', encoding='utf-8')
    # download_refs_and_run_pipeline(dataframe, SRR_dir)


if __name__ == "__main__":
    main()