
#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import os
import pbs_runners
import glob
import pandas as pd


def checkKey(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        raise Exception()


def index(csv_file, fastq_path, output_dir):
    indexing = pd.read_csv(csv_file)  # your_csv_dir
    sample_id = list(indexing.SampleID)
    index_1 = list(indexing.Index1Sequence)
    index_2 = list(indexing.Index2Sequence)

    index_dic = {}
    for i in range(len(sample_id)):
        index_dic[index_1[i] + "+" + index_2[i]] = sample_id[i]

    files = glob.glob(fastq_path)

    for f in files:
        for key in index_dic.keys():
            output = index_dic[key] + f.split("/")[-1].split(".fastq")[0].split("S")[1] + ".fastq"
            output_file = output_dir + output
            pbs_runners.script_runner("grep '%s' %s -A 3 > %s" % (key, f, output_file), alias="OST_Sample")


def trim(input_file, output_dir):
    new_f = output_dir + os.path.basename(input_file)
    new_text = ""
    with open(input_file) as infile:
        for l in infile:
            if l != "--\n":
                new_text += l
    out = open(new_f, "w")
    out.write(new_text)
    out.close()


def merge_paired_files(input_dir, out_dir):
    """
    :param input_dir: tmp directory path of the cirseq pipeline analysis
    :param out_dir: the output directory path
    :return: merged files
    """

    # print(input_dir)

    files1 = glob.glob(input_dir + "/*_*_*_R1_001.fastq")
    # print(files1)
    lst_files = []
    for fastq1 in files1:
        filename = os.path.basename(fastq1)
        # print(filename)

        #NextSeq
        # lane = filename.split("_")[1]

        sample = filename.split("_")[0] + "_" + filename.split("_")[1] + "_" + filename.split("_")[2]
        fastq2 = ("%s/%s_R2_001.fastq") % (input_dir, sample)
        if not os.path.exists(out_dir):
            out_dir = os.system("mkdir %s" % out_dir)
        output_file = "%s/%s_merged.fastq" % (out_dir, sample)

        pbs_runners.script_runner(
            "python /sternadi/home/volume3/okushnir/SternLab/scripts/merge_fastq_files.py -f %s -e %s -o %s -r 60" % (
            fastq1, fastq2, output_file), alias="merge_RV")

def main():
    # parser = OptionParser("usage: %prog [options]")
    # parser.add_option("-i", "--input_file", dest="input_file", help="fastq file")
    # parser.add_option("-o", "--output_dir", dest="output_dir", help="output dir")
    # (options, args) = parser.parse_args()
    # file = options.input_file
    # output_dir = options.output_dir

    # 1st thing to do is to index the output files from the Nextseq
    # csv_file = "/Volumes/STERNADILABHOME$/volume3/okushnir/RNAseq/180725_M04473_0026_000000000-BT9GF/index.csv"
    # fastq_path = "/Volumes/STERNADILABHOME$/volume3/okushnir/RNAseq/180725_M04473_0026_000000000-BT9GF/"
    # output_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/RNAseq/180725_M04473_0026_000000000-BT9GF/indexed"
    # index(csv_file, fastq_path, output_dir)

    # 2nd is to clean the files from --
    # trim(file, output_dir)

    # input_dir = ("/sternadi/datasets/volume1/180503_OST_FINAL_03052018/indexed/")
    # output_dir = ("/sternadi/home/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/sheri_clean/")
    # files = glob.glob(input_dir + "*.fastq")
    # for f in files:
    #     trim(f, output_dir)

    #MiSeq
    # PrimerID
    # Option 1
    # /sternadi/home/volume3/okushnir/Cluster_Scripts/extract-primer-ids20201016_I.cmd
    # /sternadi/home/volume3/okushnir/Cluster_Scripts/extract-primer-ids20201016_II.cmd
    # CHANGE THE PARAMETERS

    # Option 2
    # /sternadi/home/volume3/okushnir/Cluster_Scripts/barcode.cmd

    # 3rd merge the files

    # one by one approach
    # sample = "CVB3_p2_L001-ds.7381dd36e2ce42768db0cc15d5c79a31"
    # input_dir = "/sternadi/home/volume3/okushnir/AccuNGS/190627_RV_CV/clean/" + sample
    # merge_output_dir = "/sternadi/home/volume3/okushnir/AccuNGS/190627_RV_CV/merged/"
    # merge_paired_files(input_dir, merge_output_dir)


    # automatic approach

    # input_dir = "/sternadi/home/volume3/okushnir/AccuNGS/190807_RV_p2_p10/clean/"
    # merge_input_dir = glob.glob(input_dir+"*")
    # for dir in merge_input_dir:
    #     target = dir.split("/")[-1]
    #     merge_output_dir = "/sternadi/home/volume3/okushnir/AccuNGS/190807_RV_p2_p10/merged/" + target
    #     print(merge_output_dir)
    #     try:
    #         os.mkdir(merge_output_dir)
    #     except OSError:
    #         print("Creation of the directory %s failed" % merge_output_dir)
    #     else:
    #         print("Successfully created the directory %s " % merge_output_dir)
    #     print("merging...")
    #     merge_paired_files(dir, merge_output_dir)


    # # # 4th run pipeline:
    input_dir_RV_p7 = "/sternadi/home/volume3/okushnir/AccuNGS/190217_RV_p7/merged"
    input_dir_RV_new = "/sternadi/home/volume3/okushnir/AccuNGS/190807_RV_p2_p10/merged"
    input_dir_patients = "/sternadi/home/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/patients"
    input_dir_RV = "/sternadi/home/volume3/okushnir/AccuNGS/190627_RV_CV/merged/RVB14"
    input_dir_CV = "/sternadi/home/volume3/okushnir/AccuNGS/190627_RV_CV/merged/CVB3"
    input_dir_ATCG1 = "/sternadi/home/volume3/okushnir/AccuNGS/190627_RV_CV/merged/ATCG1"
    input_dir_FLNA = "/sternadi/home/volume3/okushnir/AccuNGS/190627_RV_CV/merged/FLNA"

    # Clinical samples
    # make_reference_from_consensus.py

    cycle = "5"
    patient_lst = ("Patient_1", "Patient_4", "Patient_5", "Patient_9", "Patient_16", "Patient_17", "Patient_20")

    ref_rv_patient = "/sternadi/home/volume3/okushnir/ref/RVA/JX025555.fasta"
    ref_rv = "/sternadi/home/volume3/okushnir/ref/RVB14/HRVB14_from_pWR3.26_1-7212.fasta"
    ref_cv = "/sternadi/home/volume3/okushnir/ref/CVB3/CVB3_from_pT7CVB3_1-7399.fasta"
    ref_atcg1 = "/sternadi/home//volume3/okushnir/ref/human/NM_001614.5_ACTG1.fa"
    ref_flna1 = "/sternadi/home/volume3/okushnir/ref/human/FLNA.fasta"

    ref_rv_patient1_consensus = "/sternadi/home/volume3/okushnir/ref/RVA/patient-1_consenX%s.fasta" % cycle
    ref_rv_patient4_consensus = "/sternadi/home/volume3/okushnir/ref/RVA/patient-4_consenX%s.fasta" % cycle
    ref_rv_patient5_consensus = "/sternadi/home/volume3/okushnir/ref/RVA/patient-5_consenX%s.fasta" % cycle
    ref_rv_patient9_consensus = "/sternadi/home/volume3/okushnir/ref/RVA/patient-9_consenX%s.fasta" % cycle
    ref_rv_patient16_consensus = "/sternadi/home/volume3/okushnir/ref/RVA/patient-16_consenX%s.fasta" % cycle
    ref_rv_patient17_consensus = "/sternadi/home/volume3/okushnir/ref/RVA/patient-17_consenX%s.fasta" % cycle
    ref_rv_patient20_consensus = "/sternadi/home/volume3/okushnir/ref/RVA/patient-20_consenX%s.fasta" % cycle

    ref_dic = {"Patient_1": "/sternadi/home/volume3/okushnir/ref/RVA/patient-1_consenX%s.fasta" % cycle,
    "Patient_4": "/sternadi/home/volume3/okushnir/ref/RVA/patient-4_consenX%s.fasta" % cycle,
    "Patient_5": "/sternadi/home/volume3/okushnir/ref/RVA/patient-5_consenX%s.fasta" % cycle,
    "Patient_9": "/sternadi/home/volume3/okushnir/ref/RVA/patient-9_consenX%s.fasta" % cycle,
    "Patient_16": "/sternadi/home/volume3/okushnir/ref/RVA/patient-16_consenX%s.fasta" % cycle,
    "Patient_17": "/sternadi/home/volume3/okushnir/ref/RVA/patient-17_consenX%s.fasta" % cycle,
    "Patient_20": "/sternadi/home/volume3/okushnir/ref/RVA/patient-20_consenX%s.fasta" % cycle}

    for patient in patient_lst:
        folders = glob.glob(input_dir_patients + "/%s" % patient)
        # print(folders)

        for d in folders:
            output_dir = (d + "/20201017_q30_consensusX%s" % cycle)
            # print(output_dir)
            try:
                os.mkdir(output_dir)
            except OSError:
                print("Creation of the directory %s failed" % output_dir)
            else:
                print("Successfully created the directory %s " % output_dir)

            cmd = "python /sternadi/home/volume1/shared/SternLab/pipeline_runner.py -i %s -o %s -r %s -NGS_or_Cirseq 2 " \
                  "-rep 2  -q 30 -b 40" % (d, output_dir, checkKey(ref_dic, patient))
            pbs_runners.script_runner(cmd, alias="pipeline_d")

    #5th analyze the freqs

    #6th run variant_caller localy to check context mutations

    # using /Users/odedkushnir/Google Drive/Studies/PhD/Python_Scripts/variant_caller_github.py (Maoz script)
    # I added pval to each position according to gamma distribution fit to RNA-Control - args = sample control -c min_coverage -o outpu_file_path
    # for example:
    # /Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/capsid/Capsid_32_Ultra/20201012_q38
    # /Capsid-32-Ultra.freqs /Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/controls/
    # IVT_3_Control/20201012_q38/IVT-3-Control.freqs -c 5000 -o /Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/
    # 20201008RV-202329127/merged/capsid/Capsid_32_Ultra/20201012_q38/Capsid-32_UltravsControl.csv

    # Using /Users/odedkushnir/Google Drive/Studies/PhD/Python_Scripts/after_variant_caller.py
    # I merged the freqs files to the variant_caller output file, to create a dataframe with the pval and Prob
    # Using /Users/odedkushnir/Google Drive/Studies/PhD/Python_Scripts/Context_analysis_RV.py
    # I created all new q38_data_mutation.csv, q38_data_UpX_by_mutation.csv, q38_data_XpA_by_mutation.csv files with new merged files.

#     Context_analysis_X.py
#     AccuNGS_analysis new_analysis_X.py

if __name__ == "__main__":
    main()
