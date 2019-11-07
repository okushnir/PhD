"""
@Author: odedkushnir

"""

import sys
import glob
import pandas as pd
import glob
import shutil
import os
from os import path
import re



def check_filename(filename, Truefile = True):
    """
    checks if filename is legit and returns its absolute path
    :param filename: input file path
    :param Truefile: is it an existing file (default: True)
    :return: absolute filename path
    """
    if filename == None:
        raise Exception("you must specify a file name")
    filename = path.abspath(filename)
    if Truefile and not path.isfile(filename):
        raise Exception("file %s is not a real file" % filename)
    return filename


def check_dirname(dirname, Truedir = True):
    """
    checks if dirname is legit and returns its absolute path
    :param dirname: input directory path
    :param Truedir: is it an existing dir (default: True)
    :return: absolute directory path
    """
    if dirname == None:
        raise Exception("you must specify a dir name")
    if not path.isabs(dirname):
        dirname = path.abspath(dirname)
    if Truedir and not path.isdir(dirname):
        raise Exception("dir name %s does not exist" % dirname)
    return dirname



def copy_talia_files(in_dir, out_dir):
    # in_dir = "/Volumes/STERNADILABTEMP$/volume1/phyVirus"
    # out_dir = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/RVB14/InVivo"

    #text_file = open("/Users/odedkushnir/Google Drive/Studies/PhD/Projects/RVB14/InVivo/rhinovirus_files.txt", "r")
    data = pd.read_csv("/Users/odedkushnir/Google Drive/Studies/PhD/Projects/RVB14/InVivo/rhinovirus_files.txt", sep="\n")
    data.rename(columns={"Files_names\n": "Files_name"}, inplace=True)
    data["New"] = data.Files_name.str[0:-2]
    file_lst = list(data["New"])

    # file_lst_new = []
    # #a = pd.DataFrame()
    # alignments_files = glob.glob(in_dir + "/final_prank/*")
    # for file in alignments_files:
    #     file_name = file.split("/")[-1].split(".")[0] + "." + file.split("/")[-1].split(".")[1] + "."
    #     # print(file_name)
    #     for i in file_lst:
    #         if i == file_name:
    #             file_lst_new.append(i)
    #             if not os.path.exists(out_dir + "/final_prank"):
    #                 os.mkdir(out_dir + "/final_prank")
    #             shutil.copyfile(in_dir + "/final_prank/" + i + "aln.best.phy", out_dir + "/final_prank/" + i + "aln.best.phy")
    #             shutil.copyfile(in_dir + "/final_prank/" + i + "aln.best.fas", out_dir + "/final_prank/" + i + "aln.best.fas")
    #
    # tree_files = glob.glob(in_dir + "/final_phyml/*")
    # for tree in tree_files:
    #     file_name = tree.split("/")[-1].split(".")[0] + "." + tree.split("/")[-1].split(".")[1] + "."
    #     # print(file_name)
    #     for i in file_lst_new:
    #         if i == file_name:
    #             if not os.path.exists(out_dir + "/final_phyml"):
    #                 os.mkdir(out_dir + "/final_phyml")
    #             shutil.copyfile(in_dir + "/final_phyml/" + i + "aln.best.phy_phyml_stats.txt", out_dir + "/final_phyml/" + i + "aln.best.phy_phyml_stats.txt")
    #             shutil.copyfile(in_dir + "/final_phyml/" + i + "aln.best.phy_phyml_tree.txt", out_dir + "/final_phyml/" + i + "aln.best.phy_phyml_tree.txt")
    #
    #
    # gi_files = glob.glob(in_dir + "/final_gi/*")
    # for gi in gi_files:
    #     file_name = gi.split("/")[-1].split(".")[0] + "." + gi.split("/")[-1].split(".")[1] + "."
    #     # print(file_name)
    #     for i in file_lst_new:
    #         if i == file_name:
    #             if not os.path.exists(out_dir + "/final_gi"):
    #                 os.mkdir(out_dir + "/final_gi")
    #             shutil.copyfile(in_dir + "/final_gi/" + i + "gi", out_dir + "/final_gi/" + i + "gi")
    #             shutil.copyfile(in_dir + "/final_gi/" + i + "gb", out_dir + "/final_gi/" + i + "gb")

    baseml_files = glob.glob(in_dir + "/*")
    for file in baseml_files:
        file_name = file.split("/")[-1].split(".")[0] + "." + file.split("/")[-1].split(".")[1] + "."
        # print(file_name)
        for i in file_lst:
            if i == file_name:
                #file_lst_new.append(i)
                if not os.path.exists(out_dir + "/baseml_model8"):
                    os.mkdir(out_dir + "/baseml_model8")
                shutil.copyfile(in_dir + i + "mlb8",
                                out_dir + "/baseml_model8/" + i + "mlb8")
                shutil.copyfile(in_dir + i + "ctl",
                                out_dir + "/baseml_model8/" + i + "ctl")


HYPHY_PROGRAM = "/Users/odedkushnir/Projects/InVivo/hyphy-2.3.13/bin/HYPHYMP"
HYPHY_BF = "/Users/odedkushnir/Projects/InVivo/hyphy-2.3.13/res/TemplateBatchFiles/AnalyzeDiNucData_talia.bf"


def run_hyphy(out_dir, model, hyphy_program=HYPHY_PROGRAM, hyphy_bf=HYPHY_BF):
    aln_dir = glob.glob(out_dir + "/final_prank/*.fas")
    for aln in aln_dir:
        base = aln.split("/")[-1].split("aln.best.fas")[0]
        tree = glob.glob(out_dir + "/final_phyml/" + base + "*tree*")[0]
        # tree = out_dir + "/final_phyml" + base
        output = out_dir + "/AG_results/" + base + "hyphy.txt"

        if not os.path.isfile(output):
            print(aln)
        os.system("%s %s %s %s %s %s" % (hyphy_program, hyphy_bf, aln, tree, output, model))



cols= ['R_TGTT', 'R_TCTT', 'R_AAAC', 'R_TCTG', 'R_CTTT',
       'R_CTGT', 'R_GAGC', 'R_GAGT', 'R_GAGG', 'R_CGTG',
       'R_CCTC', 'R_CCGC', 'R_CGCT', 'R_CGGG', 'R_GATA',
       'R_TATC', 'R_GTTT', 'R_TATG', 'R_TATT', 'R_GGTG',
       'R_GCGT', 'R_GCGG', 'R_GCTC', 'R_GGGT', 'R_CCCT',
       'R_ACAT', 'R_ACAG', 'R_ACCC', 'R_ACGC', 'R_AATA',
       'R_AAAT', 'R_AAAG', 'R_AACA', 'R_AAGA', 'R_AGAT',
       'R_CAGA', 'R_CACT', 'R_CCCG', 'R_CATA', 'R_CACG',
       'R_AGGG', 'R_AGTG', 'R_AGCG', 'R_ATGT', 'R_ATTT',
       'R_ATCT', 'R_CACC', 'AIC', 'lnL', 'rate_Normalizer',
       'filename', 'family', 'protein', 'group']



def extract_hyphy_results(files=[]):
    df = pd.DataFrame(columns=cols)
    if files==[]:
        print("no files!!")
    for f in files:
         with open(f, "r") as file:
             family = f.split("/")[-1].split("_")[0]
             protein = f.split("/")[-1].split(family + "_")[1].split(".")[0]
             group = f.split("/")[-1].split(family + "_")[1].split(".hyphy.txt")[0]
             res = {}
             res["filename"] = f
             res["protein"] = protein
             res["group"] = group
             res["family"] = family
             data = file.readlines()
             for l in data:
                     if "AIC" in l:
                             aic = float(l.split("=")[-1].strip().split("\n")[0])
                             res["AIC"] = aic
                     elif "Log Likelihood" in l:
                         lnL = float(l.split("=")[-1].strip().split("\n")[0].split(";")[0])
                         res["lnL"] = lnL
                     elif "rate_Normalizer" in l:
                         rate_Normalizer = float(l.split("=")[-1].strip().split("\n")[0])
                         res["rate_Normalizer"] = rate_Normalizer
                     elif "R_" in l:
                         name = l.split("=")[0]
                         value = float(l.split("=")[-1].strip().split("\n")[0])
                         res[name] = value
             df = df.append(res, ignore_index=True)
             df.to_csv("/Users/odedkushnir/Projects/InVivo/AG_results/results.csv", sep=",", encoding='utf-8')


def mlbs_to_df(output, mlbs=[], dirname=None):
    """
    analyzes mlb file to dataframe - extracts lnL, base frequencies and substiution matrics
    :param output: output csv file path
    :param mlbs: list of mlb files
    :param dirname: dirname that has mlb files
    :return: output file path
    """
    if mlbs == [] and dirname == None:
        raise Exception("you need to provide mlb or dirname that contains mlbs")
    if mlbs != [] and dirname != None:
        raise Exception("you need to provide only one - mlb or dirname")

    if dirname != None:
        dirname = check_dirname(dirname)
        mlbs = glob.glob(dirname + "/*.mlb")
    if mlbs != []:
        mlbs = [check_filename(m) for m in mlbs]

    output = check_filename(output, Truefile=False)

    df = pd.DataFrame(columns=["mlb_file_name", "family", "group", "model", "lnL",
                               "freq_T", "freq_C", "freq_A", "freq_G",
                               "TC", "TA", "TG", "CT", "CA", "CG", "AT",
                               "AC", "AG", "GT", "GC", "GA"])

    lnL_1 = re.compile("lnL.*")
    lnL_2 = re.compile("\-?\d*\.\d*")
    base_1 = re.compile("Base frequencies.*")
    base_2 = re.compile("0.\d+")
    rate_1 = re.compile("Rate matrix Q.*\n.*\n.*\n.*\n.*", re.IGNORECASE)
    rate_2 = re.compile("\d+.\d+")
    for mlb_file_name in mlbs:
        print(mlb_file_name)
        family = mlb_file_name.split("/")[-2]
        filename = mlb_file_name.split("/")[-1]
        if "_gtr" in filename or "_unrest" in filename:
            filename = filename.split("_gtr")[0]
            filename = filename.split("_unrest")[0]
        model = mlb_file_name.split(".mlb")[0].split("_")[-1]

        mlb = open(mlb_file_name, "r").read()
        L = lnL_1.findall(mlb)
        if len(L) != 1:
            L = None
        elif "nan" in L[0]:
            L = None
        else:
            L = float(lnL_2.findall(L[0])[0])

        B = base_1.findall(mlb)
        if len(B) != 1:
            freq_T = None;
            freq_C = None
            freq_A = None;
            freq_G = None
        elif "nan" in B[0]:
            freq_T = None;
            freq_C = None
            freq_A = None;
            freq_G = None
        else:
            B = base_2.findall(B[0])
            freq_T = float(B[0])
            freq_C = float(B[1])
            freq_A = float(B[2])
            freq_G = float(B[3])

        R = rate_1.findall(mlb)
        if len(R) != 1:
            TC = None;
            TA = None;
            TG = None;
            CT = None;
            CA = None;
            CG = None;
            AT = None;
            AC = None;
            AG = None;
            GT = None;
            GC = None;
            GA = None
        elif len(R) >= 1 and "nan" in R[0]:
            TC = None;
            TA = None;
            TG = None;
            CT = None;
            CA = None;
            CG = None;
            AT = None;
            AC = None;
            AG = None;
            GT = None;
            GC = None;
            GA = None

        else:
            R = R[0].split("\n")
            first = R[1]
            first = rate_2.findall(first)
            TC = first[1];
            TA = first[2];
            TG = first[3]
            second = R[2]
            second = rate_2.findall(second)
            CT = second[0];
            CA = second[2];
            CG = second[3]
            third = R[3]
            third = rate_2.findall(third)
            AT = third[0];
            AC = third[1];
            AG = third[3]
            fourth = R[4]
            fourth = rate_2.findall(fourth)
            GT = fourth[0];
            GC = fourth[1];
            GA = fourth[2]

        df = df.append({"mlb_file_name": mlb_file_name, "family": family, "group": filename, "model": model, "lnL": L,
                        "freq_T": freq_T, "freq_C": freq_C, "freq_A": freq_A,
                        "freq_G": freq_G, "TC": TC, "TA": TA, "TG": TG,
                        "CT": CT, "CA": CA, "CG": CG, "AT": AT,
                        "AC": AC, "AG": AG, "GT": GT, "GC": GC, "GA": GA},
                       ignore_index=True)

    df.to_csv(output)
    return (output)


def main():
    # in_dir = "/Volumes/STERNADILABTEMP$/volume1/phyVirus/baseml_model8/"
    # out_dir = "/Users/odedkushnir/Projects/InVivo"
    # copy_talia_files(in_dir, out_dir)

    # out_dir = "/Users/odedkushnir/Projects/InVivo"
    # model = "/Users/odedkushnir/Projects/InVivo/hyphy-2.3.13/res/substitution_dinuc_context_AG"

    # run_hyphy(out_dir, model)

    # files = glob.glob("/Users/odedkushnir/Projects/InVivo/AG_results/*")
    # extract_hyphy_results(files)

    output = "/Users/odedkushnir/Projects/InVivo/baseml_model8/mlb_result.csv"
    dirname = "/Volumes/STERNADILABTEMP$/volume1/phyVirus/baseml_model8/"
    mlbs = glob.glob("/Users/odedkushnir/Projects/InVivo/baseml_model8/*.mlb8")
    mlbs_to_df("/Users/odedkushnir/Projects/InVivo/baseml_model8/mlb_result.csv", mlbs=mlbs, dirname=None)


if __name__ == "__main__":
    main()