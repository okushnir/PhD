#!/powerapps/share/python-anacondan-3.2019.7/bin/python

"""
@Author: odedkushnir

"""


import pandas as pd
import os
import numpy as np
import sys, argparse
import subprocess
from Context_analysis_RV import checkKey
import rpy2.robjects as robjects
import glob
from Context_analysis_RV import checkKey
from fits_plotter import *
from fits_parameters import *
from time import sleep
from fits_fitness_united import *


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
            sleep(0.05)
        print("")
    except (subprocess.CalledProcessError):
        process = ""
    if process == "":
        status = "Done"
    return status

def submit(cmdfile):
    cmd = "/opt/pbs/bin/qsub " + cmdfile
    result = os.popen(cmd).read()
    return result.split(".")[0]

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

    # Filtering
    # all["Frequency"] = np.where(all["pval"] > 0.01, 0, all["Frequency"])
    # all["Frequency"] = np.where(all["Prob"] < 0.95, 0, all["Frequency"])

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



def main(args):
    virus = args.virus
    passages = args.passages
    from_passage = int(passages.split("p")[1].split("-")[0])
    to_passage = int(passages.split("-p")[1])
    # input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/Cirseq/PV/OPV/39mixMOI/"
    input_dir = args.input_dir # directory contains q23_data_mutation.csv
    output_dir = input_dir + "fits/input/%s/" % passages
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    mutation_lst = ["AG", "AG_adar", "AG_nonadar", "GA", "UC", "CU"]

    # 1. Create fits dataset from data_mutation.csv file
    print("Creating fits dataset from data_mutation.csv file...")
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
    # break the dataset to positions datasets
    print("breaking the dataset to positions sub-datasets")
    for mutation in mutation_lst:
        df_mutation = pd.read_table(output_dir + "final_allMutations_sorted_%s.txt" % mutation)
        for position in df_mutation["pos"].iloc[0:-1]:
            df_pos = df_mutation.groupby(["pos"]).get_group(position)
            # print(df_pos.to_string())
            try:
                os.mkdir(output_dir + "%s" % mutation)
            except OSError:
                df_pos.to_csv(output_dir + "%s/final_allMutations_sorted_%s.txt" % (mutation, position), index=False,
                              sep="\t")

    # 2. Run FITS_jobarray_mutation.cmd

    output_dir = input_dir + "fits/output/mutation/%s/" % passages
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    job_id = submit("qsub /sternadi/home/volume3/okushnir/Cluster_Scripts/FITS_jobarray_mutation_%s.cmd" % virus)
    job_id = job_id.split("[")[0]
    print("Running FITS_jobarray_mutation.cmd, job_id:%s"%job_id)
    status = check_pbs(job_id)
    if status == "Done":
    # # 3. Run fits_parameters.py
        print("Creating fitness parameters")
        fitness_parameters(input_dir=input_dir + "fits/output/mutation/%s" % passages, output_dir=input_dir +
                                                                                                "fits/input/%s" % passages)

    # 4. Run FITS_jobarray_fitness.cmd
    output_dir = input_dir + "fits/output/fitness/%s/" % passages
    for mutation in mutation_lst:
        try:
            os.mkdir(output_dir+mutation)
        except OSError:
            print("Creation of the directory %s failed" % output_dir)
        else:
            print("Successfully created the directory %s " % output_dir)
    # print("qsub -v VIRUS='%s',PASSAGES='%s' /sternadi/home/volume3/okushnir/Cluster_Scripts/FITS_jobarray_fitness.cmd" % (virus, passages))
    job_id = submit("-v VIRUS='%s',PASSAGES='%s' /sternadi/home/volume3/okushnir/Cluster_Scripts/"
                    "FITS_jobarray_fitness.cmd" % (virus, passages))
    job_id = job_id.split("[")[0]
    print("Running FITS_jobarray_fitness.cmd, job_id:%s"%job_id)
    status = check_pbs(job_id)
    if status == "Done":
    # 5. Run SumSummaries_bi_fitness_OPV2
        for mutation in mutation_lst:
                fits_fitness_united(input_dir, output_file=output_dir + mutation + "all.txt")
            # subprocess.run("python /sternadi/home/volume3/okushnir/Python_scripts/SumSummaries_bi_fitness_OPV2.py "
            #                "/sternadi/home/volume3/okushnir/AccuNGS/190627_RV_CV/merged/%s/fits/output"
            #                "/fitness/%s/%s > "
            #                "/sternadi/home/volume3/okushnir/AccuNGS/190627_RV_CV/merged/%s/fits/output"
            #                            "/fitness/%s/%s/all.txt" % (virus, passages, mutation, virus, passages, mutation), shell=True)
            # print("python /sternadi/home/volume3/okushnir/Python_scripts/SumSummaries_bi_fitness_OPV2.py "
            #                "/sternadi/home/volume3/okushnir/AccuNGS/190627_RV_CV/merged/%s/fits/output"
            #                "/fitness/%s/%s > "
            #                "/sternadi/home/volume3/okushnir/AccuNGS/190627_RV_CV/merged/%s/fits/output"
            #                "/fitness/%s/%s/all.txt" % (virus, passages, mutation, virus, passages, mutation))
    # 6. Run fits_plotter.py

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("virus", type=str, help="name of the virus, RVB14")
    parser.add_argument("passages", type=str, help="from which passages, p0-p12")
    parser.add_argument("input_dir", "--input_dir", type=str, help="the path to the directory that contains q23_data_mutation.csv")

    args = parser.parse_args(sys.argv[1:])
    main(args)