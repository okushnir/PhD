#!/powerapps/share/python-anaconda-3.2019.7/bin/python

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
from pbs_runners import array_script_runner
from pbs_jobs import create_array_pbs_cmd


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

def fits_data_construction(input_dir, output_dir, from_passage, to_passage, quality, start_position):
    all = pd.read_csv(input_dir + "%s_data_mutation.csv" % (quality))
    all = all.loc[all.Pos >= start_position]
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
    # all["pval"] = all["pval"].fillna(1)
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
    pipline_quality = args.quality
    output_dir = input_dir + "fits/input/"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    output_dir += "%s/" % passages
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    mutation_lst = ["AG", "AG_adar", "AG_nonadar", "GA", "UC", "CU"]

    # 1. Create fits dataset from data_mutation.csv file
    print("Creating fits dataset from data_mutation.csv file...")
    start_pos_dict = {"OPV": 3832, "RVB14": 3635, "CVB3": 3745} # start from 2B for all viruses
    start_position = checkKey(start_pos_dict, virus)
    fits_data_construction(input_dir,output_dir, from_passage, to_passage, pipline_quality, start_position)
    ##TODO: python script for this function
    # R-script modify 'virus', 'passage', and 'setwd' as previous
    test = robjects.r(''' library (plyr)
                        library(stringr)
                        virus = "OPV" ### MODIFY
                        passage = "p1-p7" ### MODIFY
                        setwd <- setwd(sprintf("/sternadi/home/volume3/okushnir/Cirseq/PV/%s/39mixMOI/fits/input/%s", virus, passage)) ### MODIFY
                        setwd <- setwd(sprintf("/sternadi/home/volume3/okushnir/Cirseq/PV/%s/39mixMOI/fits/input/%s", virus, passage)) ### MODIFY

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
    fits_input_dir = input_dir + "fits/input/%s/" % passages
    with open(fits_input_dir+"parameters_mutation.txt", "w") as parameter_mutation:
        parameter_mutation.write("# basic parameters\n"
                                 "N 1000000\n"
                                 "num_alleles 2\n"
                                 "num_generations 15\n"
                                 "min_log_mutation_rate0_0 -7\n"
                                 "min_log_mutation_rate0_1 -7\n"
                                 "min_log_mutation_rate1_0 -7\n"
                                 "min_log_mutation_rate1_1 -7\n"
                                 "max_log_mutation_rate0_0 -2\n"
                                 "max_log_mutation_rate0_1 -2\n"
                                 "max_log_mutation_rate1_0 -2\n"
                                 "max_log_mutation_rate1_1 -2\n"
                                 "# bottleneck\n"
                                 "bottleneck_interval 2\n"
                                 "bottleneck_size 1000000\n"
                                 "fitness_allele0 1.0\n"
                                 "fitness_allele1 1.0\n"
                                 "num_samples_from_prior 100000\n"
                                 "acceptance_rate 0.01")
    os.mkdir(input_dir + "fits/output/")
    os.mkdir(input_dir + "fits/output/mutation/")
    fits_mutation_output_dir = input_dir + "fits/output/mutation/%s/" % passages
    try:
        os.mkdir(fits_mutation_output_dir)
    except OSError:
        print("Creation of the directory %s failed" % fits_mutation_output_dir)
    else:
        print("Successfully created the directory %s " % fits_mutation_output_dir)

    job_id = submit("/sternadi/home/volume3/okushnir/Cluster_Scripts/FITS_jobarray_mutation_%s.cmd" % virus)
    job_id = job_id.split("[")[0]
    print("Running FITS_jobarray_mutation.cmd, job_id:%s"%job_id)
    status = check_pbs(job_id)
    if status == "Done":
    # 3. Run fits_parameters.py
        print("Creating fitness parameters")
        fitness_parameters(input_dir=fits_mutation_output_dir, output_dir=fits_input_dir)

    # 4. Run new_FITS_jobarray_fitness.cmd
    os.mkdir(input_dir + "fits/output/fitness/")
    cmds = "module load gcc/gcc-8.2.0\n" \
            "VIRUS=$VIRUS\n" \
           "PASSAGES=$PASSAGES\n" \
           "fits_dir=$fits_dir\n" \
           "PARAM='fitness'\n" \
           "namesj=('AG' 'UC' 'GA' 'CU' 'AG_adar' 'AG_nonadar')\n" \
           "i=$(($[PBS_ARRAY_INDEX-1-j]/${#namesj[@]}))\n" \
           "element=$(ls ${fits_dir}/input/${PASSAGES}/${namesj[j]}|grep -oh '[0-9]'*)\n" \
           "namesi=($element)\n" \
           "for element in '${namesi[@]}';do echo '$element';done\n" \
           "mkdir ${fits_dir}/output/${PARAM}/${PASSAGES}\n" \
           "mkdir ${fits_dir}/output/${PARAM}/${PASSAGES}/${namesj[j]}\n" \
           "/sternadi/home/volume1/talzinger/FITS_Analyses/FITS_bin/fits1.3.3 -$PARAM ${fits_dir}/input/${PASSAGES}/" \
            "parameters_${PARAM}_${namesj[j]}.txt ${fits_dir}/input/${PASSAGES}/${namesj[j]}/final_allMutations_sorted_" \
            "${namesi[i]}.txt${fits_dir}/output/${PARAM}/${PASSAGES}/${namesj[j]}/posterior_${PARAM}_all_${namesi[i]}.txt " \
            "${fits_dir}/output/${PARAM}/${PASSAGES}/${namesj[j]}/summary_${PARAM}_all_${namesi[i]}.txt"
    jnum = 0

    for mutation in mutation_lst:
        list = os.listdir(fits_input_dir + mutation) # dir is your directory path
        number_files = len(list)
        jnum += number_files
    jnum += 1000
    print(jnum)
    cmd_file = "/sternadi/home/volume3/okushnir/Cluster_Scripts/new_FITS_jobarray_fitness_%s.cmd" % virus
    create_array_pbs_cmd(cmd_file, jnum, alias="Fits_fitness", gmem=3, cmds=cmds)

    fitness_output_dir = input_dir + "fits/output/fitness/%s/" % passages
    try:
        os.mkdir(fitness_output_dir)
    except OSError:
        print("Creation of the directory %s failed" % fitness_output_dir)
    else:
        print("Successfully created the directory %s " % fitness_output_dir)
    for mutation in mutation_lst:
        try:
            os.mkdir(fitness_output_dir+mutation)
        except OSError:
            print("Creation of the directory %s failed" % fitness_output_dir+mutation)
        else:
            print("Successfully created the directory %s " % output_dir+mutation)
    # print("qsub -v VIRUS='%s',PASSAGES='%s' /sternadi/home/volume3/okushnir/Cluster_Scripts/FITS_jobarray_fitness.cmd" % (virus, passages))
    fits_dir = input_dir + "fits"
    # job_id = submit("-v VIRUS='%s',PASSAGES='%s, fits_dir=%s' /sternadi/home/volume3/okushnir/Cluster_Scripts/"
    #                 "FITS_jobarray_fitness_%s.cmd" % (virus, passages, fits_dir, virus))
    print("qsub -v VIRUS='%s',PASSAGES='%s', fits_dir='%s' %s" % (virus, passages, fits_dir, cmd_file))
    job_id = submit("-v VIRUS='%s',PASSAGES='%s, fits_dir=%s' %s" % (virus, passages, fits_dir, cmd_file))
    job_id = job_id.split("[")[0]
    print("Running new_FITS_jobarray_fitness_%s.cmd, job_id:%s" % (virus, job_id))
    status = check_pbs(job_id)
    if status == "Done":
    # 5. Run fits_fitness_united
        for mutation in mutation_lst:
            output_fitness_dir = input_dir + "fits/output/fitness/%s/%s/" % (passages, mutation)
            output_file = output_fitness_dir + "all.txt"
            print("Creating the fitness conjugated report of: %s" % output_fitness_dir)
            fits_fitness_united(output_fitness_dir, output_file)
    # 6. Run fits_plotter.py

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("virus", type=str, help="name of the virus, RVB14")
    parser.add_argument("passages", type=str, help="from which passages, p0-p12")
    parser.add_argument("input_dir", type=str, help="the path to the directory that contains q23_data_mutation.csv")
    parser.add_argument("quality", type=str, help="what is the prefix for the data_mutation.csv file; quality of the pipline ; for example: q38")
    args = parser.parse_args(sys.argv[1:])
    main(args)