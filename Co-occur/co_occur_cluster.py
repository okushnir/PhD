
#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import sys, argparse
from util import pbs_jobs

def blast_muation_creator(cmds, sample):
    cmd_file = "/sternadi/home/volume3/okushnir/Cluster_Scripts/all_parts.cmd"
    pbs_jobs.create_pbs_cmd(cmd_file, alias="all_parts", gmem=3, cmds=cmds, load_python=False)
    job_id = pbs_jobs.submit("-v sample='%s' %s" % (sample, cmd_file))
    print("-v sample='%s' %s" % (sample, cmd_file))
    print(job_id)
    status = pbs_jobs.check_pbs(job_id)
    if status == "Done":
        print("Done!")

def running_variants_on_the_same_read(sample, section_lst):
    for jnum in section_lst:
        cmds = "base=$sample\n" \
               "freqs=`ls ${base} | grep freqs`\n" \
               "mkdir ${base}/accungs_associations\n" \
               "python /sternadi/home/volume3/okushnir/PhD/variants_on_same_read.py ${base}/all_parts.blast.cropped ${base}/mutations_all.txt.cropped $PBS_ARRAY_INDEX ${base}/${freqs} > ${base}/accungs_associations/$PBS_ARRAY_INDEX.txt"
        cmd_file = "/sternadi/home/volume3/okushnir/Cluster_Scripts/co_occur.cmd"
        pbs_jobs.create_array_pbs_cmd(cmd_file, jnum, alias="accungs_assoc", gmem=3, cmds=cmds)
        print("qsub -v sample='%s' %s" % (sample, cmd_file))
        job_id = pbs_jobs.submit("-v sample='%s' %s" % (sample, cmd_file))
        # print(job_id)
        job_id = job_id.replace("[]", "")
        print(job_id)
        status = pbs_jobs.check_pbs(job_id)
        if status == "Done":
            print("Done %s" % (jnum))
    print("Done!!!!")

def concatenate_all_the_files(sample):
    cmds = "cd $sample/accungs_associations; cat *txt>all.txt"
    cmd_file = "/sternadi/home/volume3/okushnir/Cluster_Scripts/cat_txt.cmd"
    pbs_jobs.create_pbs_cmd(cmd_file, alias="cat_txt", gmem=3, cmds=cmds, load_python=False)
    print("qsub -v sample='%s' %s" % (sample, cmd_file))
    job_id = pbs_jobs.submit("-v sample='%s' %s" % (sample, cmd_file))
    print(job_id)
    status = pbs_jobs.check_pbs(job_id)
    if status == "Done":
        print("Done!")

def main(args):
    sample = args.sample
    blast_mutation = args.blast_mutation
    experiment = args.experiment

    """1. Create all_parts.blast, all_parts.blast.cropped, mutations_all.txt.cropped"""
    if blast_mutation == "True":
        print("Create all_parts.blast, all_parts.blast.cropped, mutations_all.txt.cropped")
        cmds_patients = "cd ${sample}; cat mutations_all.txt | grep -v ref_pos > " \
               "mutations_all.txt.cropped; rm all_parts.blast ; for file in `ls tmp/*.blast`; do cat $file >> all_parts.blast ; done ; " \
               "cat all_parts.blast | cut -f1,2,3 > all_parts.blast.cropped"
        blast_muation_creator(cmds_patients, sample)
    elif blast_mutation == "False":
        print("Running variants_on_same_read.py")
    """2. Run variants_on_same_read.py"""
    if (experiment == "passages") | (experiment == "capsid"):
        """Passages"""
        section_lst = ["3522-4999", "5000-6499", "6500-7212"]
    elif (experiment == "patients"):
        """Patients"""
        section_lst = ["336-1999", "2000-3335"]#


    running_variants_on_the_same_read(sample, section_lst)

    """3. Concatenate all the files"""

    concatenate_all_the_files(sample)

    """4. Run collect_cooccurs and merge it to freqs file"""
#     Run co_occur_local.py locally

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("blast_mutation", type=str, help="Create blast_all and mutaion_all files: True/False")
    parser.add_argument("sample", type=str, help="sample dir path")
    parser.add_argument("experiment", type=str, help="passages/capsid/patients")
    args = parser.parse_args(sys.argv[1:])
    main(args)

