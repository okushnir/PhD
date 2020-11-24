
#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import sys, argparse
from util import pbs_jobs



def blast_muation_creator(cmds):
    cmd_file = "/sternadi/home/volume3/okushnir/Cluster_Scripts/all_parts.cmd"
    pbs_jobs.create_pbs_cmd(cmd_file, alias="all_parts", gmem=3, cmds=cmds, load_python=False)
    job_id = pbs_jobs.submit(cmd_file)
    status = pbs_jobs.check_pbs(job_id)
    if status == "Done":
        print("Done!")

def running_variants_on_the_same_read(sample, section_lst):
    for jnum in section_lst:
        cmds = "base=$sample\n" \
               "freqs=`ls ${base} | grep freqs`\n" \
               "mkdir ${base}/accungs_associations\n" \
               "python /sternadi/home/volume1/maozgelbart/variants_on_same_read.py ${base}/all_parts.blast.cropped ${base}/mutations_all.txt.cropped $PBS_ARRAY_INDEX ${base}/${freqs} > ${base}/accungs_associations/$PBS_ARRAY_INDEX.txt"
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
    job_id = pbs_jobs.submit("-v sample='%s' %s" % (sample, cmd_file))
    print(job_id)
    status = pbs_jobs.check_pbs(job_id)
    if status == "Done":
        print("Done!")

def main(args):
    sample = args.sample
    blast_mutation = args.blast_mutation

    """1. Create all_parts.blast, all_parts.blast.cropped, mutations_all.txt.cropped"""
    if blast_mutation == True:
        """Patients"""
        cmds_patients = "for sample in 1 4 5 9 16 17 20; do cd /sternadi/home/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/" \
               "patients/Patient_${sample}/20201124_q30_consensusX7/; cat mutations_all.txt | grep -v ref_pos > " \
               "mutations_all.txt.cropped ; for file in `ls tmp/*.blast`; do cat $file >> all_parts.blast ; done ; " \
               "cat all_parts.blast | cut -f1,2,3 > all_parts.blast.cropped ; done"
        blast_muation_creator(cmds_patients)
    else:
        print("Running variants_on_same_read.py")
    """2. Run variants_on_same_read.py"""
    """Patients"""
    section_lst_patients = ["336-1999", "2000-3335"]

    running_variants_on_the_same_read(sample, section_lst_patients)

    """3. Concatenate all the files"""

    concatenate_all_the_files(sample)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("blast_mutation", type=bool, help="Create blast_all and mutaion_all files: True/False")
    parser.add_argument("sample", type=str, help="sample dir path")
    args = parser.parse_args(sys.argv[1:])
    main(args)

    """4. Run collect_cooccurs and merge it to freqs file"""
#     Run co_occur_patients.py locally


