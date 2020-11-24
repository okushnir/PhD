
#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import sys, argparse
import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
from util import pbs_jobs
import seaborn as sns
import numpy as np


class comutations_bubble(object):
    def __init__(self, a, b, freq, pval):
        self.nodes = set([a, b])
        self.distances = {'_'.join([str(a), str(b)]): freq, '_'.join([str(b), str(a)]): freq}
        self.pvalues = {'_'.join([str(a), str(b)]): pval, '_'.join([str(b), str(a)]): pval}
        self.meandist = np.mean(list(self.distances.values()))

    def __len__(self):
        return len(self.nodes)

    def __repr__(self):
        return ':'.join([str(a) for a in [len(self), min(self.nodes), max(self.nodes)]])

    def can_merge(self, other, distance=10):
        return len(self.nodes.intersection(other.nodes)) > 0 and (
                    (other.meandist * distance >= self.meandist) and (other.meandist / distance <= self.meandist))

    def union(self, other):
        self.nodes = self.nodes.union(other.nodes)
        self.distances.update(other.distances)
        self.pvalues.update(other.pvalues)
        self.meandist = np.mean(list(self.distances.values()))


def load_file(path, label):
    df = pd.read_csv(path, sep="\t", names=["start", "end", "fisher_stat", "pval", "variant_freq"])
    df["sample"] = label
    df = df[df["pval"] < 1]
    return df


def obtain_comutations(comutations, max_pval=10 ** -9, distance=10):
    dfs = []
    sig_positions = comutations[(comutations["pval"] < max_pval)]

    lines = sig_positions.itertuples(index=False)
    nodes = []
    for line in lines:
        a = line[0]
        b = line[1]
        pval = line[3]
        dist = line[4]
        node = comutations_bubble(a, b, dist, pval)
        nodes.append(node)

    results = []
    nodes = sorted(nodes, key=lambda item: -item.meandist)
    while len(nodes) > 0:
        bubble = nodes.pop(0)
        merged = False
        for item in nodes:
            if item.can_merge(bubble, distance):
                item.union(bubble)
                merged = True
                break
        if not merged:
            results.append(bubble)

    if results:
        # for cluster in sorted(results, key=lambda item: -item.meandist):
        #    print cluster.meandist, sorted(cluster.nodes)
        # Collate lines
        for i, cluster in zip(range(len(results)), sorted(results, key=lambda item: -item.meandist)):
            # distances=[]
            distances = map(lambda x: (int(x.split("_")[0]), int(x.split("_")[1]), cluster.distances[x]),
                            cluster.distances.keys())
            data = pd.DataFrame(distances, columns=["Pos1", "Pos2", "Freq"])
            data["Stretch"] = i

            data["Sample"] = comutations.iloc[0]['sample']
            data["meandist"] = cluster.meandist
            dfs.append(data)

    return pd.concat(dfs)


def collect_cooccurs(freqs_df, comutations_df, max_pval=10 ** -9, distance=10, add_ctx="ADAR"):
        sample_comutations = comutations_df
        sample_comutations = sample_comutations[sample_comutations["pval"] < 0.1]
        sample_comutations = obtain_comutations(sample_comutations, max_pval=max_pval, distance=distance)
        sample_comutations = sample_comutations[["Pos1", "Stretch", "meandist"]].sort_values(
            ["Pos1", "Stretch"]).drop_duplicates(['Pos1'])
        merged = freqs_df[freqs_df["Rank"] != 0]
        merged = merged.merge(sample_comutations, how="left", left_on="Pos", right_on="Pos1")

        merged["Stretch"] = "s" + merged["Stretch"].astype(str)
        merged["Stretch"] = merged["Stretch"].apply(lambda x: x.replace(".0", ""))
        merged["Stretch"] = np.where(merged["Stretch"] == "snan", '-', merged["Stretch"])
        merged["Co-occurrences_identified"] = np.where(merged["Stretch"] == "-", "No", "Yes")
        merged["ADAR_context"] = np.where(merged["Prev"].isin(["AA", "UA"]), "Yes",
                                          np.where(merged["Ref"] == "A", "No", "Not A>G"))
        merged["ADAR_reverse_context"] = np.where(merged["Next"].isin(["UU", "UA"]), "Yes",
                                                  np.where(merged["Ref"] == "U", "No", "Not U>C"))
        merged["APOBEC3F_context"] = np.where(merged["Next"].isin(["GA"]), "Yes",
                                              np.where(merged["Ref"] == "G", "No", "Not G>A"))
        # Avoid APOBEC3G
        # merged["APOBEC3G context"] = np.where(merged["next"].isin(["GG"]), "Yes",np.where(merged["Ref"]=="G","No","Not G>A"))
        merged["APOBEC3G_context"] = "No"
        merged["Editing_context"] = np.where(
            merged["ADAR_context"] == "Yes",
            "ADAR (sense)",
            np.where(
                merged["ADAR_reverse_context"] == "Yes",
                "ADAR (antisense)",
                np.where(
                    merged["APOBEC3F_context"] == "Yes",
                    "APOBEC3F",
                    np.where(
                        merged["APOBEC3G_context"] == "Yes",
                        "APOBEC3G",
                        "No editing context",
                    ),
                ),
            ),
        )

        if add_ctx == "ADAR":
            merged["ADAR"] = np.where(merged["Next"].isin(["UU", "UA"]),
                                      "Reverse complement",
                                      np.where(merged["Prev"].isin(["AA", "UA"]),
                                               "Forward",
                                               "No"
                                               )
                                      )
        return merged


def blast_muation_creator():
    cmds = "for sample in 1 4 5 9 16 17 20; do cd /sternadi/home/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/" \
           "patients/Patient_${sample}/20201124_q30_consensusX7/; cat mutations_all.txt | grep -v ref_pos > " \
           "mutations_all.txt.cropped ; for file in `ls tmp/*.blast`; do cat $file >> all_parts.blast ; done ; " \
           "cat all_parts.blast | cut -f1,2,3 > all_parts.blast.cropped ; done"
    cmd_file = "/sternadi/home/volume3/okushnir/Cluster_Scripts/all_parts.cmd"
    pbs_jobs.create_pbs_cmd(cmd_file, alias="all_parts", gmem=3, cmds=cmds, load_python=False)
    job_id = pbs_jobs.submit(cmd_file)
    status = pbs_jobs.check_pbs(job_id)
    if status == "Done":
        print("Done!")

def main(args):
    sample = args.sample
    blast_mutation = args.blast_muation

    """1. Create all_parts.blast, all_parts.blast.cropped, mutations_all.txt.cropped"""
    if blast_mutation == True:
        blast_muation_creator()
    else:
        continue
    """2. Run variants_on_same_read.py"""
    section_lst = ["336-1999", "2000-3335"] #
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
    """3. Concatenate all the files"""
    cmds = "cd $sample/accungs_associations; cat *txt>all.txt"
    cmd_file = "/sternadi/home/volume3/okushnir/Cluster_Scripts/cat_txt.cmd"
    pbs_jobs.create_pbs_cmd(cmd_file, alias="cat_txt", gmem=3, cmds=cmds, load_python=False)
    job_id = pbs_jobs.submit("-v sample='%s' %s" % (sample, cmd_file))
    print(job_id)
    status = pbs_jobs.check_pbs(job_id)
    if status == "Done":
        print("Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("blast_mutation", type=bool, help="Create blast_all and mutaion_all files: True/False")
    parser.add_argument("sample", type=str, help="sample dir path")
    args = parser.parse_args(sys.argv[1:])
    main(args)

    """4. Run collect_cooccurs and merge it to freqs file"""
#     Run co_occur_patients.py locally


