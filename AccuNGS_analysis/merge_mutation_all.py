import os
import numpy as np
import pandas as pd


def merge_mutation_all(all_mappings, mutation, freqs, output_dir):
    freqs = freqs.loc[(freqs.Base != "-") & (freqs.Ref != "-")]
    mutation = mutation.rename(columns={"Mutant": "Base"})
    mutation = mutation.merge(freqs[["Pos", "Base", "Ref", "Freq"]], on=("Pos", "Base"), how="right")
    # mutation = mutation.merge(all_mappings, on="Read_id", how="left")
    mutation[["Read_id"]] = mutation[["Read_id"]].fillna(value="remove")
    mutation = mutation[mutation["Read_id"] != "remove"]
    mutation["Base"] = mutation["Base"].astype(str)
    mutation["Read_positions"] = mutation["Read_positions"].astype(str)
    mutation["Pos"] = mutation["Pos"].astype(str)
    mutation["Read_id"] = mutation["Read_id"].astype(str)
    mutation["Ref"] = mutation["Ref"].astype(str)
    mutation["Freq"] = mutation["Freq"].astype(str)
    # mutation["Read_count"] = mutation["Read_count"].astype(str)
    # mutation["Start"] = mutation["Start"].astype(str)
    # mutation["End"] = mutation["End"].astype(str)

    grouped = mutation.groupby(['Read_id']).agg(lambda x: x.mean() if np.issubdtype(x.dtype, np.number)
    else ', '.join(x)).reset_index()
    # grouped = grouped.sort_values(by=["Pos"])
    grouped["filter"] = grouped["Pos"].str.contains(",")
    grouped = grouped[grouped["filter"] == True]

    grouped.to_csv(output_dir + "/grouped.csv", sep=",")
    return grouped

def main():
    all_blast = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/patients/Patient_4/20201124_q30_consensusX7/all_parts.blast.cropped"
    mutation_all = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/patients/Patient_4/20201124_q30_consensusX7/mutations_all.txt.cropped"
    freqs_file = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/patients/Patient_4/20201124_q30_consensusX7/Patient-4.freqs"
    output_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/patients/Patient_4/20201124_q30_consensusX7/"
    all_mappings = pd.read_csv(all_blast, names=["Read_id","Start","End"],sep="\t")
    mutation = pd.read_csv(mutation_all, names=["Pos", "Read_id", "Mutant", "Read_positions"], sep="\t")
    freqs = pd.read_csv(freqs_file, sep="\t")

    grouped_df = merge_mutation_all(all_mappings, mutation, freqs, output_dir)

if __name__ == "__main__":
    main()
