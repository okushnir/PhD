
import pandas as pd
import numpy as np
import glob


def mutation_table(mutation_rate):
    with open(mutation_rate, "r") as handle:
        file = handle.readlines()
        tmp_lst = file[15].split("    ")
        tmp_lst = tmp_lst[2:7]
        data = pd.DataFrame(np.array(tmp_lst).reshape(-1, len(tmp_lst)))
        data = data.rename(columns={0: "median", 1: "MAD", 2: "min", 3:"max", 4: "pval"})
        data["Mutation"] = mutation_rate.split("/")[-1].split("_")[-1].split(".")[0]

    return data

def main():
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/" \
                    "fits_all_pos_at_once_sampling/replica2_syn/output/mutation/p2-p12/"
    file_lst = glob.glob(input_dir + "summary_mutation_syn_*.txt")
    columns = ["median", "MAD", "min", "max", "pval", "Mutation"]
    all_data = pd.DataFrame(columns=columns)
    for file in file_lst:
        table = mutation_table(file)
        all_data = all_data.append(table, ignore_index=True)
    all_data["Mutation"] = np.where(all_data["Mutation"] == "nonadar", "A>G non_adar", all_data["Mutation"])
    all_data["Mutation"] = np.where(all_data["Mutation"] == "adar", "A>G adar", all_data["Mutation"])
    all_data["Mutation"] = np.where(all_data["Mutation"] == "UC", "U>C", all_data["Mutation"])
    all_data["Mutation"] = np.where(all_data["Mutation"] == "CU", "C>U", all_data["Mutation"])
    all_data["Mutation"] = np.where(all_data["Mutation"] == "GA", "G>A", all_data["Mutation"])
    all_data["Mutation"] = np.where(all_data["Mutation"] == "AG", "A>G", all_data["Mutation"])
    all_data.to_csv(input_dir + "mutation_rate.csv", index=False)

if __name__ == "__main__":
    main()