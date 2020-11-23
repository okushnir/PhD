import os
import numpy as np
import pandas as pd


# input_dir = "/Users/odedkushnir/Projects/fitness/CirSeq/PV/OPV/"
# mutation_lst = ["AG", "AG_adar", "AG_nonadar", "GA", "UC", "CU"]
# fits_input_dir = input_dir + "fits/input/p1-p7/"
#
# jnum = []
# for mutation in mutation_lst:
#     list = os.listdir(fits_input_dir + mutation) # dir is your directory path
#     number_files = len(list)
#     jnum.append(number_files)
# print(jnum)
# jnum = max(jnum)
# print(jnum)
# a1 = 1
# d = 6
# an = a1+d*(jnum-1)
#
# print(an)
# axes = np.array([[1, 2, 3], [4, 5, 6]], np.int32)
# print(axes)
# for i, j in range(axes):
#     print(axes[i, j])

# print(type(range))

# pdser = pd.Series(["3838, 3844", "3889, 3994", "3919, 3920, 3925, 3929, 3930", "3985, 4008", "4018, 4036, 4096", "4028, 4029" ,
#           "4050, 4051", "4064, 4070", "4080, 4090", "4081, 4084, 4091, 4096"])
# def find_haplotype(pdser):
# pdser = pd.Series(["6019, 6043, 6049, 6067, 6076, 6078, 6079, 6082, 6097", "6037, 6040, 6043, 6049", "6040, 6043, 6049",
#                    "6040, 6043, 6049, 6076, 6078, 6079, 6082, 6088, 6106", "6043, 6049"])
# for x in pdser:
#     y = x.split(", ")
#     if len(y)>2:
#         continue
#     result = pdser.str.findall(x)
# print(type(result))
# print(result)

# def compare_2_strings(string1, string2):
#     string1 = "5383, 5402, 5421, 5425"
#     string2 = "5383, 5402, 5425"
#
#     listA = string1.split(", ")
#     listB = string2.split(", ")
#
#     result = list(set(listA) & set(listB))
#     return result

# list1 = ['6043', '6049', "875"]
# list2 = ['6079']
# L=[list1, list2]
# len_l = 0
# max = 0
# for l in L:
#     len_l = len(l)
#     if len_l > max:
#         max = len_l
#         max_list = l


# print(grouped.first())


# print(data.to_string())
# lst = (1, 2, 3)
# for i in lst:
#     if i == 1:
#         print(i)
#         continue
#     elif i == 2:
#         print(i)
#         continue
#     elif i == 3:
#         print(i)
#         continue
# print("Done!")

# Access to the values of the column and the ability to remove what you want
# data_filter_ag_grouped["passage"] = data_filter_ag_grouped["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
def main():
    all_mappings = pd.read_csv("/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/patients/"
                             "Patient_1/20201017_q30_consensusX5/all_parts.blast.cropped", names=["Read_id","Start","End"],
                             sep="\t")
    mutation = pd.read_csv("/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/patients/"
                           "Patient_1/20201017_q30_consensusX5/mutations_all.txt.cropped",
                           names=["Pos", "Read_id", "Mutant", "Read_positions"], sep="\t")
    freqs = pd.read_csv("/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/patients/"
                        "Patient_1/20201017_q30_consensusX5/Patient-1.freqs",
                           sep="\t")

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

    grouped.to_csv("/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/patients/Patient_1/"
                   "20201017_q30_consensusX5/grouped.csv", sep=",")


if __name__ == "__main__":
    main()
