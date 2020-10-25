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
lst = (1, 2, 3)
for i in lst:
    if i == 1:
        print(i)
        continue
    elif i == 2:
        print(i)
        continue
    elif i == 3:
        print(i)
        continue
print("Done!")