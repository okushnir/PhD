

"""
@Author: odedkushnir

"""

import pandas as pd
import re
import numpy as np

"""
With this script one can analyze RealTime PCR result from CFXConnect.
The samples names should be as follows:
X_BioReplicate No'_TechnicalReplicate No' -> HeLa_1_1
"""
def checkKey(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        raise Exception()

def load_data_to_df(sample_file1, sample_file2, sample_file3, gene):
    label_sample1 = gene + ".plate1"
    print("loading " + sample_file1 + " as sample")
    data_RT1 = pd.read_csv(sample_file1)
    data_RT1["Source"] = label_sample1

    label_sample2 = gene + ".plate2"
    print("loading " + sample_file2 + " as sample")
    data_RT2 = pd.read_csv(sample_file2)
    data_RT2["Source"] = label_sample2

    label_sample3 = gene + ".plate3"
    print("loading " + sample_file3 + " as sample")
    data_RT3 = pd.read_csv(sample_file3)
    data_RT3["Source"] = label_sample3

    # Organize the Data
    data = pd.concat([data_RT1, data_RT2, data_RT3])
    return data



def replica_analyse(data, gene):

    columns = ["Target", "Sample", "Biological Set Name", "Cq", "Source"]
    data_filter = pd.DataFrame(data, columns=columns)
    # print(data_filter)
    # Make Replica.no Column
    pattern = re.compile("(\d|\d.\d\d|\d.\d|\d\d)$", re.MULTILINE)

    data_filter["Bio.Replica.no"] = data_filter.Sample.str[-1]

    # Make Time Column
    data_filter["Sample"] = data_filter.Sample.str[0:-2]
    data_filter["Time"] = data_filter["Sample"].str.extract(pattern, expand=True)

    # Make long data to wide data
    data_filter_wide = (data_filter.groupby(['Sample', 'Target', 'Biological Set Name', 'Cq']).agg(
        lambda x: x.mean() if np.issubdtype(x.dtype, np.number) else ', '.join(x))).reset_index()
    data_filter_wide["Time"] = data_filter_wide.Sample.str.extract(pattern, expand=True)
    data_filter_wide_rep_1 = data_filter_wide[data_filter_wide["Bio.Replica.no"] == "1"]
    data_filter_wide_rep_2 = data_filter_wide[data_filter_wide["Bio.Replica.no"] == "2"]
    data_filter_wide_rep_3 = data_filter_wide[data_filter_wide["Bio.Replica.no"] == "3"]
    data_filter_wide_rep_1["Tech.Rep.no"] = data_filter_wide_rep_1.groupby(['Sample', 'Target']).cumcount()
    data_filter_wide_rep_2["Tech.Rep.no"] = data_filter_wide_rep_2.groupby(['Sample', 'Target']).cumcount()
    data_filter_wide_rep_3["Tech.Rep.no"] = data_filter_wide_rep_3.groupby(['Sample', 'Target']).cumcount()

    df_list = (data_filter_wide_rep_1, data_filter_wide_rep_2, data_filter_wide_rep_3)
    rep_df_lst = []
    for df in df_list:
    # RV-
        cal_RV_minus_df = df[df['Biological Set Name'] == "RV-"]
        cal_RV_minus_df.rename(columns={'Biological Set Name': 'Biological_Set_Name'}, inplace=True)
        # Separate the DataFrame to two, one to each target: ADAR & GAPDH
        cal_adar = cal_RV_minus_df[cal_RV_minus_df["Target"] == gene].reset_index()
        cal_gapdh = cal_RV_minus_df[cal_RV_minus_df["Target"] == "GAPDH"].reset_index()
        # Merge the DataFrame back together
        cal_RV_minus_df_final = pd.merge(cal_adar, cal_gapdh, on=("Time", "Bio.Replica.no", "Tech.Rep.no", "Sample", "Source"))
        cal_RV_minus_df_final["dCt.Cal"] = cal_RV_minus_df_final["Cq_x"] - cal_RV_minus_df_final["Cq_y"]

        # RV+
        test_RV_plus_df = df[df['Biological Set Name'] == "RV+"]
        test_RV_plus_df.rename(columns={'Biological Set Name': 'Biological_Set_Name'}, inplace=True)
        # Separate the DataFrame to two, one to each target: ADAR & GAPDH
        test_adar = test_RV_plus_df[test_RV_plus_df["Target"] == gene].reset_index()
        test_gapdh = test_RV_plus_df[test_RV_plus_df["Target"] == "GAPDH"].reset_index()
        # Merge the DataFrame back together
        test_RV_plus_df_final = pd.merge(test_adar, test_gapdh, on=("Time", "Bio.Replica.no", "Tech.Rep.no", "Sample", "Source"))
        test_RV_plus_df_final["dCt.Test"] = test_RV_plus_df_final["Cq_x"] - test_RV_plus_df_final["Cq_y"]
        # Merge the 2 DataFrame Cal and Test
        df = pd.concat([cal_RV_minus_df_final, test_RV_plus_df_final])
        rep_df_lst.append(df)
    return rep_df_lst


def just_ddCT_replica(dct_df):

    dct_df.rename(columns={'Biological_Set_Name_x': 'Biological_Set_Name'}, inplace=True)

    ddct_df = pd.DataFrame()
    ddct_df["Biological_Set_Name"] = dct_df["Biological_Set_Name"]
    ddct_df["Time"] = dct_df["Time"]
    ddct_df["dCt_Cal"] = dct_df["dCt.Cal"]
    ddct_df["dCt_Test"] = dct_df["dCt.Test"]
    ddct_df["Bio.Replica.no"] = dct_df["Bio.Replica.no"]
    ddct_df["Tech.Rep.no"] = dct_df["Tech.Rep.no"]

    # ddct_df = (ddct_df.groupby(['Biological_Set_Name', "Time"]).agg(lambda x: x.mean() if np.issubdtype(x.dtype, np.number) else ', '.join(x))).reset_index()
    dct_minus = ddct_df[ddct_df["Biological_Set_Name"] == "RV-"].reset_index()
    dct_plus = ddct_df[ddct_df["Biological_Set_Name"] == "RV+"].reset_index()

    ddct_final = pd.merge(dct_plus, dct_minus, on=("Time", "Tech.Rep.no", "Bio.Replica.no"))
    ddct_final["ddCt"] = ddct_final["dCt_Test_x"] - ddct_final["dCt_Cal_y"]
    ddct_final["Relative Normalized Expression"] = 2**(-ddct_final["ddCt"])
    return ddct_final


def dct_cal_test_calculator(sample_file1, sample_file2, sample_file3, gene):
    """
    :param sample_file1: Quantification Cq Results_0 csv file came from the RT-PCR machine i.e plate1
    :param sample_file2: Quantification Cq Results_0 csv file came from the RT-PCR machine i.e plate2
    :param gene: Analyzed gene name
    :return:
    """

    label_sample1 = gene + ".plate1"
    print("loading " + sample_file1 + " as sample")
    data_RT1 = pd.read_csv(sample_file1)
    data_RT1["Source"] = label_sample1

    label_sample2 = gene + ".plate2"
    print("loading " + sample_file2 + " as sample")
    data_RT2 = pd.read_csv(sample_file2)
    data_RT2["Source"] = label_sample2
    
    label_sample3 = gene + ".plate3"
    print("loading " + sample_file3 + " as sample")
    data_RT3 = pd.read_csv(sample_file3)
    data_RT3["Source"] = label_sample3


    # Organize the Data
    data = pd.concat([data_RT1, data_RT2, data_RT3])
    data.to_csv("/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/Results/all_3samples_raw.csv", sep=",", encoding='utf-8')
    columns = ["Target", "Sample", "Biological Set Name", "Cq Mean", "Source"]
    data_filter = pd.DataFrame(data, columns=columns)
    print(data_filter)

    # Make Replica.no Column
    pattern = re.compile("(\d|\d.\d\d|\d.\d|\d\d)$", re.MULTILINE)

    data_filter["Replica.no"] = data_filter.Sample.str[-1]

    # Make Time Column
    data_filter["Sample"] = data_filter.Sample.str[0:-2]
    data_filter["Time"] = data_filter["Sample"].str.extract(pattern, expand=True)

    # Make long data to wide data
    data_filter_wide = (data_filter.groupby(['Sample', 'Target', 'Biological Set Name']).agg(lambda x: x.mean() if np.issubdtype(x.dtype, np.number) else ', '.join(x))).reset_index()
    data_filter_wide["Time"] = data_filter_wide.Sample.str.extract(pattern, expand=True)
    print(data_filter_wide)

    # Separate DataFrame to two, one to each biological name set: RV- & RV+ ->
    # to create columns with the values instead of rows

    # RV-
    cal_RV_minus_df = data_filter_wide[data_filter_wide['Biological Set Name'] == "RV-"]
    cal_RV_minus_df.rename(columns={'Cq Mean': 'Cq_Mean'}, inplace=True)
    cal_RV_minus_df.rename(columns={'Biological Set Name': 'Biological_Set_Name'}, inplace=True)
    # Separate the DataFrame to two, one to each target: ADAR & GAPDH
    cal_adar = cal_RV_minus_df[cal_RV_minus_df["Target"] == gene].reset_index()
    cal_gapdh = cal_RV_minus_df[cal_RV_minus_df["Target"] == "GAPDH"].reset_index()
    # Merge the DataFrame back together
    cal_RV_minus_df_final = pd.merge(cal_adar, cal_gapdh, on="Time")
    cal_RV_minus_df_final["dCt.Cal"] = cal_RV_minus_df_final["Cq_Mean_x"] - cal_RV_minus_df_final["Cq_Mean_y"]

    # RV+
    test_RV_plus_df = data_filter_wide[data_filter_wide['Biological Set Name'] == "RV+"]
    test_RV_plus_df.rename(columns={'Cq Mean': 'Cq_Mean'}, inplace=True)
    test_RV_plus_df.rename(columns={'Biological Set Name': 'Biological_Set_Name'}, inplace=True)
    # Separate the DataFrame to two, one to each target: ADAR & GAPDH
    test_adar = test_RV_plus_df[test_RV_plus_df["Target"] == gene].reset_index()
    test_gapdh = test_RV_plus_df[test_RV_plus_df["Target"] == "GAPDH"].reset_index()
    # Merge the DataFrame back together
    test_RV_plus_df_final = pd.merge(test_adar, test_gapdh, on="Time")
    test_RV_plus_df_final["dCt.Test"] = test_RV_plus_df_final["Cq_Mean_x"] - test_RV_plus_df_final["Cq_Mean_y"]
    # Merge the 2 DataFrame Cal and Test
    dct_df = pd.concat([cal_RV_minus_df_final, test_RV_plus_df_final])
    return dct_df

def just_ddCT(dCT_file):
    """
    :param dCT_file: csv file from dct_cal_test_calculator function
    :return: csv file with ddct and FoldChange columns
    """
    dct_df = pd.read_csv(dCT_file)
    pattern = re.compile("(\d|\d.\d\d|\d.\d|\d\d)$", re.MULTILINE)
    dct_df["Time"] = dct_df["Sample_x"].str.extract(pattern, expand=True)
    dct_df.rename(columns={'Biological_Set_Name_x': 'Biological_Set_Name'}, inplace=True)

    ddct_df = pd.DataFrame()
    ddct_df["Biological_Set_Name"] = dct_df["Biological_Set_Name"]
    ddct_df["Time"] = dct_df["Time"]
    ddct_df["dCt_Cal"] = dct_df["dCt.Cal"]
    ddct_df["dCt_Test"] = dct_df["dCt.Test"]

    ddct_df = (ddct_df.groupby(['Biological_Set_Name', "Time"]).agg(lambda x: x.mean() if np.issubdtype(x.dtype, np.number) else ', '.join(x))).reset_index()
    dct_minus = ddct_df[ddct_df["Biological_Set_Name"] == "RV-"].reset_index()
    dct_plus = ddct_df[ddct_df["Biological_Set_Name"] == "RV+"].reset_index()

    ddct_final = pd.merge(dct_plus, dct_minus, on="Time")
    ddct_final["ddCt"] = ddct_final["dCt_Test_x"] - ddct_final["dCt_Cal_y"]
    ddct_final["Relative Normalized Expression"] = 2**(-ddct_final["ddCt"])
    return ddct_final


def dir_path(path):
    path = path.split('/')[0:-1]
    out_dir = ''
    for i in path:
        out_dir += str(i + '/')
    return out_dir


def all_3_ADAR_genes(sample_file1, sample_file2, sample_file3):
    """
    :param sample_file1: path of ddCt analysis by just_ddCt function i.e py.ADAR1.dCt file
    :param sample_file2: path of ddCt analysis by just_ddCt function i.e py.ADAR2.dCt file
    :param sample_file3: path of ddCt analysis by just_ddCt function i.e py.ADAR3.dCt file
    :return:
    """

    # Load the Data

    label_sample1 = "ADAR1"
    print("loading " + sample_file1 + " as sample")
    data_adar1 = pd.read_csv(sample_file1)
    data_adar1["Target"] = label_sample1

    label_sample2 = "ADAR2"
    print("loading " + sample_file2 + " as sample")
    data_adar2 = pd.read_csv(sample_file2)
    data_adar2["Target"] = label_sample2

    label_sample3 = "ADAR3"
    print("loading " + sample_file3 + " as sample")
    data_adar3 = pd.read_csv(sample_file3)
    data_adar3["Target"] = label_sample3

    # Organize the Data

    data = pd.concat([data_adar1, data_adar2, data_adar3])
    filename = "all_ADAR_gene_FC.csv"
    path = dir_path(sample_file1)
    data.to_csv(path + "py." + filename, sep=",", encoding='utf-8')
    return data


def main():
    # Load the Data
    #ADAR1
    # sample_file1 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/20181224 ADAR1_Plate1/" \
    #                "Oded_2018-12-24 11-29-14_BR003827 -  Quantification Cq Results_0.csv"
    # sample_file2 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/20181225 ADAR1_Plate2/" \
    #                "Oded_2018-12-25 11-52-54_BR003827 -  Quantification Cq Results_0.csv"
    # sample_file3 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/20190414 ADAR1_Plate3/" \
    #                "Oded_2019-04-14 15-06-46_BR003827 -  Quantification Cq Results_0.csv"
    # gene = "ADAR1"
    #ADAR2
    # sample_file1 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/20190102 ADAR2_Plate1/" \
    #                "Oded_2019-01-02 13-10-55_BR003827 -  Quantification Cq Results_0.csv"
    # sample_file2 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/20190108 ADAR2_Plate2/" \
    #                "Oded_2019-01-08 11-20-07_BR003827 -  Quantification Cq Results_0.csv"
    # sample_file3 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/20190417 ADAR2_Plate3/" \
    #                "Oded_2019-04-17 11-35-52_BR003827 -  Quantification Cq Results_0.csv"
    # gene = "ADAR2"

    #ADAR3
    sample_file1 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/20181206 ADAR3_Plate1/" \
                   "Oded_2018-12-06 12-11-39_BR003827 -  Quantification Cq Results_0.csv"
    sample_file2 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/20181211 ADAR3_Plate2/" \
                   "Oded_2018-12-11 14-27-05_BR003827 -  Quantification Cq Results_0.csv"
    sample_file3 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/20190429 ADAR3_Plate3/" \
                   "Oded_2019-04-29 11-51-01_BR003827 -  Quantification Cq Results_0.csv"
    gene = "ADAR3"


    data = load_data_to_df(sample_file1, sample_file2, sample_file3, gene)
    data.to_csv("/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/Results/"+gene+"/all_3samples_raw.csv", sep=",",
                encoding='utf-8')

    input_dir = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/Results/"+gene+"/"
    df_list = replica_analyse(data, gene)

    for df in range(len(df_list)):
        # print(df_list[df])
        df_list[df].to_csv(input_dir+"Replica_"+str(df+1)+".csv", sep=",", encoding='utf-8')

    # dCT_file = input_dir+"Replica_1.csv"
        ddct_final = just_ddCT_replica(df_list[df])
        ddct_final.to_csv(input_dir+"Replica_"+str(df+1)+"ddct.csv", sep=",", encoding='utf-8')


    # dct_df = dct_cal_test_calculator(sample_file1, sample_file2, sample_file3, gene)
    # path = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/Results/ADAR2/"
    # print(path + gene + ".dCt.csv")
    # dct_df.to_csv(path + gene + ".dCt.csv", sep=",", encoding='utf-8')
    # dCT_file = path + gene + ".dCt.csv"
    # ddct_final = just_ddCT(dCT_file)
    # filename = dCT_file.split("/")[-1]
    # path = dir_path(dCT_file)
    # ddct_final.to_csv(path + "py." + filename, sep=",", encoding='utf-8')


    # sample_file1 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/py.ADAR1.dCt.csv"
    # sample_file2 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/py.ADAR2.dCt.csv"
    # sample_file3 = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/py.ADAR3.dCt.csv"

    #all_3_ADAR_genes(sample_file1, sample_file2, sample_file3)

    #Plot in R script - /Users/odedkushnir/Google Drive/Studies/PhD/R_Scripts/RealTime_analysis_ALL_ADAR.R


if __name__ == "__main__":
    main()