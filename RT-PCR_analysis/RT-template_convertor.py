#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import pandas as pd


def convert_plate_to_csv_file(input_file, output_flle):
    table = pd.read_csv(input_file, header=0)
    temp_table = table.transpose()

    temp_table_A = temp_table.iloc[:, 0:3]
    temp_table_A.drop("Row", inplace=True)
    temp_table_A["Row"] = "A"
    temp_table_A = temp_table_A[0:12]
    temp_table_A["Column"] = temp_table_A.index
    temp_table_A = temp_table_A.rename(columns={0: "*Sample Name", 1: "*Target Name", 2: "*Biological Set Name"})
    temp_table_A = temp_table_A[["Row", "Column", "*Target Name", "*Sample Name", "*Biological Set Name"]]

    temp_table_B = temp_table.iloc[:, 3:6]
    temp_table_B.drop("Row", inplace=True)
    temp_table_B["Row"] = "B"
    temp_table_B = temp_table_B[0:12]
    temp_table_B["Column"] = temp_table_B.index
    temp_table_B = temp_table_B.rename(columns={3: "*Sample Name", 4: "*Target Name", 5: "*Biological Set Name"})
    temp_table_B = temp_table_B[["Row", "Column", "*Target Name", "*Sample Name", "*Biological Set Name"]]

    temp_table_C = temp_table.iloc[:, 6:9]
    temp_table_C.drop("Row", inplace=True)
    temp_table_C["Row"] = "C"
    temp_table_C = temp_table_C[0:12]
    temp_table_C["Column"] = temp_table_C.index
    temp_table_C = temp_table_C.rename(columns={6: "*Sample Name", 7: "*Target Name", 8: "*Biological Set Name"})
    temp_table_C = temp_table_C[["Row", "Column", "*Target Name", "*Sample Name", "*Biological Set Name"]]

    temp_table_D = temp_table.iloc[:, 9:12]
    temp_table_D.drop("Row", inplace=True)
    temp_table_D["Row"] = "D"
    temp_table_D = temp_table_D[0:12]
    temp_table_D["Column"] = temp_table_D.index
    temp_table_D = temp_table_D.rename(columns={9: "*Sample Name", 10: "*Target Name", 11: "*Biological Set Name"})
    temp_table_D = temp_table_D[["Row", "Column", "*Target Name", "*Sample Name", "*Biological Set Name"]]

    temp_table_E = temp_table.iloc[:, 12:15]
    temp_table_E.drop("Row", inplace=True)
    temp_table_E["Row"] = "E"
    temp_table_E = temp_table_E[0:12]
    temp_table_E["Column"] = temp_table_E.index
    temp_table_E = temp_table_E.rename(columns={12: "*Sample Name", 13: "*Target Name", 14: "*Biological Set Name"})
    temp_table_E = temp_table_E[["Row", "Column", "*Target Name", "*Sample Name", "*Biological Set Name"]]

    temp_table_F = temp_table.iloc[:, 15:18]
    temp_table_F.drop("Row", inplace=True)
    temp_table_F["Row"] = "F"
    temp_table_F = temp_table_F[0:12]
    temp_table_F["Column"] = temp_table_F.index
    temp_table_F = temp_table_F.rename(columns={15: "*Sample Name", 16: "*Target Name", 17: "*Biological Set Name"})
    temp_table_F = temp_table_F[["Row", "Column", "*Target Name", "*Sample Name", "*Biological Set Name"]]

    temp_table_G = temp_table.iloc[:, 18:21]
    temp_table_G.drop("Row", inplace=True)
    temp_table_G["Row"] = "G"
    temp_table_G = temp_table_G[0:12]
    temp_table_G["Column"] = temp_table_G.index
    temp_table_G = temp_table_G.rename(columns={18: "*Sample Name", 19: "*Target Name", 20: "*Biological Set Name"})
    temp_table_G = temp_table_G[["Row", "Column", "*Target Name", "*Sample Name", "*Biological Set Name"]]

    temp_table_H = temp_table.iloc[:, 21:24]
    temp_table_H.drop("Row", inplace=True)
    temp_table_H["Row"] = "H"
    temp_table_H = temp_table_H[0:12]
    temp_table_H["Column"] = temp_table_H.index
    temp_table_H = temp_table_H.rename(columns={21: "*Sample Name", 22: "*Target Name", 23: "*Biological Set Name"})
    temp_table_H = temp_table_H[["Row", "Column", "*Target Name", "*Sample Name", "*Biological Set Name"]]

    out_table = pd.concat([temp_table_A, temp_table_B, temp_table_C, temp_table_D, temp_table_E, temp_table_F, temp_table_G, temp_table_H], axis=0)
    out_table.to_csv(output_flle, sep=',', encoding='utf-8', index=False)


def main():
    # experiment = "new_time_line_qRT"
    # input_file = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/CV/CVB3/%s.csv" % experiment
    # output_file = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/CV/CVB3/Plate_template_%s.csv" % experiment
    # convert_plate_to_csv_file(input_file, output_file)

    experiment = "Capsid_Quant_qRT"
    direcotry = "RV"
    virus = "RVB14"
    date = "20200722"
    input_file = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/%s/%s/%s %s.csv" % (direcotry, virus, date, experiment)
    output_file = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/%s/%s/%s Plate_template_%s.csv" % (direcotry, virus, date, experiment)
    convert_plate_to_csv_file(input_file, output_file)

if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    # parser.add_argument("virus", type=str, help="name of the virus, RVB14")
    # parser.add_argument("passages", type=str, help="from which passages, p0-p12")
    # parser.add_argument("without", type=int, help="Exclude passage no.")
    # parser.add_argument("input_dir", type=str, help="the path to the directory that contains data_mutation.csv")
    # parser.add_argument("quality", type=str, help="what is the prefix for the data_mutation.csv file; quality of the pipline ; for example: q38")
    # parser.add_argument("mutation_type", type=str, help="all/syn")
    # args = parser.parse_args(sys.argv[1:])
    main()