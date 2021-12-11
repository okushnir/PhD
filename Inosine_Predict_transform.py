#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import pandas as pd


def convert_predict_to_csv_file(input_file, output_flle):
    table = pd.read_table(input_file, header=0)
    table = table.transpose()
    table = table.reset_index()
    table.columns = table.iloc[0]
    table = table.iloc[1:]

    table.to_csv(output_flle, sep=',', encoding='utf-8', index=False)


def main():
    virus = "CVB3"
    input_file = "/Users/odedkushnir/PhD_Projects/After_review/AccuNGS/CVB3/InosinePredict_reults/%s_adar1.txt" % virus
    output_file = "/Users/odedkushnir/PhD_Projects/After_review/AccuNGS/CVB3/InosinePredict_reults/Output/%s_adar1_trans.csv" % virus
    convert_predict_to_csv_file(input_file, output_file)

if __name__ == "__main__":
    main()