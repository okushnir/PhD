#! /powerapps/share/python-anaconda-3.2019.7/bin/python3.7

"""
@Author: odedkushnir

"""

import sys, argparse
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def main():
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/190627_RV_CV/merged/RVB14"
    data_p2 = pd.read_csv(input_dir + "/RVB14_p2/20191029_q38/co_occur.csv")
    data_p5 = pd.read_csv(input_dir + "/RVB14_p5/20191029_q38/co_occur.csv")
    data_p8 = pd.read_csv(input_dir + "/RVB14_p8/20191029_q38/co_occur.csv")
    data_p10 = pd.read_csv(input_dir + "/RVB14_p10/20191029_q38/co_occur.csv")
    # data_p12 = pd.read_csv(input_dir + "/RVB14_p12/20191029_q38/co_occur.csv")
    df = pd.concat([data_p2, data_p5, data_p8, data_p10])

    #Plot
    g1 = sns.relplot(x="Pos", y="Frequency", data=df, hue="label")#, style="Stretch")
    g1.set(yscale="log")
    plt.show()
    mnkdd

if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    # parser.add_argument("virus", type=str, help="name of the virus, RVB14")
    # parser.add_argument("passages", type=str, help="from which passages, p0-p12")
    # parser.add_argument("without", type=int, help="Exclude passage no.")
    # parser.add_argument("input_dir", type=str, help="the path to the directory that contains data_mutation.csv")
    # parser.add_argument("quality", type=str, help="what is the prefix for the data_mutation.csv file; quality of the pipline ; for example: q38")
    # parser.add_argument("mutation_type", type=str, help="all/syn")
    # args = parser.parse_args(sys.argv[1:])
    # main(args)
    main()