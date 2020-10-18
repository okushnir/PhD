import sys, argparse
import pandas as pd
import numpy as np
from scipy import stats
import glob


def merege_freqs_variant(input_path, sample, out_file):

    freqs_path = input_path + "/%s.freqs" % (sample)



    freqs_data = pd.read_table(freqs_path)
    variant_path = input_path + "/%svsControl.csv" % (sample)
    variant_data = pd.read_csv(variant_path)
    variant_data = variant_data.rename(columns={"Ref_Pos": "Pos"})
    variant_data = variant_data.rename(columns={"Var": "Base"})
    variant_data = variant_data.rename(columns={"Cons": "Ref"})

    data = pd.merge(freqs_data, variant_data, how="outer", on=("Pos", "Base", "Ref"))

    return data



def main():
    virus = "RVB14"
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/patients"
    dirs = glob.glob(input_dir + "/Patient*")
    med_dir = "/20201017_q30_consensusX5/"
    for passage in dirs:
        input_path = glob.glob(passage + med_dir)
        # "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/RVB14_p5_L001-ds.b723553755a44959a7dc6debdb4558ea/q38_3UTR"
        sample = passage.split("/")[-1].replace("_", "-")#split("p")[0] + passage.split("/")[-1].split("_p")[1].split("_")[0]
        out_file = passage + med_dir + sample + ".merged.freqs"
        # print(type(sample))
        data = merege_freqs_variant(input_path[0], sample, out_file)
        data.to_csv(out_file, sep='\t', encoding='utf-8')


    # sample = "CVB3-RNA-Control"
    # input_path = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/CVB3/CVB3_RNA_Control_L001-ds.738551d760354e7fb0f07567e3ac1f70/q38_3UTR/"
    # out_file = input_path + sample + ".merged.freqs"
    # data = merege_freqs_variant(input_path, sample, out_file)
    # data.to_csv(out_file, sep='\t', encoding='utf-8')



if __name__ == "__main__":
    main()