

import os
import pbs_runners
import glob
import pandas as pd
import sys
from Bio import SeqIO

def main():
    #cluster
    csv_file = "/sternadi/home/volume3/okushnir/RNAseq/180725_M04473_0026_000000000-BT9GF/index.csv"
    fastq_path = "/sternadi/home/volume3/okushnir/RNAseq/180725_M04473_0026_000000000-BT9GF/Undetermined_S0_L001_R1_001.fastq"
    output_dir = "/sternadi/home/volume3/okushnir/RNAseq/180725_M04473_0026_000000000-BT9GF/indexed"

    #local
    # csv_file = "/Volumes/STERNADILABHOME$/volume3/okushnir/RNAseq/180725_M04473_0026_000000000-BT9GF/index.csv"
    # fastq_path = "/Volumes/STERNADILABHOME$/volume3/okushnir/RNAseq/180725_M04473_0026_000000000-BT9GF/Undetermined_S0_L001_R1_001.fastq"
    # output_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/RNAseq/180725_M04473_0026_000000000-BT9GF/indexed"

    i7_index_no = "N704"
    i5_index_no = "S504"

    fetch_index(i7_index_no, i5_index_no, csv_file, fastq_path, output_dir)


def fetch_index(i7_index_no, i5_index_no, csv_file, fastq_path, output_dir):
    index_table = pd.read_csv(csv_file)
    i7_seq = index_table[index_table["i7_Index_Name"] == i7_index_no]["i7_Bases_Adapter"].values.item()
    # print(i7_seq)
    i5_seq = index_table[index_table["i5_Index_Name"] == i5_index_no]["i5_Bases_Adapter"].values.item()
    # print(i5_seq)
    record_dict = SeqIO.to_dict(SeqIO.parse(fastq_path, "fastq"))
    for key, val in record_dict.items():
        flag7 = val.seq.find(i7_seq)
        flag5 = val.seq.find(i5_seq)
        # print(flag7)
        # print(flag5)
        if (flag7 >= 0) and (flag5 >= 0):
            with open("/Volumes/STERNADILABHOME$/volume3/okushnir/RNAseq/180725_M04473_0026_000000000-BT9GF/example.fastq", "a") as output_handle:
                SeqIO.write(val, output_handle, "fastq")



if __name__ == "__main__":
    main()