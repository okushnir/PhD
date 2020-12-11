#!/powerapps/share/python-anaconda-3.2019.7/bin/python

import sys, argparse
import os
import numpy as np
import pandas as pd
from Bio import SeqIO


def sampling(fastq_file, out_file, choice_size):
    record_lst= []
    for record in SeqIO.parse(fastq_file, "fastq"):
        record_lst.append(record.id)
    # record_df = pd.DataFrame(record_lst)
    # record_df = record_df.rename(columns={0: "record.id"})
    print(len(record_lst))
    choice = np.random.choice(record_lst, size=choice_size, replace=False)
    with open(out_file, "w") as out_handle:
        for record in SeqIO.parse(fastq_file, "fastq"):
            if record.id in choice:
                SeqIO.write(record, out_handle, "fastq")
            else:
                continue
    return choice


def main(args):
    # fastq_file = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/p10_1/test.fastq"
    # out_path = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/p10_1/out.fastq"
    # print(choice)
    # input_dir = "/sternadi/home/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/"
    # fastq_lst = ["p2_1", "p2_2", "p2_3", "p5_1", "p5_2", "p5_3", "p8_1", "p8_2", "p8_3", "p10_2", "p10_3", "p12_1",
    #              "p12_2", "p12_3"]
    # for fastq in fastq_lst:
    #     sample = fastq.replace("_", "-")
    #     # "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/p2_1/p2-1_merged.fastq"
    #     fastq_path = input_dir + fastq + "/" + sample + "_merged.fastq"
    #     out_file = sample + "_fits.fastq"
    #     out_path = input_dir + fastq + "/" + out_file
    # choice_size = 572, 872
    fastq_file = args.fastq_file
    out_path = args.out_path
    choice_size = args.choice_size
    sampling(fastq_file, out_path, choice_size=choice_size)
    print("Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_file", type=str, help="fastq_file path")
    parser.add_argument("out_path", type=str, help="ouput file path")
    parser.add_argument("-c", "--choice_size", type=int, help="number of choices", required=False, default=572872)
    args = parser.parse_args(sys.argv[1:])
    main(args)
    main()
