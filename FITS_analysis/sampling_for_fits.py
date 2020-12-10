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
    choice = np.random.choice(record_lst size=choice_size, replace=False)
    with open(out_file, "w") as out_handle:
        for record in SeqIO.parse(fastq_file, "fastq"):
            if record.id in choice:
                SeqIO.write(record, out_handle, "fastq")
            else:
                continue
    return choice


def main():
    # fastq_file = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/p10_1/test.fastq"
    # out_path = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/p10_1/out.fastq"
    # print(choice)
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/"
    fastq_lst = ["p2_1", "p2-2", "p2_3", "p5_1", "p5_2", "p5_3", "p8_1", "p8_2", "p8_3", "p10_2", "p10_3", "p12_1",
                 "p12_2", "p12_3"]
    for fastq in fastq_lst:
        sample = fastq.replace("_", "-")
        out_file = sample + "_fits.fastq"
        out_path = input_dir + out_file
        choice_size = 572,872
        choice = sampling(fastq_file, out_path, choice_size)

if __name__ == "__main__":
    main()
