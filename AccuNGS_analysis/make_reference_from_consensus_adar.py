#! /usr/local/python_anaconda/bin/python3.4
import sys, gzip
from Bio import SeqIO,Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser
import pandas as pd
import numpy as np

def replace_base(row, mut_reference):
    mut_reference[int(row["Pos"])-1]=row["Base"]

def main():
    # parser = OptionParser("usage: %prog [options]\nTry running %prog --help for more information")
    #     # parser.add_option("-f", "--fasta", dest="fasta", help="reference fasta file")
    #     # parser.add_option("-p", "--freqs", dest="freqs", help="frequency file")
    #     # parser.add_option("-o", "--output_file", dest="output_file", help="output fasta file")
    #     # parser.add_option("-c", "--min_coverage", dest="coverage", help="minimal coverage of position to be updated (default: 1000)", type="int")
    #     # (options, args) = parser.parse_args()
    #     # fasta = options.fasta
    #     # freqs = options.freqs
    #     # output_file = options.output_file
    fasta = "/Volumes/STERNADILABHOME$/volume3/okushnir/ref/RVB14/HRVB14_from_pWR3.26_1-7212.fasta"
    freqs = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/p12_2/20201012_q38/p12-2.freqs"
    data_mutation = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/q38_data_mutation.csv"
    output_file = "/Volumes/STERNADILABHOME$/volume3/okushnir/ref/RVB14/p12_consensus_hot_spot.fasta"

    mutable_reference = None
    min_coverage = 1000
    # if options.coverage:
    #     min_coverage = options.coverage

    # if options.fasta is None or options.fasta is None or options.output_file is None:
    #     parser.error("Missing file input")

    reference = None
    input_record = None

    with open(fasta, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            reference=record.seq
            input_record=record

    mutable_reference=MutableSeq(str(reference))
    data_freqs2 = pd.read_csv(data_mutation)
    data_freqs2["passage"] = data_freqs2["passage"].astype(int)
    data_freqs2 = data_freqs2[data_freqs2["passage"] == 12]
    data_freqs2["replica"] = data_freqs2["replica"].astype(int)
    data_freqs2 = data_freqs2[data_freqs2["replica"] == 2]
    data_freqs2 = data_freqs2[data_freqs2["Mutation"] == "A>G"]

    data_freqs = pd.read_table(freqs)
    data_freqs = data_freqs[data_freqs["Read_count"] >= min_coverage]

    # data_freqs["Base"] = data_freqs["Base"].astype(str)
    # data_freqs["Ref"] = data_freqs["Ref"].astype(str)
    # data_freqs["Rank"] = data_freqs["Rank"].astype(str)

    data_freqs = data_freqs[data_freqs["Rank"] == 0]
    #data_freqs = data_freqs[data_freqs["Base"] != data_freqs["Ref"]]
    data_freqs = data_freqs[data_freqs["Pos"] == np.round(data_freqs['Pos'])] #remove insertion
    data_freqs = data_freqs[data_freqs["Base"] != "-"] #remove deletion

    data_freqs = data_freqs.merge(data_freqs2, how="outer", on="Pos")
    """half of the genome ADAR-like"""
    # data_freqs["Base_x"] = np.where(data_freqs["Base_y"] == "G", "G", data_freqs["Base_x"])
    # data_freqs = data_freqs[["Pos", "Base_x", "Freq_x", "Ref_x", "Read_count_x", "Rank_x", "Prob_x"]]
    data_freqs = data_freqs.rename(columns={"Base_x": "Base", "Freq_x": "Freq", "Ref_x": "Ref", "Read_count_x": "Read_count", "Rank_x": "Rank", "Prob_x": "Prob"})
    """One stretch ADAR-like"""
    data_freqs["Pos"] = data_freqs["Pos"].astype(int)
    data_freqs["Base"] = np.where((data_freqs["Pos"] == 6019), "G", data_freqs["Base"])
    data_freqs["Base"] = np.where((data_freqs["Pos"] == 6043), "G", data_freqs["Base"])
    data_freqs["Base"] = np.where((data_freqs["Pos"] == 6049), "G", data_freqs["Base"])
    data_freqs["Base"] = np.where((data_freqs["Pos"] == 6067), "G", data_freqs["Base"])
    data_freqs["Base"] = np.where((data_freqs["Pos"] == 6076), "G", data_freqs["Base"])
    data_freqs["Base"] = np.where((data_freqs["Pos"] == 6078), "G", data_freqs["Base"])
    data_freqs["Base"] = np.where((data_freqs["Pos"] == 6079), "G", data_freqs["Base"])
    data_freqs["Base"] = np.where((data_freqs["Pos"] == 6082), "G", data_freqs["Base"])
    data_freqs["Base"] = np.where((data_freqs["Pos"] == 6097), "G", data_freqs["Base"])
    data_freqs.apply(replace_base, args=(mutable_reference,), axis=1) #updates mutable reference to hold correct consensus


    new_sequence=Seq.Seq(str(mutable_reference))
    input_record.seq=new_sequence
    with open(output_file, "w") as output_handle:
        SeqIO.write(input_record, output_handle, "fasta")

if __name__ == "__main__":
    main()