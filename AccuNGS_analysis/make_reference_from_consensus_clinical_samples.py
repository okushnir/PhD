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
    min_coverage = 1000
    cycle = 5
    patient_lst = ("Patient_1", "Patient_4", "Patient_5", "Patient_9", "Patient_16", "Patient_17", "Patient_20")
    base_dir="/Volumes/STERNADILABHOME$/volume3/okushnir/"


    for patient in patient_lst:
        fasta = base_dir + ("ref/RVA/%s_consenX%s.fasta") % (patient, str(cycle))
        freqs = base_dir + ("AccuNGS/20201008RV-202329127/merged/patients/%s/20201017_q30_consensusX%s/%s.freqs") % \
                (patient, cycle, patient.replace("_", "-"))
        output_file = base_dir + ("ref/RVA/%s_consenX%s.fasta") % (patient, str(cycle+1))

        reference = None
        input_record = None
        mutable_reference = None

        with open(fasta, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                reference=record.seq
                input_record=record

        mutable_reference = MutableSeq(str(reference))

        data_freqs = pd.read_table(freqs)
        data_freqs = data_freqs[data_freqs["Read_count"] >= min_coverage]
        data_freqs = data_freqs[data_freqs["Rank"] == 0]
        #data_freqs = data_freqs[data_freqs["Base"] != data_freqs["Ref"]]
        data_freqs = data_freqs[data_freqs["Pos"] == np.round(data_freqs['Pos'])] #remove insertion
        data_freqs = data_freqs[data_freqs["Base"] != "-"] #remove deletion
        data_freqs.apply(replace_base, args=(mutable_reference,), axis=1) #updates mutable reference to hold correct consensus

        new_sequence=Seq.Seq(str(mutable_reference))
        input_record.seq=new_sequence
        with open(output_file, "w") as output_handle:
            SeqIO.write(input_record, output_handle, "fasta")

if __name__ == "__main__":
    main()