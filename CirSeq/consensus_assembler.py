import pandas as pd
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


def main():
    file_source = open("Z:/volume1/okushnir/Cirseq/pipeline_output_all_0_ref_plasmid/RVB14p2.freqs", "rU")
    get_consensus(file_source)


def get_consensus(file_source):
    data = pd.read_table(file_source)
    data["strain"] = "HRVB14"

    consensus_data = data[['Pos','Base']][data['Rank'] == 0]
    consensus_data = consensus_data[consensus_data['Pos'] == np.round(consensus_data['Pos'])]
    #sanity for all positions in consensus
    #for i in range(1,6541):
    #    if not (i in consensus_data["Pos"].values):
    #        print "False, ",i
    data2 = consensus_data.set_index("Pos")
    consensus = ""
    for i in range(1, 6541):
        consensus += str(data2.ix[i]["Base"])[0]

    out_seq = Seq(consensus, alphabet=IUPACUnambiguousDNA)
    out_io_seq = SeqRecord(out_seq, "HRVB14p2_from_cirseq", "HRVB14p2_from_cirseq")
    SeqIO.write(out_io_seq, "Z:/volume1/okushnir/Cirseq/pipeline_output_all_0_ref_plasmid/cirseq_consensus_2.fasta", "fasta")

    return out_seq


if __name__ == "__main__":
    main()