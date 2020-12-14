
from Bio import SeqIO


def fasta_convert(file_path, output_path):
    with open(file_path, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            input_record = record

    with open(output_path, "w") as output_handle:
        SeqIO.write(input_record, output_handle, "fasta")


def main():
    file_path = "/Volumes/STERNADILABHOME$/volume3/okushnir/ref/CVB3/CVB3_from_pT7CVB3_1-7399.fasta"
    output_path = "/Volumes/STERNADILABHOME$/volume3/okushnir/ref/CVB3/new_CVB3_from_pT7CVB3_1-7399.fasta"
    fasta_convert(file_path, output_path)

if __name__ == "__main__":
    main()
