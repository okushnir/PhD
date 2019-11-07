import pandas as pd
from Bio import SeqIO

def subtract_file(path):
    with open(path, "r") as handle:
        file = handle.readlines()
        print(len(file))
        no_line = len(file)%4
        no_line = -(no_line)
        print(no_line)
        if no_line != 0:
            file = file[0:no_line]
            # print(file)
            with open(path, "w") as output_handle:
                output_handle.writelines(file)
        else:
            return path

def main():
    path = "/Users/odedkushnir/Projects/fitness/Pra_SRA/ERP014415_Entero_A/ERR2352267/fastq/ERR2352267_1 copy.fastq"
    subtract_file(path)


if __name__ == "__main__":
    main()