import sys, argparse

from .io import main

parser = argparse.ArgumentParser()
parser.add_argument("fastq", type=str, help="input reads file (fastq format)")
parser.add_argument("barcode",type=str, help="input barcode file (fasta format)")
parser.add_argument("-o","--output", type=str, help="identified barcodes file", default=None, required=False)
 
args = parser.parse_args(sys.argv[1:])
num_barcodes = main(args)
print ("found {} barcodes".format(num_barcodes))
