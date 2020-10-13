import sys, contextlib

from Bio.Seq import Seq 
from Bio import SeqIO

from .reads_to_barcodes import ReadsBarcodeExtractor

@contextlib.contextmanager
def smart_open(filename=None):
    if filename and filename != '-':
        fh = open(filename, 'w')
    else:
        fh = sys.stdout
    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()

def main(args):
    input_fastq= args.fastq
    input_fasta= args.barcode
    output_file= args.output

    barcode = str(SeqIO.read(input_fasta, "fasta").seq)
    extractor = ReadsBarcodeExtractor(barcode)

    barcodes=[]
    with open(input_fastq, "r") as fastq_handle:
        with smart_open(output_file) as output_handle:
            for seqRecord in SeqIO.parse(fastq_handle, "fastq"):
                toAlign=seqRecord.seq
                barcode=list(extractor.extract_barcodes([str(toAlign)]))
                if barcode:
                    barcodes.append(barcode[0])
                    seqRecord.description += " PID:{}".format(barcode[0])
                print(seqRecord.format("fastq"), file=output_handle, end='')
    return len(barcodes)

    

