from Bio.Seq import Seq

def barcode_to_alignment_score(barcode: str, match_score: float = 2, flanking_length: int = 5):

    #search for at least 8 N's in the provided fasta
    primer_id_template=r'(N{8,})'
       
    m = re.search(primer_id_template, str(record.seq))
    #reduce anchors to max of 20 bases from each side of the barcode
    recZone = record.seq[max(0,m.start(0)-20):min(len(record.seq)-1,m.end(0)+20)]
    m = re.search(primer_id_template, str(recZone))
    m.start(0), m.end(0)

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())