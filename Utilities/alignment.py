from Bio import pairwise2
from Bio import SeqIO
from Bio.pairwise2 import format_alignment
from Bio import AlignIO

with open("/Volumes/STERNADILABHOME$/volume3/okushnir/SRP/rvb14_ncbi.fasta", "rU") as RVB14:
    for record in SeqIO.parse(RVB14, "fasta"):
        rvb14 = (record.seq)
with open("/Volumes/STERNADILABHOME$/volume3/okushnir/SRP/JF781502.1_SRR185824.fasta", "rU") as JF781502:
    for record in SeqIO.parse(JF781502, "fasta"):
        jf781502 = (record.seq)

# alignments = pairwise2.align.localxx(rvb14, jf781502)
# with open('/Volumes/STERNADILABHOME$/volume3/okushnir/SRP/alignment.txt', 'w') as f:
#     f.write(format_alignment(*alignments[0]))

with open('/Volumes/STERNADILABHOME$/volume3/okushnir/SRP/alignment.fasta', 'w') as f:
    AlignIO.write(AlignIO.parse(alignments, "stockholm"), f, format="stockholm")
