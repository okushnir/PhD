#! /usr/local/python-anaconda-2.7//bin/python


import textwrap
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import GC


def main():
    # More generic
    # file_path = (input("Enter the file path of the sequence you want to analyze:"))
    file_path = "Z:\\volume1\\okushnir\\HRVB14\\Seq\\HRVB14_from_pWR3.26.txt"        # For PyCharm
    # For Linux
    # file_path = "/sternadi/home/volume1/okushnir/Scripts/HRVB14 from pWR3.26.txt"

    my_file = open(file_path)
    bio_seq = my_file.read()
    the_seq = textwrap.fill(bio_seq)
    my_file.close()

    # fasta file
    # for seq_record in SeqIO.parse("Z:\\volume1\\okushnir\\Scripts\\HRVB14 from pWR3.26.fasta", "fasta"):
    #     print(seq_record.id)
    #     print(repr(seq_record.seq))
    #     print(len(seq_record))

    bio_seq = Seq(bio_seq, generic_dna)

    # print(the_seq)

    total_len = the_seq.count("A") + the_seq.count("T") + the_seq.count("C") + the_seq.count("G")
    no_ct = the_seq.count("CT")
    no_ca = the_seq.count("CA")
    no_cc = the_seq.count("CC")
    no_gt = the_seq.count("GT")
    no_ga = the_seq.count("GA")
    no_gg = the_seq.count("GG")
    no_CG = the_seq.count("CG")
    no_GC = the_seq.count("GC")

    # input("Enter the paired nucleotides you want to search:")
    search_nuc = 'CG'

    if search_nuc in the_seq:
        print(the_seq.replace(search_nuc, '\033[44;33m{}\033[m'.format(search_nuc)))
    else:
        print("The nucleotide is not in the sequence")

    # Counts the paired nuc that was inserted
    no_search_nuc = the_seq.count(search_nuc)
    print "The No. of your searched paired nucleotides: %i" % no_search_nuc
    # print(highlight('The No. of your searched paired nucleotides:', color=0, bold=1),
    #       highlight(no_search_nuc, color=0, bold=1),)
    print "_____________________________________________________________________"
    print "No of CT = %i" % no_ct
    print "No of CA = %i" % no_ca
    print "No of CC = %i" % no_cc
    print "_____________________________________________________________________"
    print "No of GT = %i" % no_gt
    print "No of GA = %i" % no_ga
    print "No of GG = %i" % no_gg
    print('_____________________________________________________________________')
    print('No of CG = %i' % no_CG)
    print('No of GC = %i' % no_GC)
    """
    for python 3.5
    print(highlight('No of CG =', color='green', bold=1), highlight(no_CG, color='green', bold=1))
    print(highlight('No of GC =', color=1, bold=1), highlight(no_GC, color=1, bold=1))
    """
    print('_____________________________________________________________________')
    print "The length of your sequence is:%i%s" % (total_len, 'bp')
    print('The GC Content is:%i%c' % (GC(the_seq), '%'))
    print('_____________________________________________________________________')

    # coding_dna = bio_seq.reverse_complement()
    # print(coding_dna)
    # print('_____________________________________________________________________')
    m_rna = bio_seq.transcribe()
    # print(m_rna)
    # print('_____________________________________________________________________')

    pro_seq = str(bio_seq.translate(table=1))
    print(textwrap.fill(pro_seq))


def highlight(string, color, bold):
        attr = []
        if color == 'green':
            # green
            attr.append('32')
        else:
            # red
            attr.append('31')
        if bold:
            attr.append('1')
        return '\x1b[%sm%s\x1b[0m' % (';'.join(attr), string)


if __name__ == "__main__":
    main()


