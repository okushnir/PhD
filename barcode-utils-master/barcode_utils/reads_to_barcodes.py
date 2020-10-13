import re
from typing import Iterable

from Bio.Data.IUPACData import ambiguous_dna_letters

from .align_degenerate_dna import BarcodeAligner


class ReadsBarcodeExtractor(object):
    def __init__(self, barcode: str, flank_size: int = 5, match_len: int = None, **kwargs):
        self.barcode = barcode
        self.flank_size = flank_size
        self.match_len = match_len
        self._min_score = self.calc_min_score()
        self._aligner = BarcodeAligner(self.barcode)

    @property
    def flank_size(self):
        return self._flank_size

    @flank_size.setter
    def flank_size(self, flank_size: int):
        if not isinstance(flank_size, int):
            raise ValueError("{} must be a non-negative integer. Value provided: {}".format("flank_size", flank_size))
        if flank_size <= 0:
            raise ValueError("{} is not a valid flank size. Use non-negative values".format(flank_size))
        min_flank_length = min([len(s) for s in self.barcode.split("N"*self._len_degenerate)])
        if flank_size > min_flank_length:
            raise ValueError("{} is not a valid flank length as it can't be satisfied with provided barcode used {}. Use smaller values".format(flank_size, self.barcode))
        self._flank_size = flank_size

    @property
    def match_len(self):
        return self._match_len

    @match_len.setter
    def match_len(self, match_len: int):
        if match_len is None:
            match_len = (self.flank_size * 2 + self._len_degenerate)
        elif not isinstance(match_len, int):
            raise ValueError("{} must be a non-negative integer. Value provided: {}".format("match_len", match_len))
        elif match_len <= 0:
            raise ValueError("{} is not a valid match length. Use non-negative values".format(match_len))
        elif match_len < (self.flank_size * 2 + self._len_degenerate):
            raise ValueError(
                "{} is not a valid match length, as it is shorter than {}, which is twice the flank size {} and degenerate region length {}. Alter input values".format(
                    match_len, (self.flank_size * 2 + self._len_degenerate), self.flank_size, self._len_degenerate
                )
            )        
        if match_len > len(self.barcode):
            raise ValueError("{} is not a valid match length, as it is longer than barcode region of {}".format(match_len, len(self.barcode)))
        self._match_len = match_len

    @property
    def barcode(self):
        return self._barcode

    @barcode.setter
    def barcode(self, barcode: str):
        if not all([x in ambiguous_dna_letters for x in set(barcode)]):
            raise ValueError("{} is not a valid barcode".format(barcode))
        if barcode.count("NN") == 0:
            raise ValueError("{} does not contain a sequence of at least 2 N's".format(barcode))
        self._barcode = barcode

        n_template = r'(N{2,})'
        self._len_degenerate = len(max(re.findall(n_template, barcode)))

    def calc_min_score(self):
        return self._len_degenerate*2 + self.flank_size*2*4

    def extract_barcodes(self, reads: Iterable[str], min_score: int = None):
        if not min_score:
            min_score = self._min_score
        ext_reads = filter(lambda i: i[1] > min_score, self._aligner.align(reads))
        def extract_barcode(msa):
            gap_char=msa[0].default_gap_char
            leftmost=str(msa[0]).find("N")
            rightmost=str(msa[0]).rfind("N")
            if leftmost == -1 or leftmost==rightmost:
                return None
            #expand next to barcode to support insertions
            for i in reversed(range(leftmost)):
                if str(msa[0])[i] == gap_char:
                    leftmost -= 1 
                else:
                    break
            if leftmost < self.flank_size or rightmost + self.flank_size > len(msa[0]): #Enforce barcode flanking
                return None
            if "N" in str(msa[1])[leftmost:rightmost+1]:
                return None
            return str(msa[1])[leftmost:rightmost+1].replace(gap_char,"")

        ext_barcodes = filter(lambda x: x is not None, (extract_barcode(i[0]) for i in ext_reads))
        return ext_barcodes