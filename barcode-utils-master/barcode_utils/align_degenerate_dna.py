import warnings
from typing import Iterable
from skbio.util import EfficiencyWarning

from skbio.sequence import DNA
from skbio.alignment import local_pairwise_align

from .degenerate_dna import DegenerateDNA
from .utils import reverse_complement

def sub_mat(
    grammared_sequence, match_score, mismatch_score
):

    alphabet = list(
        grammared_sequence.definite_chars
    ) + list(grammared_sequence.degenerate_map.keys())
    is_definite = {
        x: (x in grammared_sequence.definite_chars)
        for x in alphabet
    }

    result = {}
    for c1 in alphabet:
        row = {}
        c1_set = (
            set(c1)
            if is_definite[c1]
            else grammared_sequence.degenerate_map[c1]
        )

        for c2 in alphabet:
            c2_set = (
                set(c2)
                if is_definite[c2]
                else grammared_sequence.degenerate_map[c2]
            )
            row[c2] = (
                (match_score if is_definite[c1] else (match_score - 0.001))
                if c1_set.intersection(c2_set)
                else mismatch_score
            )
        result[c1] = row
    return result


def degenerate_dna_pairwise_aligner(seq1, seq2: str, substitution_matrix):
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=EfficiencyWarning)
        forward = DegenerateDNA(seq2)
        lpa1 = local_pairwise_align(
            seq1, forward, 3, 1, substitution_matrix
        )
        reverse = DegenerateDNA(reverse_complement(seq2))
        lpa2 = local_pairwise_align(
            seq1, reverse, 3, 1, substitution_matrix
        )        
        if lpa1[1] >= lpa2[1]:
            return lpa1[0], lpa1[1]  # TabularMSA, score
        return lpa2[0], lpa2[1]


class BarcodeAligner(object):
    def __init__(
        self,
        barcode,
        match_score=2,
        mismatch_score=-2,
        gap_penalty=2,
    ):
        self.barcode = DegenerateDNA(barcode)
        self.substitution_mat = sub_mat(
            DegenerateDNA,
            match_score=match_score,
            mismatch_score=mismatch_score,
        )

    def align(self, reads: Iterable[str]):
        for read in reads:
            yield degenerate_dna_pairwise_aligner(
                self.barcode, read, self.substitution_mat
            )
