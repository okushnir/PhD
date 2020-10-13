from skbio.util import classproperty
from skbio.sequence import GrammaredSequence
import Bio.Data.IUPACData

class DnaGrammaredSequence(GrammaredSequence):
    @classproperty
    def definite_chars(cls):
        return set(
            Bio.Data.IUPACData.unambiguous_dna_letters
        )

    @classproperty
    def default_gap_char(cls):
        return "-"

    @classproperty
    def gap_chars(cls):
        return set("-.")

class DegenerateDNA(DnaGrammaredSequence):
    @classproperty
    def degenerate_map(cls):
        return {
            key: set(value)
            for (
                key,
                value,
            ) in Bio.Data.IUPACData.ambiguous_dna_values.items()
            if key
            not in Bio.Data.IUPACData.unambiguous_dna_letters
        }


