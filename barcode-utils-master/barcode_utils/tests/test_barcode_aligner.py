import re, pytest

from Bio.Seq import Seq

from ..align_degenerate_dna import BarcodeAligner
from ..reads_to_barcodes import ReadsBarcodeExtractor


class TestBarcodeAligner(object):

    def test_barcode_score(self):

        seq_with_barcode="GTTCACTGCCATCCTGGTCGAGCTACCCATACGATGTTCCAGATTACGTAAATGATGTGACGGGCAAATGCAAACATGTCCAAAAG"
        barcode_in_seq="TAAATGATGTGACGG"
        barcode_degenerate_sequence="TACCCATACGATGTTCCAGATTACGNNNNNNNNNNNNNNNGCAAATGCAAACATGTCCAAAA"

        primer_id_template=r'(N{8,})'
        m = re.search(primer_id_template, barcode_degenerate_sequence)
        #reduce anchors to max of 20 bases from each side of the barcode
        recZone = barcode_degenerate_sequence[max(0,m.start(0)-20):min(len(barcode_degenerate_sequence)-1,m.end(0)+20)]

        aligner = BarcodeAligner(recZone)
        (aln, score) = list(aligner.align([seq_with_barcode]))[0]
        assert round(score) == 110

        aligner = BarcodeAligner(recZone)
        (aln, score) = list(aligner.align(["GTATAAGCACGATGATATGAATAATAGATGATTGAGAGTGTCCCAGGTCCTCTCTAGTGTCTTGGCGGTGCGTTGGTCCTTGGTTTTGGACATGGTTGCATTTGCATCCCCTGCTGAGCTCGTAATCTGGAAC"]))[0]
        assert score > 55

class TestReadsBarcodeExtractor(object):
    
    invalid_barcode="ABCDEFG"
    short_invalid_barcode="TTTTTTTT"
    barcode="TAAATGATGTGACGG"
    pre_barcode="TACCCATACGATGTTCCAGATTACG"
    post_barcode="GCAAATGCAAACATGTCCAAAA"
    valid_barcode=pre_barcode+("N"*15)+post_barcode


    def test_flank_size(self):
        extractor = ReadsBarcodeExtractor(self.valid_barcode)
        assert extractor.flank_size == 5 #default value for flank_size

        extractor = ReadsBarcodeExtractor(self.valid_barcode, flank_size = 3)
        assert extractor.flank_size == 3 

        with pytest.raises(ValueError):
            extractor = ReadsBarcodeExtractor(self.valid_barcode, flank_size = 0)
                
        with pytest.raises(ValueError):
            extractor = ReadsBarcodeExtractor(self.valid_barcode, flank_size = -6)
                
        with pytest.raises(ValueError):
            extractor = ReadsBarcodeExtractor(self.valid_barcode, flank_size = 8.2)
                
        with pytest.raises(ValueError):
            extractor = ReadsBarcodeExtractor(self.valid_barcode, flank_size = "str")

        # ensure reads without enough flank_size are not identified
        read="GGATATTTTAGAGCAAATGCAACCATGTCCAAAAAATCCACCAAAAAAAGATGATTACCATTTTGAAGTGTTCAACTTTGTTCCCTATAGT"
        barcode="ATACGATGTTCCAGATTACGNNNNNNNNNNNNNNNGCAAATGCAAMCATGTCCAAAA"
        extractor = ReadsBarcodeExtractor(barcode, flank_size = 1)
        assert len(list(extractor.extract_barcodes([read], 30))) == 0

        read="ATACGATGTTCCAGATTACGGGCGGATATTTTAGAGCAAATGCAACCATGTCCAAAAAATCCACCAAAAAAAGATGATTACCATTTTGAAGTGTTCAACTTTGTTCCCTATAGT"
        extractor = ReadsBarcodeExtractor(barcode, flank_size = 5)
        assert len(list(extractor.extract_barcodes([read], 30))) == 1

        with pytest.raises(ValueError):
            extractor = ReadsBarcodeExtractor(barcode, flank_size = 25)


    def test_match_length(self):
        extractor = ReadsBarcodeExtractor(self.valid_barcode)
        assert extractor.match_len == 25 #calculated value for match_len with this barcode and default flank size

        extractor = ReadsBarcodeExtractor(self.valid_barcode, flank_size = 10)
        assert extractor.match_len == 35 #calculated value for match_len with this barcode and flank size

        for match_len in [0, -6, 8.2, "str"]: #all sorts of illegal sizes
            with pytest.raises(ValueError):
                extractor = ReadsBarcodeExtractor(self.valid_barcode, match_len = match_len)

        with pytest.raises(ValueError):
            extractor = ReadsBarcodeExtractor(self.valid_barcode, match_len = 540)


    def test_barcode_validity(self):

        with pytest.raises(ValueError):
            ReadsBarcodeExtractor(self.invalid_barcode)

        with pytest.raises(ValueError):
            ReadsBarcodeExtractor(self.short_invalid_barcode)

        extractor = ReadsBarcodeExtractor(self.valid_barcode)
        assert extractor.barcode == self.valid_barcode

        with pytest.raises(ValueError):
            extractor.barcode=self.invalid_barcode


    def test_barcode_extraction_filter(self):
        extractor = ReadsBarcodeExtractor(self.valid_barcode)
        reads = ["GTTCACTGCCATCCTGGTCGAGCTACCCATACGATGTTCCAGATTACGTAAATGATGTGACGGGCAAATGCAAACATGTCCAAAAG"]
        # expected score on this barcode, upon perfect match is 124
        assert len(list(extractor.extract_barcodes(reads, 100))) == 1
        assert len(list(extractor.extract_barcodes(reads, 130))) == 0

        reads.append("ACGT"*4)
        assert len(list(extractor.extract_barcodes(reads, 100))) == 1
        assert len(list(extractor.extract_barcodes(reads, 130))) == 0


    def test_extracted_barcodes(self):
        extractor = ReadsBarcodeExtractor(self.valid_barcode)
        # one barcoded sequence, one without
        reads = [self.pre_barcode + self.barcode + self.post_barcode, "ACGTACGTACGTACGT"]
        assert list(extractor.extract_barcodes(reads, 100)) == [self.barcode]

	# 15-bases barcode to 15N template, 14, 16 and 17 bases as well
        for barcode in [self.barcode, "TAAAGATGTGACGG", "TAAATGCATGTGACGG", "TAAATGCATGTGATCGG"]:
            reads = [self.pre_barcode + barcode + self.post_barcode]
            assert list(extractor.extract_barcodes(reads, 100)) == [barcode]
        
        # Test reverse complement
        read = str(Seq(self.pre_barcode + self.barcode + self.post_barcode).reverse_complement())
        assert list(extractor.extract_barcodes([read], 100)) == [self.barcode]

        read = str(Seq("GTATAAGCACGATGATATGAATAATAGATGATTGAGAGTGTCCCAGGTCCTCTCTAGTGTCTTGGCGGTGCGTTGGTCCTTGGTTTTGGACATGGTTGCATTTGCATCCCCTGCTGAGCTCGTAATCTGGAAC"))
        assert list(extractor.extract_barcodes([read], 75)) == ["AGCTCAGCAGGGGAT"]

        read = str(Seq("GTATAAGCACGATGATATGAATAATAGATGATTGAGAGTGTCCCAGGTCCTCTCTAGTGTCTTGGCGGTGCGTTGGTCCTTGGTTTTGGACATGGTTGCATTTGCATCCCCTGCTGAGCTCGTA"))
        assert list(extractor.extract_barcodes([read], 75)) == []

        read = str(Seq("CTAGTGTCTTGGCGGTGCGTTGGTCCTTGGTTTTGGACATGGTTGCATTTGCGCAGATTGTCCTCATCGTAATCT"))
        assert list(extractor.extract_barcodes([read], 75)) == ["ATGAGGACAATCTGC"]

        read = str(Seq("CTAGTGTCTTGGCGGTGCGTTGGTCCTTGGTTTTGGACATGGTTGCATTTGCGCAGATTGTCCTCATCGTAATCT"))
        assert list(extractor.extract_barcodes([read])) == ["ATGAGGACAATCTGC"]


