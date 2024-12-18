import pytest
from bio_functions import reverse_complement, translate

def test_reverse_complement_standard_sequences():
    # Test basic reverse complement functionality
    assert reverse_complement("ATGC") == "GCAT", "Failed to reverse complement 'ATGC'"
    assert reverse_complement("GCTAGC") == "GCTAGC", "Failed to reverse complement 'GCTAGC'"

def test_reverse_complement_with_invalid_characters():
    # Test reverse complement with invalid characters in sequence
    with pytest.raises(ValueError, match="DNA sequence contains invalid characters."):
        reverse_complement("ATXGCC")  # 'X' is invalid
    with pytest.raises(ValueError, match="DNA sequence contains invalid characters."):
        reverse_complement("ATGGNC")  # 'N' is invalid

def test_reverse_complement_with_empty_sequence():
    # Test reverse complement with an empty string
    assert reverse_complement("") == "", "Failed to handle empty sequence"

def test_translate_standard_sequences():
    # Test basic translation with complete codons
    assert translate("ATGGCTTCCTCCGAAGACGTTATCAAAGAGTTCATG") == "MASSEDVIKEFM", "Failed to translate mRFP no stop"
    assert translate("ATGGCTTCCTCCGAAGACGTTATCAAAGAGTTCATGTAA") == "MASSEDVIKEFM_", "Failed to translate mRFP with stop"
    
def test_translate_with_invalid_characters():
    # Test translate with invalid characters in sequence
    with pytest.raises(ValueError, match="DNA sequence contains invalid characters."):
        translate("ATGXCC")  # 'X' is invalid
    with pytest.raises(ValueError, match="DNA sequence contains invalid characters."):
        translate("ATGNCC")  # 'N' is invalid

def test_translate_with_non_multiple_of_three():
    # Test translate when sequence length is not a multiple of three
    with pytest.raises(ValueError, match="Length of DNA sequence is not a multiple of three"):
        translate("ATGGC")  # Incomplete codon

def test_translate_with_stop_codon():
    # Test for sequences containing stop codons
    assert translate("ATGTAA") == "M_", "Failed to translate stop codon 'TAA' in 'ATGTAA'"
    assert translate("TGGATATAG") == "WI_", "Failed to translate stop codon 'TAG' in 'TGGATTAG'"

def test_translate_with_empty_sequence():
    # Test with an empty sequence
    assert translate("") == "", "Failed to handle empty sequence"

def test_translate_with_lowercase():
    # Test translation with lowercase input
    assert translate("atggcc".upper()) == "MA", "Failed to handle lowercase input 'atggcc'"
