import pytest
from bio_functions import reverse_complement

def test_reverse_complement_standard_sequences():
    # Test basic reverse complement functionality
    assert reverse_complement("ATGC") == "GCAT", "Failed to reverse complement 'ATGC'"
    assert reverse_complement("GCTAGC") == "GCTAGC", "Failed to reverse complement 'GCTAGC'"
    assert reverse_complement("AATTCCGG") == "CCGGAATT", "Failed to reverse complement 'AATTCCGG'"

def test_reverse_complement_with_invalid_characters():
    # Test reverse complement with invalid characters in sequence
    with pytest.raises(ValueError, match="DNA sequence contains invalid characters."):
        reverse_complement("ATXGCC")  # 'X' is invalid
    with pytest.raises(ValueError, match="DNA sequence contains invalid characters."):
        reverse_complement("ATGGNC")  # 'N' is invalid

def test_reverse_complement_with_empty_sequence():
    # Test reverse complement with an empty string
    assert reverse_complement("") == "", "Failed to handle empty sequence"

def test_reverse_complement_with_lowercase_input():
    # Test lowercase sequence input
    assert reverse_complement("atgc".upper()) == "GCAT", "Failed to handle lowercase input 'atgc'"
    assert reverse_complement("aattccgg".upper()) == "CCGGAATT", "Failed to handle lowercase input 'aattccgg'"

def test_reverse_complement_with_palindromic_sequence():
    # Test palindromic sequences (reverse complement should be identical to original)
    assert reverse_complement("GATATC") == "GATATC", "Failed to handle palindromic sequence 'GATATC'"
