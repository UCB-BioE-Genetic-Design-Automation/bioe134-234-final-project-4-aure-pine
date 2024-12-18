import pytest
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, AllEnzymes
from create_gb_and_cf.build_gb import find_restriction_sites, VALID_ENZYMES

# Helper to retrieve enzymes as RestrictionBatch
def get_valid_enzymes(strategy):
    return RestrictionBatch([enzyme for enzyme in VALID_ENZYMES.get(strategy, []) if enzyme in AllEnzymes])

def build_seq():
    seq = ["ATGGCTAATGCGGAATT"]
    enzymes =[]
    for enz in VALID_ENZYMES.values():
        enzymes.extend(enz)

    rb = RestrictionBatch(enzymes)

    for enz in rb:
        seq.append(enz.site)

    seq.append("CCCTAGCGGAATTCCGGA")
    return ''.join(seq)

seq = build_seq()

# Fixtures
@pytest.fixture
def test_sequence():
    return Seq(seq)

@pytest.fixture
def no_enzymes():
    return []


# Test Cases
def test_find_restriction_sites_valid_golden_gate(test_sequence):
    """
    Test valid restriction site detection for GoldenGate strategy.
    """
    result = find_restriction_sites(test_sequence, "GoldenGate", VALID_ENZYMES["GoldenGate"])
    valid_enzymes = get_valid_enzymes("GoldenGate")

    assert isinstance(result, dict), "The result should be a dictionary."
    for enzyme in valid_enzymes:
        assert enzyme in result, f"{enzyme} should be in the result for GoldenGate enzymes."


def test_find_restriction_sites_valid_restriction_ligation(test_sequence):
    """
    Test restriction site detection for RestrictionLigation strategy.
    """
    result = find_restriction_sites(test_sequence, "RestrictionLigation", VALID_ENZYMES["RestrictionLigation"])
    valid_enzymes = get_valid_enzymes("RestrictionLigation")

    assert isinstance(result, dict), "The result should be a dictionary."
    for enzyme in valid_enzymes:
        assert enzyme in result, f"{enzyme} should be in the result for RestrictionLigation enzymes."

def test_find_restriction_sites_no_cut_sites():
    """
    Test behavior when the sequence does not contain any restriction sites for the default enzymes.
    """
    # Define a sequence that does not contain any known restriction sites
    no_cut_site_sequence = Seq("ATGGCTAATGCGGAAAGTCGACGAGTACCTAGGATCCAGTACCG")
    
    # Get restriction sites for GoldenGate enzymes
    result = find_restriction_sites(no_cut_site_sequence, "GoldenGate", [])
    valid_enzymes = get_valid_enzymes("GoldenGate")

    # Verify that no cut sites are found
    assert result == {}, "Result should be an empty dictionary when no restriction sites are found."


def test_find_restriction_sites_no_enzymes_provided(test_sequence):
    """
    Test behavior when no user-specified enzymes are provided.
    """
    result = find_restriction_sites(test_sequence, "GoldenGate", [])
    valid_enzymes = get_valid_enzymes("GoldenGate")

    # Verify that default GoldenGate enzymes are present
    for enzyme in valid_enzymes:
        assert enzyme in result, f"{enzyme} should be in the result as a default GoldenGate enzyme."


def test_find_restriction_sites_empty_valid_enzymes(test_sequence):
    """
    Test early exit when no valid enzymes are specified.
    """
    result = find_restriction_sites(test_sequence, "Gibson", [])  # Gibson has no valid enzymes
    assert result == {}, "Result should be empty when no valid enzymes are provided."


def test_find_restriction_sites_invalid_cloning_strategy(test_sequence):
    """
    Test behavior when an invalid cloning strategy is provided.
    """
    result = find_restriction_sites(test_sequence, "InvalidStrategy", [])
    assert result == {}, "Result should be empty for an unsupported cloning strategy."


def test_find_restriction_sites_duplicate_enzymes(test_sequence):
    """
    Test behavior when duplicate enzymes are provided.
    """
    duplicate_enzymes = ["BsaI", "BsaI", "SapI"]
    result = find_restriction_sites(test_sequence, "GoldenGate", duplicate_enzymes)
    valid_enzymes = get_valid_enzymes("GoldenGate")

    for enzyme in valid_enzymes:
        assert enzyme in result, f"{enzyme} should appear in the result."
    assert len(result) == len(valid_enzymes), "Result should not have duplicate enzymes."
