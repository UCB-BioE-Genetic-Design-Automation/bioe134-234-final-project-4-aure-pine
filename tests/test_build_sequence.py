import pytest
from Bio.Seq import Seq
from create_gb_and_cf.create_gb_and_cf import build_sequence

@pytest.fixture
def valid_input_data():
    """Provide valid input data for build_sequence."""
    return {
        "sequences": {
            "UTR5": "GCTA",
            "promoter": "GGTCA",
            "cds": "ATGGCTA",
            "UTR3": "AGGCT",
            "terminator": "TGCCT"
        },
        "overhangs": {
            "5_prime": "AATT",
            "3_prime": "GGCC"
        }
    }

@pytest.fixture
def invalid_input_data():
    """Provide invalid input data for build_sequence."""
    return {
        "sequences": {
            "UTR5": "GC123",
            "promoter": "GGTCA",
            "cds": "ATGGCTA",
            "UTR3": "AGGCT",
            "terminator": "TGCCT"
        },
        "overhangs": {
            "5_prime": "AATT",
            "3_prime": "GG12"
        }
    }

# Test Cases

def test_build_sequence_valid_inputs(valid_input_data):
    """Test building a sequence with valid components and overhangs."""
    sequence = build_sequence(valid_input_data)
    expected_sequence = Seq("AATTGCTAGGTCAATGGCTAAGGCTTGCCTGGCC")
    assert str(sequence) == str(expected_sequence), f"Expected {expected_sequence}, got {sequence}"

def test_build_sequence_missing_cds(valid_input_data):
    """Test behavior when 'cds' is missing from components."""
    del valid_input_data["sequences"]["cds"]
    with pytest.raises(ValueError, match="Component 'cds' is required but missing."):
        build_sequence(valid_input_data)

def test_build_sequence_invalid_bases(invalid_input_data):
    """Test building a sequence with invalid bases in components."""
    with pytest.raises(ValueError, match="Component 'UTR5' is missing a valid sequence or contains invalid characters"):
        build_sequence(invalid_input_data)

def test_build_sequence_invalid_overhangs(valid_input_data):
    """Test building a sequence with invalid bases in overhangs."""
    valid_input_data["overhangs"] = {
        "5_prime": "AATT",
        "3_prime": "GG12"
    }
    with pytest.raises(ValueError, match="3' overhang contains invalid characters"):
        build_sequence(valid_input_data)

