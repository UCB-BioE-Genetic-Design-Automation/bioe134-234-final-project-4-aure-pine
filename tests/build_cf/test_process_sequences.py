import pytest
from create_gb_and_cf.build_cf import process_sequences

def test_process_sequences_valid():
    """Test that valid sequences are processed correctly."""
    input_sequences = {
        "Seq1": {
            "sequence": "ATCG",
            "ext5": "A",
            "ext3": "T",
            "is_double_stranded": True,
            "is_circular": False,
            "mod_ext5": "phosphate",
            "mod_ext3": "hydroxyl"
        },
        "Seq2": {
            "sequence": "GGCC",
            "ext5": None,
            "ext3": None,
            "is_double_stranded": False,
            "is_circular": True,
            "mod_ext5": None,
            "mod_ext3": None
        }
    }
    expected_output = {
        "Seq1": {
            "sequence": "ATCG",
            "ext5": "A",
            "ext3": "T",
            "is_double_stranded": True,
            "is_circular": False,
            "mod_ext5": "phosphate",
            "mod_ext3": "hydroxyl"
        },
        "Seq2": {
            "sequence": "GGCC",
            "ext5": None,
            "ext3": None,
            "is_double_stranded": False,
            "is_circular": True,
            "mod_ext5": None,
            "mod_ext3": None
        }
    }
    assert process_sequences(input_sequences) == expected_output

def test_process_sequences_missing_required_key():
    """Test that missing required keys in a sequence raises a ValueError."""
    input_sequences = {
        "Seq1": {
            "sequence": "ATCG",
            # Missing 'is_double_stranded'
            "is_circular": False
        }
    }
    with pytest.raises(ValueError, match="Sequence data missing required key: is_double_stranded"):
        process_sequences(input_sequences)

def test_process_sequences_empty():
    """Test that an empty sequence dictionary returns an empty dictionary."""
    input_sequences = {}
    assert process_sequences(input_sequences) == {}

def test_process_sequences_default_values():
    """Test that missing optional keys in a sequence get default values."""
    input_sequences = {
        "Seq1": {
            "sequence": "ATCG",
            "is_double_stranded": True,
            "is_circular": False
        }
    }
    expected_output = {
        "Seq1": {
            "sequence": "ATCG",
            "ext5": None,
            "ext3": None,
            "is_double_stranded": True,
            "is_circular": False,
            "mod_ext5": None,
            "mod_ext3": None
        }
    }
    assert process_sequences(input_sequences) == expected_output
