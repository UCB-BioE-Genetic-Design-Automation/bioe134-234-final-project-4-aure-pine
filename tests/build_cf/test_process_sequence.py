import pytest
from create_gb_and_cf.build_cf import process_sequence

def test_process_sequence_valid():
    """Test processing a valid sequence."""
    input_sequence = {
        "sequence": "ATCG",
        "ext5": "A",
        "ext3": "T",
        "is_double_stranded": True,
        "is_circular": False,
        "mod_ext5": "phosphate",
        "mod_ext3": "hydroxyl"
    }
    expected_output = {
        "sequence": "ATCG",
        "ext5": "A",
        "ext3": "T",
        "is_double_stranded": True,
        "is_circular": False,
        "mod_ext5": "phosphate",
        "mod_ext3": "hydroxyl"
    }
    assert process_sequence(input_sequence) == expected_output

def test_process_sequence_missing_required_key():
    """Test that missing a required key raises a ValueError."""
    input_sequence = {
        "sequence": "ATCG",
        # Missing 'is_double_stranded'
        "is_circular": False
    }
    with pytest.raises(ValueError, match="Sequence data missing required key: is_double_stranded"):
        process_sequence(input_sequence)

def test_process_sequence_defaults():
    """Test that optional keys get default values when missing."""
    input_sequence = {
        "sequence": "ATCG",
        "is_double_stranded": True,
        "is_circular": False
    }
    expected_output = {
        "sequence": "ATCG",
        "ext5": None,
        "ext3": None,
        "is_double_stranded": True,
        "is_circular": False,
        "mod_ext5": None,
        "mod_ext3": None
    }
    assert process_sequence(input_sequence) == expected_output

def test_process_sequence_empty_sequence():
    """Test that an empty sequence raises a ValueError."""
    input_sequence = {
        "sequence": "",
        "is_double_stranded": True,
        "is_circular": False
    }
    expected_output = {
        "sequence": "",
        "ext5": None,
        "ext3": None,
        "is_double_stranded": True,
        "is_circular": False,
        "mod_ext5": None,
        "mod_ext3": None
    }
    assert process_sequence(input_sequence) == expected_output

def test_process_sequence_extra_keys():
    """Test that extra keys in the input are ignored."""
    input_sequence = {
        "sequence": "ATCG",
        "is_double_stranded": True,
        "is_circular": False,
        "extra_key": "ExtraValue"
    }
    expected_output = {
        "sequence": "ATCG",
        "ext5": None,
        "ext3": None,
        "is_double_stranded": True,
        "is_circular": False,
        "mod_ext5": None,
        "mod_ext3": None
    }
    assert process_sequence(input_sequence) == expected_output
