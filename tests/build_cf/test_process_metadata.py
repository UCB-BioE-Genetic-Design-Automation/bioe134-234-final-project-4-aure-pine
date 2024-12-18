import pytest
from create_gb_and_cf.build_cf import process_metadata

def test_process_metadata_valid():
    """Test processing valid metadata with all fields provided."""
    input_metadata = {
        "experiment_name": "Test Experiment",
        "date": "2024-12-17",
        "author": "Ian",
        "description": "A test experiment."
    }
    expected_output = {
        "experiment_name": "Test Experiment",
        "date": "2024-12-17",
        "author": "Ian",
        "description": "A test experiment."
    }
    assert process_metadata(input_metadata) == expected_output

def test_process_metadata_missing_fields():
    """Test that missing fields in metadata are filled with default values."""
    input_metadata = {
        "experiment_name": "Partial Metadata"
    }
    expected_output = {
        "experiment_name": "Partial Metadata",
        "date": None,
        "author": None,
        "description": ""
    }
    assert process_metadata(input_metadata) == expected_output

def test_process_metadata_empty():
    """Test that an empty metadata dictionary returns only default values."""
    input_metadata = {}
    expected_output = {
        "experiment_name": "Unknown Experiment",
        "date": None,
        "author": None,
        "description": ""
    }
    assert process_metadata(input_metadata) == expected_output

def test_process_metadata_extra_keys():
    """Test that extra keys in the metadata are ignored."""
    input_metadata = {
        "experiment_name": "Extra Keys Test",
        "extra_key": "ExtraValue"
    }
    expected_output = {
        "experiment_name": "Extra Keys Test",
        "date": None,
        "author": None,
        "description": ""
    }
    assert process_metadata(input_metadata) == expected_output
