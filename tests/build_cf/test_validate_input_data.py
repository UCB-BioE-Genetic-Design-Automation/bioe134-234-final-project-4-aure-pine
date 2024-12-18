import pytest
from create_gb_and_cf.build_cf import validate_input_data

def test_validate_input_data_valid():
    """Test that valid input data passes without raising an error."""
    valid_data = {
        "steps": [],
        "sequences": {}
    }
    # Should not raise any exceptions
    validate_input_data(valid_data)

def test_validate_input_data_non_dict():
    """Test that non-dictionary input raises a ValueError."""
    invalid_data = []  # Not a dictionary
    with pytest.raises(ValueError, match="Input data must be a dictionary."):
        validate_input_data(invalid_data)

def test_validate_input_data_missing_steps():
    """Test that missing 'steps' key raises a ValueError."""
    invalid_data = {
        "sequences": {}
    }
    with pytest.raises(ValueError, match="Input data must include a 'steps' list."):
        validate_input_data(invalid_data)

def test_validate_input_data_steps_not_list():
    """Test that 'steps' key not being a list raises a ValueError."""
    invalid_data = {
        "steps": "not a list",
        "sequences": {}
    }
    with pytest.raises(ValueError, match="Input data must include a 'steps' list."):
        validate_input_data(invalid_data)

def test_validate_input_data_missing_sequences():
    """Test that missing 'sequences' key raises a ValueError."""
    invalid_data = {
        "steps": []
    }
    with pytest.raises(ValueError, match="Input data must include a 'sequences' dictionary."):
        validate_input_data(invalid_data)

def test_validate_input_data_sequences_not_dict():
    """Test that 'sequences' key not being a dictionary raises a ValueError."""
    invalid_data = {
        "steps": [],
        "sequences": "not a dictionary"
    }
    with pytest.raises(ValueError, match="Input data must include a 'sequences' dictionary."):
        validate_input_data(invalid_data)
