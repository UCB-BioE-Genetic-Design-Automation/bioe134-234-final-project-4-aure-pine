import pytest
from create_gb_and_cf.create_gb_and_cf import validate_inputs

@pytest.fixture
def valid_input_data():
    """Provide valid input data for validate_inputs."""
    return {
        "metadata": {
            "cloning_strategy": "GoldenGate",
            "enzymes": ["BsaI", "SapI"]
        },
        "sequences": {
            "cds": "ATGGCTA",
            "promoter": "ATGCG",
            "terminator": "TCCCG"
        }
    }

# Test Cases

def test_validate_inputs_valid_data(valid_input_data):
    """
    Test if validate_inputs works with valid data.
    """
    components, enzymes = validate_inputs(valid_input_data)
    assert components == valid_input_data["sequences"], "Components validation failed with valid data."
    assert enzymes == valid_input_data["metadata"]["enzymes"], "Enzyme validation failed with valid data."

def test_validate_inputs_missing_cds(valid_input_data):
    """
    Test if validate_inputs raises an error when 'cds' is missing in components.
    """
    del valid_input_data["sequences"]["cds"]
    with pytest.raises(ValueError, match="must contain at least a 'cds'"):
        validate_inputs(valid_input_data)

def test_validate_inputs_invalid_cloning_strategy(valid_input_data):
    """
    Test if validate_inputs raises an error for an invalid cloning strategy.
    """
    valid_input_data["metadata"]["cloning_strategy"] = "InvalidStrategy"
    with pytest.raises(ValueError, match="Unsupported cloning strategy"):
        validate_inputs(valid_input_data)

def test_validate_inputs_invalid_component_sequence(valid_input_data):
    """
    Test if validate_inputs raises an error for invalid characters in components.
    """
    valid_input_data["sequences"]["cds"] = "ATGGCTA123"
    with pytest.raises(ValueError, match="contains invalid characters"):
        validate_inputs(valid_input_data)

def test_validate_inputs_invalid_enzymes(valid_input_data):
    """
    Test if validate_inputs raises an error for invalid restriction enzymes.
    """
    valid_input_data["metadata"]["enzymes"] = ["InvalidEnzyme"]
    with pytest.raises(ValueError, match="Invalid enzymes"):
        validate_inputs(valid_input_data)

def test_validate_inputs_no_enzymes(valid_input_data):
    """
    Test if validate_inputs works when no restriction enzymes are provided.
    """
    valid_input_data["metadata"]["enzymes"] = None
    components, enzymes = validate_inputs(valid_input_data)
    assert components == valid_input_data["sequences"], "Components validation failed when no enzymes provided."
    assert enzymes == [], "Enzyme validation failed when no enzymes provided."