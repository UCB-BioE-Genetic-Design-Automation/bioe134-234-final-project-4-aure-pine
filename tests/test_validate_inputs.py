import pytest
from create_gb_and_cf.create_gb_and_cf import validate_inputs

# Fixtures
@pytest.fixture
def valid_components():
    return {
        "cds": "ATGGCTA",
        "promoter": "ATGCG",
        "terminator": "TCCCG"
    }

@pytest.fixture
def valid_cloning_strategy():
    return "GoldenGate"

@pytest.fixture
def valid_enzymes():
    return ["BsaI", "SapI"]

# Test Cases

def test_validate_inputs_valid_data(valid_components, valid_cloning_strategy, valid_enzymes):
    """
    Test if validate_inputs works with valid data.
    """
    components, enzymes = validate_inputs(valid_cloning_strategy, valid_components, valid_enzymes)
    assert components == valid_components, "Components validation failed with valid data."
    assert enzymes == valid_enzymes, "Enzyme validation failed with valid data."


def test_validate_inputs_missing_cds(valid_components, valid_cloning_strategy):
    """
    Test if validate_inputs raises an error when 'cds' is missing in components.
    """
    invalid_components = valid_components.copy()
    del invalid_components["cds"]
    with pytest.raises(ValueError, match="must contain at least a 'cds'"):
        validate_inputs(valid_cloning_strategy, invalid_components)


def test_validate_inputs_invalid_cloning_strategy(valid_components):
    """
    Test if validate_inputs raises an error for an invalid cloning strategy.
    """
    invalid_strategy = "InvalidStrategy"
    with pytest.raises(ValueError, match="Unsupported cloning strategy"):
        validate_inputs(invalid_strategy, valid_components)


def test_validate_inputs_invalid_component_sequence(valid_components, valid_cloning_strategy):
    """
    Test if validate_inputs raises an error for invalid characters in components.
    """
    invalid_components = valid_components.copy()
    invalid_components["cds"] = "ATGGCTA123"
    with pytest.raises(ValueError, match="contains invalid characters"):
        validate_inputs(valid_cloning_strategy, invalid_components)


def test_validate_inputs_invalid_enzymes(valid_components, valid_cloning_strategy):
    """
    Test if validate_inputs raises an error for invalid restriction enzymes.
    """
    invalid_enzymes = ["InvalidEnzyme"]
    with pytest.raises(ValueError, match="Invalid enzymes"):
        validate_inputs(valid_cloning_strategy, valid_components, invalid_enzymes)


def test_validate_inputs_no_enzymes(valid_components, valid_cloning_strategy):
    """
    Test if validate_inputs works when no restriction enzymes are provided.
    """
    components, enzymes = validate_inputs(valid_cloning_strategy, valid_components, None)
    assert components == valid_components, "Components validation failed when no enzymes provided."
    assert enzymes == [], "Enzyme validation failed when no enzymes provided."
