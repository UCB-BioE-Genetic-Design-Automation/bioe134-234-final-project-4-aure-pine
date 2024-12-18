import copy
import pytest
from unittest.mock import patch
from create_gb_and_cf.create_gb_and_cf import create_gb_and_cf

@pytest.fixture
def valid_input_data():
    """Provide valid JSON-like input data for tests."""
    return {
        "metadata": {
            "experiment_name": "TestConstruct",
            "cloning_strategy": "GoldenGate",
            "antibiotic": "Amp",
            "strain": "S. cerevisiae",
            "gb_filename": "test_gb_output",
            "cf_filename": "test_cf_output"
        },
        "sequences": {
            "promoter": {"sequence": "ATGCGT", "is_double_stranded": True, "is_circular": False},
            "cds": {"sequence": "TTGACC", "is_double_stranded": True, "is_circular": False},
            "terminator": {"sequence": "GGCTAA", "is_double_stranded": True, "is_circular": False}
        },
        "steps": [
            {
                "operation": "Assemble",
                "input": ["promoter", "cds", "terminator"],
                "output": "TestPlasmid",
                "enzymes": ["BsaI"]
            }
        ]
    }

def test_create_gb_and_cf_valid(valid_input_data):
    """Test create_gb_and_cf with valid inputs."""
    with patch("create_gb_and_cf.build_gb.build_genbank_file") as mock_gb, \
     patch("create_gb_and_cf.build_cf.build_construction_file") as mock_cf:

        mock_gb.return_value = "test_gb_output.gb"
        mock_cf.return_value = "test_cf_output.cf"  # Return the filename

        gb_file, cf_file = create_gb_and_cf(valid_input_data)

        assert gb_file == "test_gb_output.gb"
        assert cf_file == "test_cf_output.cf"

def test_create_gb_and_cf_missing_cds(valid_input_data):
    """Test create_gb_and_cf raises error when 'cds' is missing in components."""
    invalid_input_data = valid_input_data.copy()
    del invalid_input_data["sequences"]["cds"]

    with pytest.raises(ValueError, match="The 'components' dictionary must contain at least a 'cds'"):
        create_gb_and_cf(invalid_input_data)

def test_create_gb_and_cf_invalid_enzymes(valid_input_data):
    """Test create_gb_and_cf raises error for invalid restriction enzymes."""
    invalid_input_data = copy.deepcopy(valid_input_data)
    invalid_input_data["steps"][0]["enzymes"] = ["InvalidEnzyme"]

    with pytest.raises(ValueError, match="Invalid enzymes: \\['InvalidEnzyme'\\] for cloning strategy 'GoldenGate'"):
        create_gb_and_cf(invalid_input_data)

def test_create_gb_and_cf_empty_inputs():
    """Test create_gb_and_cf raises error for empty inputs."""
    empty_input_data = {
        "metadata": {},
        "sequences": {},
        "steps": []
    }

    with pytest.raises(ValueError, match="Unsupported cloning strategy"):
        create_gb_and_cf(empty_input_data)


def test_create_gb_and_cf_default_file_names(valid_input_data):
    """Test create_gb_and_cf with default file names."""
    # Deep copy the fixture to prevent unintended modifications
    test_input_data = copy.deepcopy(valid_input_data)

    # Remove file name fields
    test_input_data["metadata"].pop("gb_filename", None)
    test_input_data["metadata"].pop("cf_filename", None)

    with patch("create_gb_and_cf.build_gb.build_genbank_file") as mock_gb, \
     patch("create_gb_and_cf.build_cf.build_construction_file") as mock_cf:

        mock_gb.return_value = "default_gb.gb"
        mock_cf.return_value = "construction_file_output.cf"

        gb_file, cf_file = create_gb_and_cf(test_input_data)

        assert gb_file == "default_gb.gb"
        assert cf_file == "default_cf.cf"


def test_create_gb_and_cf_unsupported_cloning_strategy(valid_input_data):
    """Test create_gb_and_cf raises error for unsupported cloning strategy."""
    invalid_input_data = valid_input_data.copy()
    invalid_input_data["metadata"]["cloning_strategy"] = "UnknownStrategy"

    with pytest.raises(ValueError, match="Unsupported cloning strategy: 'UnknownStrategy'"):
        create_gb_and_cf(invalid_input_data)

def test_create_gb_and_cf_invalid_dna_characters(valid_input_data):
    """Test create_gb_and_cf raises error for invalid DNA characters in components."""
    invalid_input_data = valid_input_data.copy()
    invalid_input_data["sequences"]["cds"]["sequence"] = "TTGACCX"  # Invalid character 'X'

    with pytest.raises(ValueError, match="Component 'cds' contains invalid characters"):
        create_gb_and_cf(invalid_input_data)
