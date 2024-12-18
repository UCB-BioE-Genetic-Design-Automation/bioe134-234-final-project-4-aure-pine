import pytest
from create_gb_and_cf.build_cf import build_construction_file

def test_build_cf_with_mock_data():
    """Test the build_cf function with realistic mock data."""
    mock_data = {
        "metadata": {
            "experiment_name": "Mock Experiment",
            "date": "2024-12-17",
            "author": "Ian",
            "description": "A mock test experiment."
        },
        "sequences": {
            "Oligo1": {
                "sequence": "ATCG",
                "ext5": None,
                "ext3": None,
                "is_double_stranded": False,
                "is_circular": False,
                "mod_ext5": "phosphate",
                "mod_ext3": "hydroxyl"
            },
            "TemplateDNA": {
                "sequence": "GGCC",
                "ext5": "G",
                "ext3": "C",
                "is_double_stranded": True,
                "is_circular": False,
                "mod_ext5": None,
                "mod_ext3": None
            }
        },
        "steps": [
            {
                "operation": "PCR",
                "output": "PCRProduct",
                "forward_oligo": "Oligo1",
                "reverse_oligo": "Oligo1",
                "template": "TemplateDNA",
                "product_size": 100
            },
            {
                "operation": "Digest",
                "output": "DigestFragment",
                "dna": "PCRProduct",
                "enzymes": ["EcoRI", "BamHI"],
                "fragSelect": 1
            },
            {
                "operation": "Ligate",
                "output": "LigatedProduct",
                "dnas": ["DigestFragment", "TemplateDNA"]
            }
        ]
    }

    expected_output = {
        "metadata": {
            "experiment_name": "Mock Experiment",
            "date": "2024-12-17",
            "author": "Ian",
            "description": "A mock test experiment."
        },
        "sequences": {
            "Oligo1": {
                "sequence": "ATCG",
                "ext5": None,
                "ext3": None,
                "is_double_stranded": False,
                "is_circular": False,
                "mod_ext5": "phosphate",
                "mod_ext3": "hydroxyl"
            },
            "TemplateDNA": {
                "sequence": "GGCC",
                "ext5": "G",
                "ext3": "C",
                "is_double_stranded": True,
                "is_circular": False,
                "mod_ext5": None,
                "mod_ext3": None
            }
        },
        "steps": [
            {
                "operation": "PCR",
                "output": "PCRProduct",
                "forward_oligo": "Oligo1",
                "reverse_oligo": "Oligo1",
                "template": "TemplateDNA",
                "product_size": 100
            },
            {
                "operation": "Digest",
                "output": "DigestFragment",
                "dna": "PCRProduct",
                "enzymes": ["EcoRI", "BamHI"],
                "fragSelect": 1
            },
            {
                "operation": "Ligate",
                "output": "LigatedProduct",
                "dnas": ["DigestFragment", "TemplateDNA"]
            }
        ]
    }

    filename = build_construction_file(mock_data)
    assert filename == "default_cf.cf"

def test_build_cf_empty_input():
    """Test build_cf with an empty input dictionary."""
    mock_data = {}
    with pytest.raises(ValueError, match="Input data must include a 'steps' list."):
        build_construction_file(mock_data)

def test_build_cf_missing_steps():
    """Test build_cf with missing steps key."""
    mock_data = {
        "metadata": {},
        "sequences": {}
    }
    with pytest.raises(ValueError, match="Input data must include a 'steps' list."):
        build_construction_file(mock_data)

def test_build_cf_steps_referencing_nonexistent_sequences():
    """Test build_cf with steps referencing non-existent sequences."""
    mock_data = {
        "metadata": {},
        "sequences": {
            "Oligo1": {
                "sequence": "ATCG",
                "is_double_stranded": False,
                "is_circular": False
            }
        },
        "steps": [
            {
                "operation": "PCR",
                "output": "PCRProduct",
                "forward_oligo": "NonExistentOligo",
                "reverse_oligo": "Oligo1",
                "template": "NonExistentTemplate",
                "product_size": 100
            }
        ]
    }
    result = build_construction_file(mock_data)
    # Validate that steps are included but external validation may flag missing sequences
    assert result["steps"] == mock_data["steps"]

def test_build_cf_large_data():
    """Test build_cf with a large number of sequences and steps."""
    sequences = {
        f"Seq{i}": {
            "sequence": "ATCG" * 10,
            "is_double_stranded": (i % 2 == 0),
            "is_circular": (i % 3 == 0)
        }
        for i in range(1, 101)
    }
    steps = [
        {
            "operation": "PCR",
            "output": f"Product{i}",
            "forward_oligo": f"Seq{i}",
            "reverse_oligo": f"Seq{i+1}",
            "template": f"Seq{i+2}",
            "product_size": 100 + i
        }
        for i in range(1, 51)
    ]
    mock_data = {
        "metadata": {
            "experiment_name": "Large Data Test"
        },
        "sequences": sequences,
        "steps": steps
    }
    result = build_construction_file(mock_data)
    assert len(result["sequences"]) == 100
    assert len(result["steps"]) == 50
