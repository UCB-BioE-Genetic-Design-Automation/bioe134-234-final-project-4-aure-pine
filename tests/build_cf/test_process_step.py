import pytest
from create_gb_and_cf.build_cf import process_step

def test_process_step_valid_pcr():
    """Test processing a valid PCR step."""
    step = {
        "operation": "PCR",
        "output": "Product1",
        "forward_oligo": "Oligo1",
        "reverse_oligo": "Oligo2",
        "template": "TemplateDNA",
        "product_size": 500
    }
    sequences = {
        "Oligo1": {},
        "Oligo2": {},
        "TemplateDNA": {}
    }
    assert process_step(step, sequences) == step

def test_process_step_missing_operation():
    """Test that a step missing 'operation' raises a ValueError."""
    step = {
        "output": "Product1"
    }
    sequences = {}
    with pytest.raises(ValueError, match="Step data missing 'operation' key."):
        process_step(step, sequences)

def test_process_step_invalid_operation():
    """Test that an invalid operation type is returned as is (if allowed)."""
    step = {
        "operation": "InvalidOperation",
        "output": "Product1"
    }
    sequences = {}
    # If process_step doesn't validate operation, it should pass as-is
    assert process_step(step, sequences) == step

def test_process_step_missing_input_sequence():
    """Test that missing input sequences raises no error but can validate externally."""
    step = {
        "operation": "PCR",
        "output": "Product1",
        "forward_oligo": "Oligo1",
        "reverse_oligo": "Oligo2",
        "template": "MissingTemplate"
    }
    sequences = {
        "Oligo1": {},
        "Oligo2": {}
    }
    # No exception should be raised unless specifically checking for sequence presence
    assert process_step(step, sequences) == step

def test_process_step_extra_keys():
    """Test that extra keys in a step are preserved."""
    step = {
        "operation": "Digest",
        "output": "Fragment1",
        "dna": "InputDNA",
        "enzymes": ["EcoRI", "BamHI"],
        "extra_key": "ExtraValue"
    }
    sequences = {
        "InputDNA": {}
    }
    assert process_step(step, sequences) == step
