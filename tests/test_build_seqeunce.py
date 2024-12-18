import pytest
from Bio.Seq import Seq
from create_gb_and_cf.create_gb_and_cf import build_sequence

# Fixtures
@pytest.fixture
def valid_components():
    return {
        "UTR5": "GCTA",
        "promoter": "GGTCA",
        "cds": "ATGGCTA",
        "UTR3": "AGGCT",
        "terminator": "TGCCT"
    }

@pytest.fixture
def valid_overhangs():
    return {
        "5_prime": "AATT",
        "3_prime": "GGCC"
    }

@pytest.fixture
def empty_overhangs():
    return {}

@pytest.fixture
def invalid_components():
    return {
        "UTR5": "GC123",
        "promoter": "GGTCA",
        "cds": "ATGGCTA",
        "UTR3": "AGGCT",
        "terminator": "TGCCT"
    }

# Test Cases

def test_build_sequence_valid_inputs(valid_components, valid_overhangs):
    """
    Test building a sequence with valid components and overhangs.
    """
    sequence = build_sequence(valid_components, valid_overhangs)
    expected_sequence = Seq("AATTGCTAGGTCAATGGCTAAGGCTTGCCTGGCC")  # 5' + components + 3'
    assert str(sequence) == str(expected_sequence), f"Expected {expected_sequence}, got {sequence}"


def test_build_sequence_no_overhangs(valid_components, empty_overhangs):
    """
    Test building a sequence without homologous overhangs.
    """
    sequence = build_sequence(valid_components, empty_overhangs)
    expected_sequence = Seq("GCTAGGTCAATGGCTAAGGCTTGCCT")  # components only
    assert str(sequence) == str(expected_sequence), f"Expected {expected_sequence}, got {sequence}"


def test_build_sequence_missing_cds(valid_components):
    """
    Test behavior when 'cds' is missing from components.
    """
    invalid_components = valid_components.copy()
    del invalid_components["cds"]
    with pytest.raises(ValueError, match="cds"):
        build_sequence(invalid_components, {})


def test_build_sequence_invalid_bases(invalid_components, valid_overhangs):
    """
    Test building a sequence with invalid bases in components.
    """
    with pytest.raises(ValueError, match="contains invalid characters"):
        build_sequence(invalid_components, valid_overhangs)


def test_build_sequence_invalid_overhangs(valid_components):
    """
    Test building a sequence with invalid bases in overhangs.
    """
    invalid_overhangs = {
        "5_prime": "AATT",
        "3_prime": "GG12"
    }
    with pytest.raises(ValueError, match="overhang contains invalid characters"):
        build_sequence(valid_components, invalid_overhangs)
