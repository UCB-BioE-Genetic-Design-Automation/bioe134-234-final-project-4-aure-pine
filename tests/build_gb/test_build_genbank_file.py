import pytest
from Bio.Seq import Seq
from Bio import SeqIO
from create_gb_and_cf.build_gb import build_genbank_file
import os

# Fixtures
@pytest.fixture
def test_sequence():
    return Seq("ATGCGCAGCTGTTCAGAGTTCATGAATTCTTATAATGGTAATTTCAAAGCAAGATTGGGAGCACAATTGTCAAATCTATGATGAATTCTATTTTGTTGCGTTATACCCAGCTTTTGGGCTTTTTCGAGAACGATCTCCCT")

@pytest.fixture
def test_components():
    return {
        "cds": "CAGCTGTTCAGAGTTCATGAATTCTTATAATGGTAATTTCAAAGCAAGATTGGGAGCACAATTGTCAAATCTATGATGAATTCTATTTTGTTGCGTTATACCCAGCTTTTGGGCTTTTTCGAGAACGATC",
        "promoter": "ATGCG",
        "terminator": "TCCCT"
    }

@pytest.fixture
def output_filename():
    return "test_genbank_file"

@pytest.fixture
def id():
    return "TestID"

@pytest.fixture
def construct_name():
    return "TestConstruct"

@pytest.fixture
def cloning_strategy():
    return "GoldenGate"

# Test Cases

def test_build_genbank_file_valid_input(id, construct_name, test_sequence, test_components, cloning_strategy, output_filename):
    """
    Test valid GenBank file creation with annotations.
    """
    genbank_filename = build_genbank_file(id, construct_name, test_sequence, test_components, cloning_strategy, output_filename)

    # Check file creation
    assert os.path.exists(genbank_filename), f"GenBank file {genbank_filename} was not created."

    # Validate file content
    record = SeqIO.read(genbank_filename, "genbank")
    assert str(record.seq) == str(test_sequence), "The sequence in the GenBank file is incorrect."
    assert record.id == id, "ID in the GenBank file is incorrect."
    assert record.name == construct_name, "Construct name in the GenBank file is incorrect."
    assert len(record.features) == len(test_components), "Number of features in the GenBank file does not match components."

    # Cleanup
    os.remove(genbank_filename)


def test_build_genbank_file_missing_components(id, construct_name, test_sequence, cloning_strategy, output_filename):
    """
    Test behavior when some components are not found in the sequence.
    """
    components = {"cds": "INVALIDSEQ"}  # Component not in the sequence
    genbank_filename = build_genbank_file(id, construct_name, test_sequence, components, cloning_strategy, output_filename)

    # Validate file content
    record = SeqIO.read(genbank_filename, "genbank")
    assert len(record.features) == 0, "No features should be annotated for missing components."

    # Cleanup
    os.remove(genbank_filename)


def test_build_genbank_file_empty_components(id, construct_name, test_sequence, cloning_strategy, output_filename):
    """
    Test behavior when no components are provided.
    """
    genbank_filename = build_genbank_file(id, construct_name, test_sequence, {}, cloning_strategy, output_filename)

    # Validate file content
    record = SeqIO.read(genbank_filename, "genbank")
    assert len(record.features) == 0, "GenBank file should not contain any features."

    # Cleanup
    os.remove(genbank_filename)


def test_build_genbank_file_invalid_output_filename(id, construct_name, test_sequence, test_components, cloning_strategy):
    """
    Test behavior when an invalid output filename is provided.
    """
    with pytest.raises(OSError, match="Invalid filename provided"):
        build_genbank_file(id, construct_name, test_sequence, test_components, cloning_strategy, "")


def test_build_genbank_file_complex_sequence(id, construct_name, cloning_strategy, output_filename):
    """
    Test behavior with complex sequences and overlapping components.
    """
    sequence = Seq("ATGCATGCATGCATGC")
    components = {
        "cds1": "ATGCATGC",
        "cds2": "GCATGCAT"
    }
    genbank_filename = build_genbank_file(id, construct_name, sequence, components, cloning_strategy, output_filename)

    # Validate file content
    record = SeqIO.read(genbank_filename, "genbank")
    assert len(record.features) == len(components), "Number of features does not match components."
    os.remove(genbank_filename)
