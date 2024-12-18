import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from create_gb_and_cf.build_gb import write_genbank_file
import os

# Fixtures
@pytest.fixture
def test_sequence():
    return Seq("ATGCGCAGCTGTTCAGAGTTCATGAATTCTTATAATGGTAATTTCAAAGCAAGATTGGGAGCACAATTGTCAAATCTATGATGAATTCTATTTTGTTGCGTTATACCCAGCTTTTGGGCTTTTTCGAGAACGATC")

@pytest.fixture
def test_features():
    return [
        SeqFeature(FeatureLocation(0, 6), type="promoter", qualifiers={"label": "promoter"}),
        SeqFeature(FeatureLocation(6, 20), type="CDS", qualifiers={"label": "cds"})
    ]

@pytest.fixture
def output_filename():
    return "test_genbank_output"

@pytest.fixture
def construct_name():
    return "TestConstruct"

@pytest.fixture
def cloning_strategy():
    return "GoldenGate"

@pytest.fixture
def test_id():
    return "TestID"

# Test Cases

def test_write_genbank_file_valid_input(test_id, construct_name, test_sequence, test_features, output_filename, cloning_strategy):
    """
    Test valid GenBank file creation.
    """
    genbank_filename = write_genbank_file(test_id, construct_name, test_sequence, test_features, output_filename, cloning_strategy)

    # Check that the file is created
    assert os.path.exists(genbank_filename), f"GenBank file {genbank_filename} was not created."

    # Validate content
    record = SeqIO.read(genbank_filename, "genbank")

    # Check sequence
    assert str(record.seq) == str(test_sequence), "The sequence in the GenBank file is incorrect."

    # Check construct name
    assert record.id == test_id, "Construct ID in the GenBank file is incorrect."
    assert record.name == construct_name, "Construct name in the GenBank file is incorrect."

    # Check features
    assert len(record.features) == len(test_features), "Number of features in the GenBank file does not match."
    for i, feature in enumerate(record.features):
        assert feature.type == test_features[i].type, f"Feature type mismatch at index {i}."
        assert feature.location.start == test_features[i].location.start, f"Feature start position mismatch at index {i}."
        assert feature.location.end == test_features[i].location.end, f"Feature end position mismatch at index {i}."

    # Cleanup
    os.remove(genbank_filename)


def test_write_genbank_file_empty_features(test_id, construct_name, test_sequence, output_filename, cloning_strategy):
    """
    Test GenBank file creation with no features.
    """
    genbank_filename = write_genbank_file(test_id, construct_name, test_sequence, [], output_filename, cloning_strategy)

    # Check that the file is created
    assert os.path.exists(genbank_filename), f"GenBank file {genbank_filename} was not created."

    # Validate content
    record = SeqIO.read(genbank_filename, "genbank")

    # Check sequence
    assert str(record.seq) == str(test_sequence), "The sequence in the GenBank file is incorrect."

    # Check construct name
    assert record.id == test_id, "Construct ID in the GenBank file is incorrect."
    assert record.name == construct_name, "Construct name in the GenBank file is incorrect."

    # Check features
    assert len(record.features) == 0, "GenBank file should not contain any features."

    # Cleanup
    os.remove(genbank_filename)


def test_write_genbank_file_invalid_filename(test_id, construct_name, test_sequence, test_features, cloning_strategy):
    """
    Test behavior with an invalid filename.
    """
    with pytest.raises(OSError, match="Invalid filename provided"):
        write_genbank_file(test_id, construct_name, test_sequence, test_features, "", cloning_strategy)
