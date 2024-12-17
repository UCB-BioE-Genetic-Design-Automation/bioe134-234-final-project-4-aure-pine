import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from create_gb_and_cf.build_gb import find_restriction_sites, annotate_features, write_genbank_file
from Bio import SeqIO
import os

# Fixtures
@pytest.fixture
def test_sequence():
    return Seq("ATGGCTAATGCGGAATTCCCTAGCGGAATTC")

@pytest.fixture
def test_components():
    return {
        "cds": "ATGGCTA",
        "promoter": "ATGCG",
        "terminator": "TCCCT"
    }

@pytest.fixture
def test_restriction_enzymes():
    return ["EcoRI", "BsaI"]  # EcoRI: GAATTC site exists twice

@pytest.fixture
def test_output_filename():
    return "test_output"

@pytest.fixture
def test_plasmid():
    return "pUC1345"

@pytest.fixture
def test_cloning_strategy():
    return "RestrictionLigation"

# Test find_restriction_sites
def test_find_restriction_sites_valid_enzymes(test_sequence, test_cloning_strategy, test_restriction_enzymes):
    """
    Test if restriction sites are found correctly for given enzymes.
    """
    sites = find_restriction_sites(test_sequence, test_cloning_strategy, test_restriction_enzymes)
    assert "EcoRI" in sites, "EcoRI should be found in the sequence."
    assert sites["EcoRI"] == [10, 25], "EcoRI sites are not at the expected positions."

def test_find_restriction_sites_no_enzymes(test_sequence, test_cloning_strategy):
    """
    Test behavior when no restriction enzymes are provided.
    """
    sites = find_restriction_sites(test_sequence, test_cloning_strategy, req_restr_enzymes=[])
    assert sites == {}, "No enzymes provided; restriction sites dictionary should be empty."

def test_find_restriction_sites_invalid_strategy(test_sequence):
    """
    Test behavior when an unsupported cloning strategy is provided.
    """
    sites = find_restriction_sites(test_sequence, "InvalidStrategy", req_restr_enzymes=["EcoRI"])
    assert "EcoRI" not in sites, "No sites should be found for an invalid strategy."

# Test annotate_features
def test_annotate_features_valid_input(test_sequence, test_components, test_restriction_enzymes):
    """
    Test if features (components and restriction sites) are annotated correctly.
    """
    # Simulate restriction sites
    restriction_sites = {"EcoRI": [10, 25]}

    features = annotate_features(test_sequence, test_components, restriction_sites, test_restriction_enzymes)

    # Verify component annotations
    cds_feature = next((f for f in features if f.qualifiers["label"] == "cds"), None)
    assert cds_feature is not None, "CDS feature is missing."
    assert cds_feature.location.start == 0 and cds_feature.location.end == 7, "CDS location is incorrect."

    # Verify restriction site annotations
    ecoRI_features = [f for f in features if f.qualifiers["label"] == "EcoRI"]
    assert len(ecoRI_features) == 2, "EcoRI restriction sites should be annotated twice."
    assert ecoRI_features[0].location.start == 10, "First EcoRI site location is incorrect."
    assert ecoRI_features[1].location.start == 25, "Second EcoRI site location is incorrect."

def test_annotate_features_no_restriction_sites(test_sequence, test_components):
    """
    Test annotation behavior when no restriction sites are provided.
    """
    features = annotate_features(test_sequence, test_components, restriction_sites={}, enzymes=[])
    restriction_site_features = [f for f in features if f.type == "restriction_site"]
    assert len(restriction_site_features) == 0, "There should be no restriction site annotations."

# Test write_genbank_file
def test_write_genbank_file_valid_input(test_sequence, test_output_filename, test_plasmid):
    """
    Test writing a valid GenBank file and its content.
    """
    # Create a mock list of features
    features = [
        SeqFeature(FeatureLocation(0, 7), type="CDS", qualifiers={"label": "cds"}),
        SeqFeature(FeatureLocation(10, 16), type="restriction_site", qualifiers={"label": "EcoRI"})
    ]

    # Write the file
    genbank_filename = write_genbank_file(test_sequence, features, test_output_filename, test_plasmid, "GoldenGate")

    # Check file creation
    assert os.path.exists(genbank_filename), f"GenBank file {genbank_filename} was not created."

    # Verify content
    record = SeqIO.read(genbank_filename, "genbank")
    assert len(record.features) == 2, "GenBank file should contain 2 features."
    assert record.features[0].type == "CDS", "First feature type should be 'CDS'."
    assert record.features[1].type == "restriction_site", "Second feature type should be 'restriction_site'."

    # Cleanup
    os.remove(genbank_filename)

def test_write_genbank_file_empty_features(test_sequence, test_output_filename, test_plasmid):
    """
    Test writing a GenBank file when no features are provided.
    """
    genbank_filename = write_genbank_file(test_sequence, [], test_output_filename, test_plasmid, "GoldenGate")

    # Check file creation
    assert os.path.exists(genbank_filename), f"GenBank file {genbank_filename} was not created."

    # Verify no features
    record = SeqIO.read(genbank_filename, "genbank")
    assert len(record.features) == 0, "GenBank file should have no features."

    # Cleanup
    os.remove(genbank_filename)