import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from create_gb_and_cf.build_gb import annotate_features

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

# Test Cases

def test_annotate_features_standard_components(test_sequence, test_components):
    """
    Test annotation of standard components with valid mappings.
    """
    features = annotate_features(test_sequence, test_components)

    # Verify CDS annotation
    cds_feature = next((f for f in features if f.qualifiers.get("label") == "cds"), None)
    assert cds_feature is not None, "CDS feature is missing."
    assert cds_feature.type == "CDS", "CDS feature type is incorrect."
    assert cds_feature.location.start == test_sequence.find(test_components["cds"]), "CDS start position is incorrect."
    assert cds_feature.location.end == test_sequence.find(test_components["cds"]) + len(test_components["cds"]), "CDS end position is incorrect."

    # Verify Promoter annotation
    promoter_feature = next((f for f in features if f.qualifiers.get("label") == "promoter"), None)
    assert promoter_feature is not None, "Promoter feature is missing."
    assert promoter_feature.type == "promoter", "Promoter feature type is incorrect."
    assert promoter_feature.location.start == test_sequence.find(test_components["promoter"]), "Promoter start position is incorrect."
    assert promoter_feature.location.end == test_sequence.find(test_components["promoter"]) + len(test_components["promoter"]), "Promoter end position is incorrect."


def test_annotate_features_missing_components(test_sequence):
    """
    Test handling of components that are not found in the sequence.
    """
    invalid_components = {"cds": "INVALIDSEQ"}
    features = annotate_features(test_sequence, invalid_components)

    # Ensure no features are annotated
    assert len(features) == 0, "No features should be annotated for invalid components."


def test_annotate_features_custom_component_type(test_sequence):
    """
    Test annotation of components that do not map to predefined types.
    """
    custom_components = {"custom_element": "GAGTTCATGAATTCTTATA"}
    features = annotate_features(test_sequence, custom_components)

    # Verify custom component is annotated as misc_feature
    custom_feature = next((f for f in features if f.qualifiers.get("label") == "custom_element"), None)
    assert custom_feature is not None, "Custom element is missing."
    assert custom_feature.type == "misc_feature", "Custom elements should be annotated as misc_feature."


def test_annotate_features_combined_components_and_sites(test_sequence, test_components):
    """
    Test combined annotation of components and restriction sites.
    """
    features = annotate_features(test_sequence, test_components)

    # Verify total feature count matches components + restriction sites
    expected_feature_count = len(test_components)
    assert len(features) == expected_feature_count, "Total number of features does not match expectation."
