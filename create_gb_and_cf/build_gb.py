from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def annotate_features(seq: Seq, components: dict):
    """
    Description:
    Annotates sequence components on a DNA sequence with GenBank-compliant feature annotations.

    Input:
    - seq (Seq): A Bio.Seq object representing the DNA sequence to annotate.
    - components (dict): A dictionary of component names (e.g., "cds") mapped to their sequences (str).

    Output:
    - list: A list of Bio.SeqFeature objects representing the annotated features.

    Raises:
    - ValueError: If a component's sequence contains invalid DNA characters.
    """

    features = []
    
    # Map component names to GenBank standard types
    component_type_map = {
        "cds": "CDS",
        "promoter": "promoter",
        "terminator": "terminator",
        "utr5": "5'UTR",
        "utr3": "3'UTR",
    }

    # Annotate components
    for name, component_seq in components.items():
        start = seq.find(component_seq)
        if start != -1:
            end = start + len(component_seq)
            feature_type = component_type_map.get(name.lower(), "misc_feature")
            feature = SeqFeature(
                FeatureLocation(start, end),
                type=feature_type,
                qualifiers={"label": name, "sequence": component_seq}
            )
            features.append(feature)
        else:
            print(f"Warning: Component '{name}' not found in the sequence.")

    return features

def write_genbank_file(id, construct_name: str, seq: Seq, features: list, output_filename: str, cloning_strategy: str):
    """
    Description:
    Writes a GenBank file with annotated sequence features.

    Input:
    - id (str): Identifier for the sequence record.
    - construct_name (str): Name for the construct in the GenBank file.
    - seq (Seq): A Bio.Seq object representing the full sequence.
    - features (list): A list of Bio.SeqFeature objects for annotation.
    - output_filename (str): Base name for the GenBank file (without extension).
    - cloning_strategy (str): The cloning strategy used (e.g., "GoldenGate").

    Output:
    - str: The name of the generated GenBank file.

    Raises:
    - OSError: If the output filename is empty or invalid.
    """
    if not output_filename or not output_filename.strip():
        raise OSError("Invalid filename provided")

    # Create a SeqRecord with metadata
    record = SeqRecord(
        seq=seq,
        id=id,
        name=construct_name,
        description=f"Annotated construct for {construct_name} using {cloning_strategy}",
        annotations={"molecule_type": "DNA"}  # Add molecule_type annotation
    )
    record.features.extend(features)  # Add the annotated features

    # Define output filename
    genbank_filename = f"{output_filename}.gb"
    
    # Write the SeqRecord to a GenBank file
    SeqIO.write(record, genbank_filename, "genbank")
    print(f"GenBank file successfully written: {genbank_filename}")
    
    return genbank_filename

def build_genbank_file(id: str, construct_name: str, sequence: Seq, components: dict, cloning_strategy: str, output_filename: str):
    """
    Description:
    Builds and writes a GenBank file with annotated sequence components.

    Input:
    - id (str): Identifier for the sequence record.
    - construct_name (str): Name for the construct in the GenBank file.
    - sequence (Seq): A Bio.Seq object representing the DNA sequence.
    - components (dict): Dictionary of components to annotate.
    - cloning_strategy (str): The chosen cloning strategy.
    - output_filename (str): Base name for the GenBank file.

    Output:
    - str: The name of the generated GenBank file.

    Raises:
    - ValueError: If a component contains invalid DNA characters.
    - OSError: If the output filename is empty or invalid.
    """
    # Step 1: Annotate the features
    features = annotate_features(sequence, components)

    # Step 2: Write the GenBank file
    filename = write_genbank_file(id, construct_name, sequence, features, output_filename, cloning_strategy)
    
    return filename