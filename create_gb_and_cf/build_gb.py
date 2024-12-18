from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def annotate_features(seq: Seq, components: dict):
    """
    Annotates components on a sequence.

    Parameters:
    - seq (Seq): The sequence to annotate.
    - components (dict): Dictionary of components with their sequences (e.g., {"cds": "ATGCGT"}).

    Returns:
    - list: List of SeqFeature annotations for components.
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
    Writes a GenBank file with annotated sequence features.

    Parameters:
    - id (str):
    - construct_name (str):
    - seq (Seq): The full sequence to write to the GenBank file.
    - features (list): List of SeqFeature annotations for components and restriction sites.
    - output_filename (str): Base name for the output GenBank file.
    - plasmid (str): Name of the backbone plasmid (used for metadata).
    - cloning_strategy (str): The chosen cloning strategy (used for metadata).

    Returns:
    - str: The name of the generated GenBank file.
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
    Builds and writes a GenBank file of the created insert with annotated components.

    Parameters:
    - sequence (Seq): The built sequence from components.
    - components (dict): Dictionary of components with their sequences.
    - cloning_strategy (str): Cloning strategy to determine valid restriction enzymes.
    - plasmid (str): Name of the backbone plasmid for metadata.
    - output_filename (str): Name for the output GenBank file.
    """
    # Step 1: Annotate the features
    features = annotate_features(sequence, components)

    # Step 2: Write the GenBank file
    filename = write_genbank_file(id, construct_name, sequence, features, output_filename, cloning_strategy)
    
    return filename