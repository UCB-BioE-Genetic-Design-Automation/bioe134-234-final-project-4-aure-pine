import json
from Bio.Seq import Seq
from create_gb_and_cf.build_gb import build_genbank_file
from create_gb_and_cf.build_cf import build_construction_file

VALID_ENZYMES = {
    "GoldenGate": ['AarI', 'BbsI', 'BsaI', 'BsmBI', 'SapI', 'BseRI'],
    "Gibson": [],  # Gibson does not use restriction enzymes.
    "RestrictionLigation": ['BamHI', 'BglII', 'EcoRI', 'XhoI', 'SpeI', 
                            'XbaI', 'PstI', 'HindIII', 'NotI', 'XmaI', 'SmaI', 'KpnI', 'SacI', 'SalI']
}

def create_gb_and_cf(input_data: dict):
    """
    Creates both a GenBank file and a Construction File from input data.

    Parameters:
    - input_data (dict): JSON-like dictionary containing metadata, sequences, and steps.

    Returns:
    - tuple: Names of the created GenBank and Construction Files.
    """
    # Step 1: Validate inputs
    components, enzymes = validate_inputs(input_data)

    # Step 2: Build the full sequence
    sequence = build_sequence(input_data["sequences"])

    # Preprocess components to extract raw sequences
    raw_components = {key: value["sequence"] for key, value in input_data["sequences"].items()}

    # Step 3: Build the GenBank file
    print("Calling build_genbank_file...")  # Debugging statement
    gb_file = build_genbank_file(
        id=input_data["metadata"]["experiment_name"],
        construct_name=input_data["metadata"]["experiment_name"],
        sequence=sequence,
        components=raw_components,  # Use raw sequences
        cloning_strategy=input_data["metadata"].get("cloning_strategy", ""),
        output_filename=input_data["metadata"].get("gb_filename", "default_gb")
    )

    # Step 4: Build the Construction File
    # cf_file = build_construction_file(input_data)

    return gb_file, #cf_file

def validate_inputs(input_data: dict):
    """
    Validates the input data structure for the 'create_gb_and_cf' function.

    Parameters:
    - input_data (dict): The JSON-like dictionary containing input metadata, sequences, and steps.

    Returns:
    - tuple: Validated components dictionary and restriction enzymes list.

    Raises:
    - ValueError: For unsupported cloning strategies or invalid DNA characters.
    """
    # Extract metadata and sequences
    metadata = input_data.get("metadata", {})
    sequences = input_data.get("sequences", {})

    # Validate cloning strategy
    cloning_strategy = metadata.get("cloning_strategy", "")
    if cloning_strategy not in VALID_ENZYMES:
        raise ValueError(f"Unsupported cloning strategy: '{cloning_strategy}'. "
                         f"Valid options: {list(VALID_ENZYMES.keys())}")

    # Validate components
    if "cds" not in sequences:
        raise ValueError("The 'sequences' dictionary must contain at least a 'cds' (coding sequence).")

    valid_bases = {"A", "T", "C", "G"}
    for name, seq_dict in sequences.items():
        seq = seq_dict['sequence']
        if not seq or not set(seq.upper()).issubset(valid_bases):
            raise ValueError(f"Component '{name}' contains invalid characters. Only A, T, C, G are allowed.")

    # Validate restriction enzymes
    restriction_enzymes = metadata.get("enzymes", []) or []  # Ensure None is replaced with an empty list
    valid_enzymes = VALID_ENZYMES[cloning_strategy]
    invalid_enzymes = [enz for enz in restriction_enzymes if enz not in valid_enzymes]
    if invalid_enzymes:
        raise ValueError(f"Invalid enzymes: {invalid_enzymes} for cloning strategy '{cloning_strategy}'. "
                         f"Allowed enzymes: {valid_enzymes}")

    return sequences, restriction_enzymes


def build_sequence(input_data: dict):
    """
    Constructs a DNA sequence by concatenating components and adding homologous overhangs.

    Parameters:
    - input_data (dict): Dictionary containing sequences and optional overhangs.

    Returns:
    - Seq: A Bio.Seq object representing the constructed sequence.
    """
    sequences = input_data.get("sequences", {})
    homologous_overhangs = input_data.get("overhangs", {})
    valid_bases = {"A", "T", "C", "G"}
    sequence_parts = []

    # Ensure 'cds' is present
    if "cds" not in sequences:
        raise ValueError("Component 'cds' is required but missing.")

    # Start with optional homologous 5' overhang
    if "5_prime" in homologous_overhangs:
        overhang_5 = homologous_overhangs["5_prime"]
        if not set(overhang_5.upper()).issubset(valid_bases):
            raise ValueError("5' overhang contains invalid characters. Only A, T, C, and G are allowed.")
        sequence_parts.append(overhang_5.upper())

    # Process each sequence component
    for name, seq in sequences.items():
        if not seq or not set(seq.upper()).issubset(valid_bases):
            raise ValueError(f"Component '{name}' is missing a valid sequence or contains invalid characters.")
        sequence_parts.append(seq.upper())

    # Add optional homologous 3' overhang
    if "3_prime" in homologous_overhangs:
        overhang_3 = homologous_overhangs["3_prime"]
        if not set(overhang_3.upper()).issubset(valid_bases):
            raise ValueError("3' overhang contains invalid characters. Only A, T, C, and G are allowed.")
        sequence_parts.append(overhang_3.upper())

    # Combine all parts into a single sequence
    full_sequence = "".join(sequence_parts)
    return Seq(full_sequence)

