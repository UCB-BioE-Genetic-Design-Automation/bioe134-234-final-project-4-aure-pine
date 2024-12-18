from Bio.Seq import Seq
from create_gb_and_cf.build_gb import build_genbank_file
from create_gb_and_cf.build_cf import build_construction_file

VALID_ENZYMES = {
    "GoldenGate": ['AarI', 'BbsI', 'BsaI', 'BsmBI', 'SapI', 'BseRI'],
    "Gibson": [],  # Gibson does not use restriction enzymes.
    "RestrictionLigation": ['BamHI', 'BglII', 'EcoRI', 'XhoI', 'SpeI', 
                            'XbaI', 'PstI', 'HindIII', 'NotI', 'XmaI', 'SmaI', 'KpnI', 'SacI', 'SalI']
}

def validate_inputs(cloning_strategy: str, components: dict, restriction_enzymes: list[str] = None):
    """
    Description:
    Validates the inputs for the `create_gb_and_cf` function.

    Input:
    - cloning_strategy (str): Cloning strategy (e.g., "GoldenGate").
    - components (dict): Dictionary of component names mapped to DNA sequences.
    - restriction_enzymes (list[str], optional): List of restriction enzymes.

    Output:
    - tuple: Validated components dictionary and list of restriction enzymes.

    Raises:
    - ValueError: If the cloning strategy is unsupported.
    - ValueError: If the components dictionary is missing a CDS or contains invalid DNA characters.
    - ValueError: If any restriction enzyme is not valid for the cloning strategy.
    """
    # 1. Validate cloning strategy
    if cloning_strategy not in VALID_ENZYMES:
        raise ValueError(f"Unsupported cloning strategy: '{cloning_strategy}'. "
                         f"Valid options: {list(VALID_ENZYMES.keys())}")
    
    # 2. Validate components
    if "cds" not in components:
        raise ValueError("The 'components' dictionary must contain at least a 'cds' (coding sequence).")
    
    valid_bases = {"A", "T", "C", "G"}
    for name, seq in components.items():
        if not seq or not set(seq.upper()).issubset(valid_bases):
            raise ValueError(f"Component '{name}' contains invalid characters. Only A, T, C, G are allowed.")

    # 3. Validate restriction enzymes
    valid_enzymes_for_strategy = VALID_ENZYMES[cloning_strategy]
    validated_enzymes = []
    if restriction_enzymes:
        invalid_enzymes = [enz for enz in restriction_enzymes if enz not in valid_enzymes_for_strategy]
        if invalid_enzymes:
            raise ValueError(f"Invalid enzymes: {invalid_enzymes} for cloning strategy '{cloning_strategy}'. "
                             f"Allowed enzymes: {valid_enzymes_for_strategy}")
        validated_enzymes = restriction_enzymes

    return components, validated_enzymes

def build_sequence(components: dict, homologous_overhangs: dict = None):
    """
    Description:
    Constructs a DNA sequence by concatenating components and adding homologous overhangs.

    Input:
    - components (dict): Dictionary of component names mapped to DNA sequences.
    - homologous_overhangs (dict, optional): Dictionary of 5' and 3' overhang sequences.

    Output:
    - Seq: A Bio.Seq object representing the constructed sequence.

    Raises:
    - ValueError: If the components dictionary is missing a CDS.
    - ValueError: If any component or overhang contains invalid DNA characters.
    """
    valid_bases = {"A", "T", "C", "G"}

    # Ensure the components dictionary contains at least a CDS
    if "cds" not in components:
        raise ValueError("The 'components' dictionary must contain at least a 'cds' (coding sequence).")
    
    # Predefined order of components
    component_order = ["UTR5", "promoter", "cds", "UTR3", "terminator"]

    # Start with optional homologous 5' overhang
    sequence_parts = []
    if homologous_overhangs and "5_prime" in homologous_overhangs:
        overhang_5 = homologous_overhangs["5_prime"]
        if not set(overhang_5.upper()).issubset(valid_bases):
            raise ValueError("5' overhang contains invalid characters. Only A, T, C, and G are allowed.")
        sequence_parts.append(overhang_5.upper())

    # Add components in defined order if they exist
    for key in component_order:
        if key in components:
            seq = components[key].upper()
            if not set(seq).issubset(valid_bases):
                raise ValueError(f"Component '{key}' contains invalid characters. Only A, T, C, and G are allowed.")
            sequence_parts.append(seq)

    # Add optional homologous 3' overhang
    if homologous_overhangs and "3_prime" in homologous_overhangs:
        overhang_3 = homologous_overhangs["3_prime"]
        if not set(overhang_3.upper()).issubset(valid_bases):
            raise ValueError("3' overhang contains invalid characters. Only A, T, C, and G are allowed.")
        sequence_parts.append(overhang_3.upper())

    # Combine all parts into a single sequence
    full_sequence = "".join(sequence_parts)
    return Seq(full_sequence)

def create_gb_and_cf(construct_id, construct_name, cloning_strategy:str, plasmid:str, components:dict, restriction_enzymes:list[str], 
                     antibiotic:str = 'Amp', strain:str = 'S. cerevisiae' , gb_filename:str = 'gb_output', cf_filename:str = 'cf_output'):
    """
    Description:
    Creates both a GenBank file and a Construction File for a specified cloning strategy.

    Input:
    - construct_id (str): Identifier for the construct.
    - construct_name (str): Name of the construct.
    - cloning_strategy (str): Cloning strategy (e.g., "GoldenGate").
    - plasmid (str): Name of the backbone plasmid.
    - components (dict): Dictionary of components to annotate.
    - restriction_enzymes (list[str]): User-specified restriction enzymes.
    - antibiotic (str, optional): Antibiotic resistance marker (default: "Amp").
    - strain (str, optional): Strain name for the construction file (default: "S. cerevisiae").
    - gb_filename (str, optional): Desired name for the output GenBank file.
    - cf_filename (str, optional): Desired name for the output Construction File.

    Output:
    - tuple: The names of the created GenBank and Construction File.

    Raises:
    - ValueError: If inputs fail validation.
    """
    # Step 1: Validate inputs
    components, enzymes = validate_inputs(cloning_strategy, components, restriction_enzymes)

    # Step 2: Build the full sequence
    sequence = build_sequence(components)

    # Step 3: Build the GenBank file
    gb_file = build_genbank_file(construct_id, construct_name, sequence, components, cloning_strategy, gb_filename)

    # Step 4: Placeholder for Construction File creation
    cf_file = build_construction_file(components, cloning_strategy, enzymes, plasmid, antibiotic, strain, cf_filename)

    return gb_file, cf_file