from Bio.Seq import Seq
from build_gb import build_genbank_file
from build_cf import build_construction_file

VALID_ENZYMES = {
    "GoldenGate": ['AarI', 'BbsI', 'BsaI', 'BsmBI', 'SapI', 'BseRI'],
    "Gibson": [],  # Gibson does not use restriction enzymes.
    "RestrictionLigation": ['BamHI', 'BglII', 'EcoRI', 'XhoI', 'SpeI', 
                            'XbaI', 'PstI', 'HindIII', 'NotI', 'XmaI', 'SmaI', 'KpnI', 'SacI', 'SalI']
}

def validate_inputs(cloning_strategy: str, components: dict, req_restr_enzymes: list[str] = None):
    """
    Validates the inputs for the create_gb_and_cf function.

    Checks:
    1. Cloning strategy is supported.
    2. Components contain only valid DNA bases (A, T, C, G) and must include at least a CDS.
    3. Restriction enzymes (if provided) are valid for the cloning strategy.

    Parameters:
    - cloning_strategy (str): Cloning strategy (e.g., "GoldenGate", "Gibson").
    - components (dict): Dictionary of component names to DNA sequences.
    - req_restr_enzymes (list[str], optional): User-defined restriction enzymes.

    Returns:
    - dict: Cleaned and validated components.
    - list[str]: Validated restriction enzymes.
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
    if req_restr_enzymes:
        invalid_enzymes = [enz for enz in req_restr_enzymes if enz not in valid_enzymes_for_strategy]
        if invalid_enzymes:
            raise ValueError(f"Invalid enzymes: {invalid_enzymes} for cloning strategy '{cloning_strategy}'. "
                             f"Allowed enzymes: {valid_enzymes_for_strategy}")
        validated_enzymes = req_restr_enzymes

    return components, validated_enzymes

def build_sequence(components: dict, homologous_overhangs: dict = None):
    """
    Builds a sequence by concatenating components and adding homologous overhangs.

    Parameters:
    - components (dict): Dictionary of components with their sequences.
        Must include at least 'cds'. Example:
        {"cds": "ATGGCTA", "promoter": "GGTCA", "UTR5": "GCTA", "terminator": "TGCCT"}
    - homologous_overhangs (dict, optional): Dictionary of overhangs to add to the sequence.
        Example: {"5_prime": "AATT", "3_prime": "GGCC"}

    Returns:
    - Seq: A Bio.Seq object representing the built sequence.
    """
    valid_bases = {"A", "T", "C", "G"}

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

def create_gb_and_cf(cloning_strategy:str, plasmid:str, components:dict, req_restr_enzymes:list[str] = None, antibiotic:str = 'Amp', gb_filename:str = 'gb_output', cf_filename:str = 'cf_output'):
    """
    Main function to create GenBank and Construction File.
    1) Checks the inputs
        - These are if the component sequences are made out of the proper bases, the cloning strategy is one of the supported ones, 
        and if the optionally provided enzymes work with the provided cloning strategy
    2) Takes the components and builds a sequence
        - This includes all elements in the "components" dictionary and the homologous overhangs if they are given
    3) Builds a GenBank file
        1) Finds restriction sites on the built sequence based on the cloning strategy and/or given required restriction enzymes (req_restr_enzymes)
        2) Annotates the components, using the GenBank feature annotations
        3) Annotates restriction sites using the length of the binding sites for the restriction enzymes and the proper GenBank feature annotation categories.
        4) Writing a GenBank file with the name from output_filename (if given, otherwise just 'output')
    4) Creates a ConstructionFile text file
        - To be implemented and thought out later. 

    Parameters:
    - cloning_strategy (str): Cloning strategy (e.g., "GoldenGate", "Gibson").
    - plasmid (str): Name of the backbone plasmid. This is NOT a sequence for the backbone, but a name eg. "pUC1345". Used for 	building the Construction file.
    - components (dict): Dictionary of components to annotate, e.g., 
        {"cds": "GCTA", "promoter": "ATCG", "UTR5": "GATC", "terminator": "CTAG", "UTR3": "GCTA"}. These are mappings of specific sequence features to 
        sequence strings of A,T,G, or C characters. It must contain at least a CDS.
    - req_restr_enzymes (list[str], optional): List of user given restriction enzymes that they want to use in the experiment. 
        These are validated to work with the desired cloning strategy by comparing them to the restriction enzyme lists.
    - output_filename (str, optional): The desired name for the output GenBank file.
    """
    # Step 1: Validate inputs
    components, enzymes = validate_inputs(cloning_strategy, components, req_restr_enzymes)

    # Step 2: Build the full sequence
    sequence = build_sequence(components)

    # Step 3: Build the GenBank file
    gb_file = build_genbank_file(sequence, components, cloning_strategy, enzymes, plasmid, gb_filename)

    # Step 4: Placeholder for Construction File creation
    cf_file = build_construction_file(components, cloning_strategy, enzymes, plasmid, antibiotic, cf_filename)

    return gb_file, cf_file