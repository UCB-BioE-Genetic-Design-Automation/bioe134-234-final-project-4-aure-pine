from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import RestrictionBatch, Analysis, RestrictionType
from Bio import SeqIO

VALID_ENZYMES = {
    "GoldenGate": ['AarI', 'BbsI', 'BsaI', 'BsmBI', 'SapI', 'BseRI'],
    "Gibson": [],  # Gibson does not use restriction enzymes.
    "RestrictionLigation": ['BamHI', 'BglII', 'EcoRI', 'XhoI', 'SpeI', 
                            'XbaI', 'PstI', 'HindIII', 'NotI', 'XmaI', 'SmaI', 'KpnI', 'SacI', 'SalI']
}

def find_restriction_sites(seq: Seq, cloning_strategy: str, req_restr_enzymes: list[str] = None):
    """
    Finds restriction enzyme cut sites on a sequence based on the cloning strategy and/or user-specified enzymes.

    Parameters:
    - seq (Seq): The sequence to analyze for restriction enzyme sites.
    - cloning_strategy (str): The chosen cloning strategy (e.g., "GoldenGate", "Gibson").
    - req_restr_enzymes (list[str], optional): User-specified restriction enzymes.

    Returns:
    - dict: A dictionary mapping enzymes to their cut site positions on the sequence.
    """
    # Combine strategy-specific enzymes and user-defined enzymes
    strategy_enzymes = VALID_ENZYMES.get(cloning_strategy, [])
    all_enzymes = set(strategy_enzymes + (req_restr_enzymes or []))
    
    if not all_enzymes:
        print("No restriction enzymes specified for this strategy. Skipping restriction site analysis.")
        return {}

    # Create a RestrictionBatch for the combined enzymes
    restriction_batch = RestrictionBatch(list(all_enzymes))

    # Analyze the sequence for cut sites
    analysis = Analysis(restriction_batch, seq)
    cut_sites = analysis.full()  # Dictionary mapping enzyme names to cut positions

    return cut_sites

def annotate_features(seq: Seq, components: dict, restriction_sites: dict, enzymes: list[str]):
    """
    Annotates components and restriction sites on a sequence.

    Parameters:
    - seq (Seq): The sequence to annotate.
    - components (dict): Dictionary of components with their sequences (e.g., {"cds": "ATGCGT"}).
    - restriction_sites (dict): Dictionary of enzymes to their cut site positions.
    - enzymes (list[str]): List of enzyme names to dynamically retrieve restriction site lengths.

    Returns:
    - list: List of SeqFeature annotations for components and restriction sites.
    """
    from Bio.Restriction import RestrictionBatch

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

    # Annotate restriction sites with variable lengths
    restriction_batch = RestrictionBatch(enzymes)
    for enzyme_name, sites in restriction_sites.items():
        enzyme = restriction_batch.get(enzyme_name)
        if enzyme and isinstance(enzyme, RestrictionType):
            site_length = len(enzyme)  # Dynamically retrieve site length
        else:
            site_length = 6  # Default to 6 bp if enzyme site length is unavailable

        for site in sites:
            end = site + site_length
            feature = SeqFeature(
                FeatureLocation(site, end),
                type="restriction_site",
                qualifiers={"label": enzyme_name}
            )
            features.append(feature)

    return features

def write_genbank_file(seq: Seq, features: list, output_filename: str, plasmid: str, cloning_strategy: str):
    """
    Writes a GenBank file with annotated sequence features.

    Parameters:
    - seq (Seq): The full sequence to write to the GenBank file.
    - features (list): List of SeqFeature annotations for components and restriction sites.
    - output_filename (str): Base name for the output GenBank file.
    - plasmid (str): Name of the backbone plasmid (used for metadata).
    - cloning_strategy (str): The chosen cloning strategy (used for metadata).

    Returns:
    - str: The name of the generated GenBank file.
    """
    # Create a SeqRecord with metadata
    record = SeqRecord(
        seq=seq,
        id="constructed_sequence",
        name=output_filename,
        description=f"Annotated construct for cloning strategy: {cloning_strategy}, Backbone: {plasmid}"
    )
    record.features = features  # Add the annotated features

    # Define output filename
    genbank_filename = f"{output_filename}.gb"
    
    # Write the SeqRecord to a GenBank file
    SeqIO.write(record, genbank_filename, "genbank")
    print(f"GenBank file successfully written: {genbank_filename}")
    
    return genbank_filename

def build_genbank_file(sequence: Seq, components: dict, cloning_strategy: str, req_restr_enzymes: list[str], plasmid: str, output_filename: str):
    """
    Builds and writes a GenBank file with annotated components and restriction sites.

    Parameters:
    - sequence (Seq): The built sequence from components.
    - components (dict): Dictionary of components with their sequences.
    - cloning_strategy (str): Cloning strategy to determine valid restriction enzymes.
    - req_restr_enzymes (list[str]): Validated restriction enzymes.
    - output_filename (str): Name for the output GenBank file.
    - plasmid (str): Name of the backbone plasmid for metadata.
    """
    # Step 1: Find restriction sites
    restriction_sites = find_restriction_sites(sequence, cloning_strategy, req_restr_enzymes)

    # Step 2: Annotate the features
    features = annotate_features(sequence, components, restriction_sites, req_restr_enzymes)

    # Step 3: Write the GenBank file
    filename = write_genbank_file(sequence, features, output_filename, plasmid, cloning_strategy)
    
    return filename