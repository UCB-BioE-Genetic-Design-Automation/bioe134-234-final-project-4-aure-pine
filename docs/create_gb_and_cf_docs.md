Main function to create GenBank and Construction File.
    1) Checks the inputs
        - These are if the component sequences are made out of the proper bases, the cloning strategy is one of the supported ones, 
        and if the optionally provided enzymes work with the provided cloning strategy
    2) Takes the components and builds a sequence
        - This includes all elements in the "components" dictionary and the homologous overhangs if they are given
    3) Builds a GenBank file (build_gb.py)
        1) Annotates the components, using the GenBank feature annotations
        2) Writing a GenBank file with the name construct_name
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