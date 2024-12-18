def build_construction_file(data):
    """
    Constructs a Construction File (CF) from the provided data.

    Parameters:
        data (dict): A dictionary containing the steps, sequences, and metadata.

    Returns:
        dict: A dictionary representing the constructed CF.
    """
    # Validate the structure of the input data
    validate_input_data(data)

    # Initialize the CF structure
    cf = {
        "steps": [],
        "sequences": {},
        "metadata": {}
    }

    # Process metadata
    cf["metadata"] = process_metadata(data.get("metadata", {}))

    # Process sequences
    cf["sequences"] = process_sequences(data.get("sequences", {}))

    # Process steps
    for step in data.get("steps", []):
        processed_step = process_step(step, cf["sequences"])
        cf["steps"].append(processed_step)

    return cf

def validate_input_data(data):
    """
    Validates the input data structure to ensure it meets the expected format.

    Parameters:
        data (dict): The input data to validate.

    Raises:
        ValueError: If the data structure is invalid.
    """
    if not isinstance(data, dict):
        raise ValueError("Input data must be a dictionary.")
    if "steps" not in data or not isinstance(data["steps"], list):
        raise ValueError("Input data must include a 'steps' list.")
    if "sequences" not in data or not isinstance(data["sequences"], dict):
        raise ValueError("Input data must include a 'sequences' dictionary.")

def process_metadata(metadata):
    """
    Processes the metadata section of the input data.

    Parameters:
        metadata (dict): The metadata dictionary.

    Returns:
        dict: Processed metadata.
    """
    default_metadata = {
        "experiment_name": "Unknown Experiment",
        "date": None,
        "author": None,
        "description": ""
    }
    merged_metadata = {**default_metadata, **metadata}
    return {key: merged_metadata[key] for key in default_metadata}

def process_sequences(sequences):
    """
    Processes the sequences section of the input data.

    Parameters:
        sequences (dict): The sequences dictionary.

    Returns:
        dict: Processed sequences.
    """
    processed_sequences = {}
    for name, seq_data in sequences.items():
        processed_sequences[name] = process_sequence(seq_data)
    return processed_sequences

def process_sequence(seq_data):
    """
    Processes an individual sequence.

    Parameters:
        seq_data (dict): The sequence data dictionary.

    Returns:
        dict: Processed sequence data.
    """
    required_keys = ["sequence", "is_double_stranded", "is_circular"]
    for key in required_keys:
        if key not in seq_data:
            raise ValueError(f"Sequence data missing required key: {key}")

    return {
        "sequence": seq_data.get("sequence", ""),
        "ext5": seq_data.get("ext5"),
        "ext3": seq_data.get("ext3"),
        "is_double_stranded": seq_data.get("is_double_stranded", False),
        "is_circular": seq_data.get("is_circular", False),
        "mod_ext5": seq_data.get("mod_ext5"),
        "mod_ext3": seq_data.get("mod_ext3")
    }

def process_step(step, sequences):
    """
    Processes an individual step in the CF.

    Parameters:
        step (dict): The step data dictionary.
        sequences (dict): The dictionary of sequences for validation.

    Returns:
        dict: Processed step data.
    """
    if "operation" not in step:
        raise ValueError("Step data missing 'operation' key.")

    operation = step["operation"]

    # Add operation-specific processing logic here if needed
    return step
