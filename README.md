# DNA Cloning Utilities Project

## Project Overview
My group and I worked on the Prototyping and Optimization project. We implemented methods to streamline pathway or circuit design in Yeast.

---

## Scope of Work
What I specifically worked on was the generation of Construction and GenBank files from a JSON description. This was intended to automate the creation of annotated GenBank and Construction files for DNA cloning experiments. The project supports multiple cloning strategies, including GoldenGate, Gibson, and Restriction-Ligation, offering robust input validation and compatibility checks.

I used JSON as the input so an LLM could theoretically parse user input and send this suite of functions data in easily usable format. 

The goal of this project is to streamline the process of constructing annotated DNA sequences for cloning experiments. The key functionalities include:

1. **Annotation of Sequence Features**: Identifies and annotates components such as CDS, promoters, and UTRs on a given DNA sequence.
2. **GenBank File Generation**: Creates GenBank files with annotated sequence features.
3. **Construction File Creation**: Generates construction files for DNA assembly protocols.
4. **Validation of Cloning Inputs**: Ensures that user-provided DNA components, cloning strategies, and restriction enzymes meet the specified criteria.

Each function is implemented in **Python** with error handling and input validation to ensure compatibility and accuracy.

---

## Function Descriptions

### 1. Annotate Features (`annotate_features`)
- **Description**: Annotates DNA components on a sequence with GenBank-compliant feature annotations.
- **Input**:
  - `seq` (string): The DNA sequence to annotate.
  - `components` (dict): Dictionary mapping component names (e.g., "cds", "promoter") to their sequences.
- **Output**: A list of annotated features.

---

### 2. Write GenBank File (`write_genbank_file`)
- **Description**: Generates a GenBank file with annotated features.
- **Input**:
  - `id` (string): Identifier for the sequence.
  - `construct_name` (string): Name of the construct.
  - `seq` (string): DNA sequence of the construct.
  - `features` (list): Annotated features for the GenBank file.
  - `output_filename` (string): Desired GenBank file name.
  - `cloning_strategy` (string): Cloning strategy used (e.g., "GoldenGate").
- **Output**: Name of the generated GenBank file.

---

### 3. Build GenBank File (`build_genbank_file`)
- **Description**: Combines sequence building, annotation, and file writing into one step for creating annotated GenBank files.
- **Input**:
  - `id` (string): Sequence identifier.
  - `construct_name` (string): Construct name.
  - `sequence` (string): Built DNA sequence.
  - `components` (dict): Dictionary of components to annotate.
  - `cloning_strategy` (string): Cloning strategy.
  - `output_filename` (string): Output GenBank file name.
- **Output**: Name of the generated GenBank file.

---

### 4. Validate Inputs (`validate_inputs`)
- **Description**: Validates cloning inputs, including components, cloning strategies, and restriction enzymes.
- **Input**:
  - `cloning_strategy` (string): Cloning strategy to validate.
  - `components` (dict): Dictionary of DNA components.
  - `restriction_enzymes` (list, optional): Restriction enzymes to validate.
- **Output**: Validated components and restriction enzymes.

---

### 5. Create GenBank and Construction File (`create_gb_and_cf`)
- **Description**: Validates inputs, builds a DNA sequence, and generates both GenBank and Construction files.
- **Input**:
  - `construct_id` (string): Unique construct identifier.
  - `construct_name` (string): Construct name.
  - `cloning_strategy` (string): Cloning strategy.
  - `plasmid` (string): Backbone plasmid name.
  - `components` (dict): Dictionary of DNA components.
  - `restriction_enzymes` (list, optional): User-specified restriction enzymes.
  - `antibiotic` (string, optional): Antibiotic resistance marker (default: "Amp").
  - `strain` (string, optional): Host strain (default: "S. cerevisiae").
  - `gb_filename` (string, optional): GenBank file name.
  - `cf_filename` (string, optional): Construction file name.
- **Output**: Filenames of the generated GenBank and Construction files.

---

## Error Handling

### Annotate Features
- Raises `ValueError` if the `components` dictionary contains sequences with invalid DNA characters (non-A, T, C, G).
- Logs warnings if a component is not found in the provided sequence.

### Write GenBank File
- Raises `OSError` if the `output_filename` is empty or invalid.

### Build GenBank File
- Propagates errors from `annotate_features` and `write_genbank_file`.

### Validate Inputs
- Raises `ValueError` if:
  - An unsupported cloning strategy is provided.
  - The `components` dictionary does not include a CDS.
  - Component sequences contain invalid DNA characters.
  - Restriction enzymes provided are not valid for the chosen cloning strategy.

### Build Sequence
- Raises `ValueError` if the `components` dictionary does not contain a CDS.
- Raises `ValueError` if any sequence in the `components` or `homologous_overhangs` contains invalid DNA characters.


### Create GenBank and Construction File
- Propagates errors from `validate_inputs`, `build_sequence`, and `build_genbank_file`.

---

## Testing

The functions have been thoroughly tested for:
- Valid DNA sequences and components.
- Missing or invalid components.
- Unsupported cloning strategies.
- Restriction enzyme compatibility.

Testing was performed using **pytest**, with comprehensive test cases covering standard, edge, and error scenarios.

---

## Usage Instructions

Clone the repository and install the required dependencies listed in `requirements.txt`. Functions can be imported from the respective Python modules.

**Example**:
```bash
pip install -r requirements.txt
```

```python
from create_gb_and_cf import create_gb_and_cf

components = {
    "promoter": "GGTCA",
    "cds": "ATGGCTA",
    "terminator": "TGCCT"
}

result = create_gb_and_cf(
    construct_id="1",
    construct_name="ExampleConstruct",
    cloning_strategy="GoldenGate",
    plasmid="pUC19",
    components=components,
    restriction_enzymes=["BsaI"]
)

print("GenBank File:", result[0])
print("Construction File:", result[1])
```

## Conclusion
This project provides a robust pipeline for generating annotated GenBank and Construction files for various cloning strategies, streamlining DNA sequence assembly and validation for bioengineering applications.