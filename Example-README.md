
# BioE 134 Final Project Submission

## Project Overview

This project provides two core bioinformatics utilities: 

1. **Reverse Complement (revcomp)**: Calculates the reverse complement of a DNA sequence.
2. **Translate**: Translates a DNA sequence into a protein sequence according to the standard genetic code.

These functions are implemented in **Python** and are part of the broader bioinformatics toolset aimed at automating genetic sequence analysis tasks.

---

## Scope of Work

As part of the final project for BioE 134, I developed two functions that are foundational for sequence analysis:

1. **Reverse Complement**: This function returns the reverse complement of a DNA sequence, which is an essential task in many genetic analysis pipelines.
   
2. **Translate**: This function translates a DNA sequence into a corresponding protein sequence by converting each codon into its corresponding amino acid, based on the standard genetic code.

Both functions include input validation and error handling to ensure proper use. The reverse complement function raises an error for sequences containing invalid characters, while the translate function raises an error for sequences not divisible by three, as well as sequences containing invalid characters.

---

## Function Descriptions

### 1. Reverse Complement (`reverse_complement`)

- **Description**: This function takes a DNA sequence and returns its reverse complement. Only valid nucleotides (A, T, C, G) are allowed. The function raises a `ValueError` if invalid characters are found.
- **Input**: A string representing the DNA sequence.
- **Output**: A string representing the reverse complement of the input DNA sequence.

**Example**:
```python
reverse_complement("ATGC")
# Returns: "GCAT"
```

### 2. Translate (`translate`)

- **Description**: This function translates a DNA sequence into a corresponding protein sequence. The input sequence must be divisible by 3. If it contains invalid characters or is not a multiple of three, the function raises a `ValueError`. Stop codons are represented as underscores (`_`).
- **Input**: A string representing the DNA sequence.
- **Output**: A string representing the translated protein sequence.

**Example**:
```python
translate("ATGGCC")
# Returns: "MA"
```

---

## Error Handling

### Reverse Complement
- Raises `ValueError` if invalid characters (anything other than A, T, C, G) are present in the DNA sequence.

### Translate
- Raises `ValueError` if the sequence contains invalid characters or if the sequence length is not a multiple of three.

---

## Testing

Both functions have been tested with standard, edge, and invalid input cases. A comprehensive suite of tests has been implemented using **pytest**.

- **Test File**: `tests/test_bio_functions.py`

The tests include:
- Valid sequences
- Sequences containing invalid characters
- Sequences with lengths not divisible by three (for the translate function)
- Palindromic sequences (for reverse complement)
- Lowercase input handling

---

## Usage Instructions

Clone the repository and install the required dependencies listed in `requirements.txt`. The functions can be imported from the `bio_functions.py` module.

**Example**:

```bash
pip install -r requirements.txt
```

Once installed, you can use the functions as follows:

```python
from bio_functions import reverse_complement, translate

# Example DNA sequence
dna_sequence = "ATGC"

# Reverse complement
print(reverse_complement(dna_sequence))

# Translate
print(translate("ATGGCC"))
```

---

## Conclusion

These two functions provide foundational operations for working with DNA sequences in bioinformatics pipelines. They have been tested and documented, ensuring proper error handling and robust functionality.
