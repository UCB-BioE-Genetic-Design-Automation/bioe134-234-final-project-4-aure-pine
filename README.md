
# BioE 134 Final Project Submission

## Project Overview


## Scope of Work

---

## Function Descriptions

### 1. asdf (`asdf`)
- **Description**: This function takes a DNA sequence and returns its reverse complement. Only valid nucleotides (A, T, C, G) are allowed. The function raises a `ValueError` if invalid characters are found.
- **Input**: A string representing the DNA sequence.
- **Output**: A string representing the reverse complement of the input DNA sequence.

**Example**:
```python
reverse_complement("ATGC")
# Returns: "GCAT"
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
