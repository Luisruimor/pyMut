# PyMutation - Main Class for Mutation Analysis

The **PyMutation** class is the central object of pyMut that encapsulates mutation data and provides methods for analysis and visualization.

## What is PyMutation?

PyMutation is the main class that represents a mutation dataset. It contains data in structured format, metadata about the data origin, and methods to perform analysis and generate visualizations.

## Object Structure

### Main Attributes

```python
class PyMutation:
    def __init__(self, data: pd.DataFrame, metadata: MutationMetadata, samples: List[str]):
        self.data = data           # DataFrame with mutations in VCF-like format
        self.samples = samples     # List of samples in the dataset
        self.metadata = metadata   # Metadata about origin and configuration
```

### `data` Attribute (pd.DataFrame)
Contains mutations in VCF-like format with standard columns:

```
CHROM | POS | ID | REF | ALT | QUAL | FILTER | SAMPLE_001 | SAMPLE_002 | ANNOTATIONS | ...
chr1  | 100 | .  | A   | G   | .    | .      | A|G        | A|A        | ...         | ...
chr2  | 200 | .  | C   | T   | .    | .      | C|C        | C|T        | ...         | ...
```

### `samples` Attribute (List[str])
List of sample identifiers:
```python
['SAMPLE_001', 'SAMPLE_002', 'SAMPLE_003', ...]
```

### `metadata` Attribute (MutationMetadata)
Information about origin and configuration:
```python
metadata.source_format    # "MAF" OR "VCF"
metadata.file_path       # Path to original file
metadata.loaded_at       # Loading timestamp
metadata.filters         # Applied filters
metadata.fasta          # Reference FASTA file
metadata.notes          # Comments from original file
```

## Creating PyMutation Objects

### From MAF file
```python
from pyMut.io import read_maf

# Recommended method
py_mut = read_maf("mutations.maf")
```

### From VCF file
Not implemented yet, but will be similar to MAF loading.