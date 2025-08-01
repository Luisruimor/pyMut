### read_vcf

#### Short description

High-performance reader that converts a VCF (or VCF-GZ) file into a `PyMutation` object, with PyArrow acceleration, automatic caching and optional Tabix indexing.&#x20;

#### Signature

```python
def read_vcf(
    path: str | pathlib.Path,
    assembly: str,
    create_index: bool = False,
    cache_dir: str | pathlib.Path | None = None
) -> PyMutation:
```

#### Parameters

| Parameter      | Type                  | Required | Description                                                                                                   |
| -------------- | --------------------- | -------- | ------------------------------------------------------------------------------------------------------------- |
| `path`         | `str \| Path`         | **Yes**  | Path to the VCF (plain or `*.vcf.gz`).                                                                        |
| `assembly`     | `str`                 | **Yes**  | Genome build identifier, **must** be `"37"` or `"38"`.                                                        |
| `create_index` | `bool`                | No       | If `True`, create a Tabix (`.tbi`) index if missing (requires `tabix` in `PATH`). Default `False`.            |
| `cache_dir`    | `str \| Path \| None` | No       | Directory where parsed Parquet caches are stored. `None` (default) writes next to the VCF in `.pymut_cache/`. |

#### Return value

`PyMutation` — a wide-format table of variants plus metadata and sample columns, ready for downstream analysis.

#### Exceptions

* `FileNotFoundError` – VCF file does not exist.
* `ValueError` – invalid assembly, missing required VCF columns, or header problems.
* `Exception` – any other I/O or parsing error (e.g. broken compression, PyArrow failure).

#### Minimal usage example

```python
from pymutation.io import read_vcf

pymut = read_vcf(
    "tumour.vcf.gz",
    assembly="38",
    create_index=True          # build Tabix if needed
)

print(pymut.data.shape)
```
### Standard VCF Columns

```
CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO | FORMAT | SAMPLE_001 | SAMPLE_002
chr1  | 100 | .  | A   | G   | 60   | PASS   | ...  | GT:DP  | 0/1:30    | 1/1:25
```

### Conversion to **pyMut** Format

```
CHROM | POS | ID | REF | ALT | QUAL | FILTER | SAMPLE_001 | SAMPLE_002 | INFO_parsed
chr1  | 100 | .  | A   | G   | 60   | PASS   | A|G        | G|G        | {...}
```

## Complete Example

```python
from pyMut.input import read_vcf
import logging

# Enable logging to monitor progress
logging.basicConfig(level=logging.INFO)

# Load a VCF file with all options enabled
py_mut = read_vcf(
    path="src/pyMut/data/examples/ALL.chr10.vcf.gz",
    fasta="reference/hg38.fasta",
    create_index=True,
    cache_dir="cache/"
)

# Verify that the file was loaded successfully
print(f"Loaded samples: {len(py_mut.samples)}")
print(f"Total variants: {len(py_mut.data)}")
print(f"Unique chromosomes: {py_mut.data['CHROM'].unique()}")

# Genotype information
print(f"Sample columns: {py_mut.samples[:5]}...")  # First 5 samples

# Check metadata
print(f"Source format: {py_mut.metadata.source_format}")
print(f"FASTA file: {py_mut.metadata.fasta}")
```

## Genotype Handling

The function automatically converts VCF genotypes into pyMut’s allelic format:

### VCF Genotypes (input)

```
FORMAT: GT:DP:GQ
SAMPLE_001: 0/1:30:99    # Heterozygous
SAMPLE_002: 1/1:25:99    # Homozygous alternate
SAMPLE_003: 0/0:35:99    # Homozygous reference
```

### pyMut Format (output)

```
SAMPLE_001: A|G    # REF|ALT
SAMPLE_002: G|G    # ALT|ALT  
SAMPLE_003: A|A    # REF|REF
```

## Cache System

```python
# First load — cache is created
py_mut1 = read_vcf("large_file.vcf.gz", cache_dir="cache/")

# Second load — cache is used (much faster)
py_mut2 = read_vcf("large_file.vcf.gz", cache_dir="cache/")
```
