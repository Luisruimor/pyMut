### read_maf

#### Short description

Reads a MAF file and returns a PyMutation object with automatic caching and performance optimizations.

#### Signature

```python
def read_maf(path: str | Path, assembly: str, cache_dir: Optional[str | Path] = None, consolidate_variants: bool = True) -> PyMutation:
```

#### Parameters

| Parameter | Type | Required | Description |
| --------- | ---- | -------- | ----------- |
| `path` | `str \| Path` | Yes | Path to the MAF file (.maf or .maf.gz). |
| `assembly` | `str` | Yes | Version of the genome assembly. Must be either "37" or "38". |
| `cache_dir` | `str \| Path` | No | Directory for caching processed files. If None, uses a .pymut_cache directory next to the input file. |
| `consolidate_variants` | `bool` | No | If True (default), consolidates identical variants that appear in multiple samples into a single row. If False, maintains the original behavior where each MAF row becomes a separate DataFrame row. |

#### Return value

Returns a `PyMutation` object containing the mutations read from the MAF file, converted to wide-format. Includes metadata with information about the source file, comments, and configuration.

#### Exceptions

* `FileNotFoundError`: if the specified MAF file path does not exist.
* `ValueError`: if the MAF file is missing required columns or has an invalid format, or if the assembly parameter is not "37" or "38".
* `ImportError`: if the 'pyarrow' library cannot be imported and the 'c' engine alternative also fails.
* `Exception`: for any other errors encountered while reading or processing the file.

#### Minimal usage example

```python
>>> from pyMut import read_maf
>>> mutations = read_maf("data/mutations.maf", assembly="38")
>>> print(mutations.data.shape)
(1000, 25)
```