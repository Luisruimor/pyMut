### to_maf

#### Short description

Exports a PyMutation object back to MAF (Mutation Annotation Format) file format.

#### Signature

```python
def to_maf(self, output_path: str | Path) -> None:
```

#### Parameters

| Parameter     | Type         | Required | Description                                                                 |
| ------------- | ------------ | -------- | --------------------------------------------------------------------------- |
| `output_path` | `str | Path` | Yes      | Path where the MAF file will be written. Supports .maf and .maf.gz formats. |

#### Return value

Returns `None`. The method writes the MAF file to the specified path.

#### Exceptions

List only those the user should handle:

* `ValueError`: if the PyMutation object doesn't contain the necessary data for MAF export or if required VCF-style columns are missing.

#### Minimal usage example

```python
>>> from pyMut.input import read_maf
>>> py_mut = read_maf("input.maf")
>>> py_mut.to_maf("output.maf")
```