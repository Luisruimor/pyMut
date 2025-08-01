### to_vcf

#### Short description

Exports a PyMutation object to VCF format with proper headers and metadata information.

#### Signature

```python
def to_vcf(self, output_path: str | Path) -> None:
```

#### Parameters

| Parameter     | Type         | Required | Description                                                                 |
| ------------- | ------------ | -------- | --------------------------------------------------------------------------- |
| `output_path` | `str | Path` | Yes      | Path where the VCF file will be written. Can be a string or Path object.   |

#### Return value

Returns `None`. The method writes the VCF file to the specified output path.

#### Exceptions

List only those the user should handle:

* `ValueError`: if the PyMutation object doesn't contain the necessary data for VCF export (missing required VCF-style columns like CHROM, POS, REF, ALT, ID or missing sample columns).

#### Minimal usage example

```python
>>> pymutation_obj.to_vcf("output.vcf")
>>> # Or using Path object
>>> from pathlib import Path
>>> pymutation_obj.to_vcf(Path("data/output.vcf"))
```