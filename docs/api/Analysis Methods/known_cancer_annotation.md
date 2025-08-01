### knownCancer

#### Short description

Annotates mutations with COSMIC and OncoKB cancer-related annotations to identify known cancer genes.

#### Signature

```python
def knownCancer(self, annotation_table, output_path=None, compress_output=True, join_column="Hugo_Symbol", oncokb_table=None, in_place=False):
```

#### Parameters

| Parameter         | Type           | Required | Description                                                                                           |
| ----------------- | -------------- | -------- | ----------------------------------------------------------------------------------------------------- |
| `annotation_table`| `str | Path`   | Yes      | Path to the COSMIC annotation table (.tsv or .tsv.gz format).                                        |
| `output_path`     | `str | Path`   | No       | Output file path. If not provided, saves with default naming convention.                             |
| `compress_output` | `bool`         | No       | Whether to compress the output file with gzip (default: True).                                       |
| `join_column`     | `str`          | No       | Column name to use for joining (default: "Hugo_Symbol").                                             |
| `oncokb_table`    | `str | Path`   | No       | Path to the OncoKB cancer gene list table (.tsv). Adds OncoKB annotations if provided.              |
| `in_place`        | `bool`         | No       | If True, replaces self.data with annotated data. If False, returns annotated DataFrame (default: False). |

#### Return value

Returns `pd.DataFrame` if `in_place=False`, containing the annotated mutation data with cancer gene annotations. Returns `None` if `in_place=True` and updates `self.data` directly.

#### Exceptions

List only those the user should handle:

* `FileNotFoundError`: if annotation files don't exist.
* `ValueError`: if join column is not found in DataFrame.

#### Minimal usage example

```python
>>> annotated_data = py_mut.knownCancer("cosmic_cancer_genes.tsv")
>>> print(f"Annotated {len(annotated_data)} mutations")
```