### annotate_pfam

#### Short description

Annotate PyMutation data with PFAM domains using database annotation or VEP data extraction.

#### Signature

```python
def annotate_pfam(self, db_conn: Optional[duckdb.DuckDBPyConnection] = None, *, aa_column: str = 'aa_pos', auto_extract: bool = True, prefer_database: bool = True):
```

#### Parameters

| Parameter | Type | Required | Description |
| --------- | ---- | -------- | ----------- |
| `db_conn` | `Optional[duckdb.DuckDBPyConnection]` | No | DuckDB connection for PFAM database. If None, will create one automatically. |
| `aa_column` | `str` | No | Name of the column containing amino acid positions. Default is 'aa_pos'. |
| `auto_extract` | `bool` | No | If True, automatically extract uniprot/aa_pos from VEP data when missing. Default is True. |
| `prefer_database` | `bool` | No | If True, prefer database annotation over VEP parsing when both are available. Default is True. |

#### Return value

Returns a new PyMutation object with PFAM domain annotations added as additional columns. Can return `None` if annotation fails.

#### Exceptions

* `PfamAnnotationError`: if database connection fails or annotation process encounters errors.
* `ValueError`: if required columns are missing or data format is invalid.

#### Minimal usage example

```python
>>> from pyMut.io import read_maf
>>> py_mut = read_maf("mutations_vep_annotated.maf")
>>> annotated = py_mut.annotate_pfam()
>>> print(annotated.data.columns)
```