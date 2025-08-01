### calculate_tmb_analysis

#### Short description

Calculate Tumor Mutation Burden (TMB) analysis for each sample in a PyMutation object.

#### Signature

```python
def calculate_tmb_analysis(self, variant_classification_column: Optional[str] = None, genome_size_bp: int = 60456963, output_dir: str = ".", save_files: bool = True) -> Dict[str, pd.DataFrame]:
```

#### Parameters

| Parameter | Type | Required | Description |
| --------- | ---- | -------- | ----------- |
| `variant_classification_column` | `Optional[str]` | No | Name of the column containing variant classification information. If None, will automatically detect variant classification columns. |
| `genome_size_bp` | `int` | No | Size of the interrogated region in base pairs for TMB normalization. Default 60,456,963 bp (WES standard). Use ~3,000,000,000 bp for WGS. |
| `output_dir` | `str` | No | Directory where output files will be saved. Default is current directory. |
| `save_files` | `bool` | No | Whether to save the results to TSV files. Default is True. |

#### Return value

Returns a dictionary with two DataFrames: `{'analysis': pd.DataFrame, 'statistics': pd.DataFrame}`. The 'analysis' DataFrame contains per-sample TMB metrics, and the 'statistics' DataFrame contains global TMB statistics (mean, median, quartiles, etc.).

#### Exceptions

* `ValueError`: if PyMutation object is invalid (missing 'data' or 'samples' attributes).
* `ValueError`: if PyMutation data is empty.
* `ValueError`: if no samples found in PyMutation object.
* `ValueError`: if missing required columns (REF, ALT) in PyMutation data.
* `ValueError`: if provided variant classification column is not found in data.

#### Minimal usage example

```python
>>> from pyMut.io import read_maf
>>> py_mut = read_maf("mutations.maf")
>>> tmb_results = py_mut.calculate_tmb_analysis()
>>> print(tmb_results['analysis'].head())
```