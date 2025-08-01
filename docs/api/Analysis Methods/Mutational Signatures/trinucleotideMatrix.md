### trinucleotideMatrix

#### Short description

Builds the **96-trinucleotide context matrix** for all single-nucleotide variants (SNVs) in a `PyMutation` object, enriching the original data with context annotations to enable downstream signature analysis.

#### Signature

```python
def trinucleotideMatrix(
    self,
    fasta_file: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
```

#### Parameters

| Parameter    | Type         | Required | Description                                                                                       |
| ------------ | ------------ | -------- | ------------------------------------------------------------------------------------------------- |
| `fasta_file` | `str` (path) | **Yes**  | Path to the reference-genome FASTA file used to look up trinucleotide contexts (must be indexed). |

#### Return value

A tuple `(contexts_df, enriched_data)`:

* **`contexts_df`** – `pd.DataFrame` of shape **96 × N-samples**. Each row corresponds to one of the 96 canonical trinucleotide mutation classes; each column contains the raw counts for a sample.
* **`enriched_data`** – original mutation table with three extra columns:

  * `trinuc` – the reference trinucleotide (e.g. `"ACA"`).
  * `class96` – class label in the form `"A[C>T]A"`.
  * `idx96` – integer index (0–95) into the standard context order.

Both DataFrames are never `None`.

#### Exceptions

* `ImportError` – *pyfaidx* is missing.
* `ValueError` – required columns are absent, no valid SNVs are found, or the FASTA cannot be read.

#### Minimal usage example

```python
from pymutation import PyMutation           # hypothetical package
mut = PyMutation.from_maf("my_cohort.maf")  # load variants
ctx, enriched = mut.trinucleotideMatrix("GRCh38.fa")
print(ctx.head())

```
