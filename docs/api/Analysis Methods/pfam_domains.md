### pfam_domains

#### Short description

Summarise already-annotated Pfam protein-domain information in a `PyMutation` object, returning the most frequent domains or hotspots.&#x20;

#### Signature

```python
def pfam_domains(
    *,
    aa_column: str = "aa_pos",
    summarize_by: str = "PfamDomain",
    top_n: int = 10,
    include_synonymous: bool = False
) -> pandas.DataFrame:
```

#### Parameters

| Parameter            | Type   | Required | Description                                                                      |
| -------------------- | ------ | -------- | -------------------------------------------------------------------------------- |
| `aa_column`          | `str`  | No       | Column with amino-acid positions. Default `"aa_pos"`.                            |
| `summarize_by`       | `str`  | No       | Either `"PfamDomain"` (aggregate by domain) or `"AAPos"` (per-residue hotspots). |
| `top_n`              | `int`  | No       | Number of rows to keep in the output (sorted by variant count). Default `10`.    |
| `include_synonymous` | `bool` | No       | If `False` (default) silent variants are excluded from the summary.              |

#### Return value

`pandas.DataFrame` containing, for each reported line, counts of **variants** and **unique genes**:

* When `summarize_by="PfamDomain"` → columns `pfam_id`, `pfam_name`, `n_genes`, `n_variants`.
* When `summarize_by="AAPos"` → columns depend on grouping (`uniprot`, `aa_pos`, Pfam columns) plus `n_variants`, `n_genes`.

#### Exceptions

* `PfamAnnotationError` – Pfam columns (`pfam_id`, `pfam_name`) missing, or `Hugo_Symbol` absent when required.
* `ValueError` – `summarize_by` is not `"PfamDomain"` nor `"AAPos"`.

#### Minimal usage example

```python
# Assuming you have already run annotate_pfam()
summary = pymut.pfam_domains(
    summarize_by="PfamDomain",
    top_n=20
)

print(summary.head())
```