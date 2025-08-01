### compare_signatures

#### Short description

Matches **extracted signatures** against the **COSMIC catalogue** via cosine similarity, reporting the best COSMIC hit for each and optionally the full similarity matrix.

#### Signature

```python
def compare_signatures(
    W: np.ndarray,
    cosmic_path: str,
    min_cosine: float = 0.6,
    return_matrix: bool = False
) -> Dict:
```

#### Parameters

| Parameter       | Type         | Required | Description                                                                          |
| --------------- | ------------ | -------- | ------------------------------------------------------------------------------------ |
| `W`             | `np.ndarray` | **Yes**  | 96 × k matrix of normalised signatures (columns sum to 1) from `extract_signatures`. |
| `cosmic_path`   | `str`        | **Yes**  | File path to COSMIC SBS catalogue (tab-separated, first column = context labels).    |
| `min_cosine`    | `float`      | No       | Threshold below which a COSMIC match is reported as “No match”. *Default = 0.6*.     |
| `return_matrix` | `bool`       | No       | If `True`, include the complete k × N cosine-similarity matrix in the output.        |

#### Return value

`dict` with keys:

| Key             | Type                | Meaning                                                                              |
| --------------- | ------------------- | ------------------------------------------------------------------------------------ |
| `summary_df`    | `pd.DataFrame`      | One row per column in `W` with `Signature_W`, `Best_COSMIC`, `Cosine`, `Aetiology`.  |
| `cosine_matrix` | `np.ndarray` *opt.* | k × N matrix of pairwise cosine similarities (only present if `return_matrix=True`). |

`summary_df` is never empty; if no COSMIC signature reaches `min_cosine`, `Best_COSMIC` is set to *“No match”*.

#### Exceptions

* `FileNotFoundError` – the COSMIC file cannot be located.
* `ValueError` – mismatched shapes/contexts, empty catalog, or malformed data.
* `ImportError` – *scikit-learn* missing.

#### Minimal usage example

```python
# W is the signature matrix from extract_signatures
cmp = compare_signatures(W, "data/COSMIC_v3.4_SBS_GRCh38.txt", min_cosine=0.7)

print(cmp["summary_df"].head())
```
