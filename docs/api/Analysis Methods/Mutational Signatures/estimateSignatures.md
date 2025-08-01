### estimateSignatures

#### Short description

Determines the **optimal number of mutational signatures** in a 96 × samples matrix by running multiple non-negative matrix factorization (NMF) decompositions and evaluating stability metrics.

#### Signature

```python
def estimateSignatures(
    contexts_df: pd.DataFrame,
    nMin: int = 2,
    nTry: int = 6,
    nrun: int = 5,
    parallel: int = 4,
    pConstant: Optional[float] = None
) -> Dict:
```

#### Parameters

| Parameter     | Type           | Required | Description                                                       |                                                                                        |
| ------------- | -------------- | -------- | ----------------------------------------------------------------- | -------------------------------------------------------------------------------------- |
| `contexts_df` | `pd.DataFrame` | **Yes**  | 96 × samples count matrix produced by `trinucleotideMatrix`.      |                                                                                        |
| `nMin`        | `int`          | No       | Smallest number of signatures (`k`) to test. *Default = 2*.       |                                                                                        |
| `nTry`        | `int`          | No       | Largest `k` to test (inclusive). *Default = 6*.                   |                                                                                        |
| `nrun`        | `int`          | No       | Independent NMF runs per `k` to assess robustness. *Default = 5*. |                                                                                        |
| `parallel`    | `int`          | No       | CPU threads used for parallel NMF fits. *Default = 4*.            |                                                                                        |
| `pConstant`   | \`float        | None\`   | No                                                                | Small positive value added if the matrix is extremely sparse; leave `None` to disable. |

#### Return value

`dict` with the keys:

| Key                 | Type           | Meaning                                                                                        |
| ------------------- | -------------- | ---------------------------------------------------------------------------------------------- |
| `metrics`           | `pd.DataFrame` | Stability statistics for each tested `k` (mean RSS, dispersion, cophenetic correlation, etc.). |
| `models`            | `list`         | All successful NMF model results (`W`, `H`, RSS, run index…).                                  |
| `optimal_k`         | `int`          | Suggested best number of signatures based on cophenetic drop-off.                              |
| `normalized_matrix` | `np.ndarray`   | Input matrix after column-wise frequency normalisation (used for NMF).                         |
| `original_matrix`   | `np.ndarray`   | Raw count matrix (same values as `contexts_df.values`).                                        |

None of these items are ever `None`; if every NMF fit fails, the function raises instead of returning.

#### Exceptions

* `ImportError` – *scikit-learn* or *scipy* not installed.
* `ValueError` – invalid inputs (wrong shape, impossible `nMin/nTry`, all decompositions fail).

#### Minimal usage example

```python
# ctx is the 96 × samples matrix from trinucleotideMatrix
results = estimateSignatures(ctx, nMin=2, nTry=8, nrun=10, parallel=6)

print("Optimal k =", results["optimal_k"])
print(results["metrics"][["k", "mean_rss", "cophenetic_corr"]])

```
