### extract_signatures

#### Short description

Decomposes a **96 × samples** mutation-context matrix into mutational **signature profiles (W)** and their **contributions per sample (H)** using repeated non-negative matrix factorisation (NMF) for stability.

#### Signature

```python
def extract_signatures(
    contexts_df: pd.DataFrame,
    k: int,
    nrun: int = 30,
    pseudocount: float = 1e-4,
    random_seed: Optional[int] = None
) -> Dict:
```

#### Parameters

| Parameter     | Type           | Required | Description                                                                                       |
|---------------|----------------|----------|---------------------------------------------------------------------------------------------------|
| `contexts_df` | `pd.DataFrame` | **Yes**  | 96 × N matrix of raw counts produced by `trinucleotideMatrix`.                                    |
| `k`           | `int`          | **Yes**  | Number of signatures to extract (usually the `optimal_k` from `estimateSignatures`).              |
| `nrun`        | `int`          | No       | Independent NMF runs to perform; the best run (lowest KL divergence) is returned. *Default = 30*. |
| `pseudocount` | `float`        | No       | Small value added to avoid zeros before normalising to frequencies. *Default = 1e-4*.             |
| `random_seed` | \`int          | None\`   | Base seed for reproducibility; if `None`, each run uses a fresh random seed.                      |

#### Return value

`dict` with keys:

| Key                    | Type          | Meaning                                                          |
| ---------------------- | ------------- | ---------------------------------------------------------------- |
| `W`                    | `np.ndarray`  | **96 × k** signature matrix; each column sums to 1.              |
| `H`                    | `np.ndarray`  | **k × samples** exposure matrix scaled back to raw counts.       |
| `reconstruction_error` | `float`       | KL-divergence between input and reconstruction for the best run. |
| `best_run`             | `int`         | Index of the run that achieved the lowest error.                 |
| `all_errors`           | `list[float]` | Error from every run (length = `nrun`).                          |

None of these entries is ever `None`; the function raises if every run fails.

#### Exceptions

* `ImportError` – *scikit-learn* missing.
* `ValueError` – invalid arguments (e.g. `k` > samples), wrong input shape, or all NMF runs fail.

#### Minimal usage example

```python
# ctx is the matrix from trinucleotideMatrix
sig = extract_signatures(ctx, k=4, nrun=50, random_seed=123)

print("Chosen run:", sig["best_run"])
print("W shape:", sig["W"].shape)
```