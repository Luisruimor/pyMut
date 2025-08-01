### combine_pymutations

#### Short description

Combines two **PyMutation** objects into a single instance, deduplicating variants, merging sample columns, and ensuring both input cohorts share the same genome assembly.

#### Signature

```python
def combine_pymutations(pymut1: PyMutation, pymut2: PyMutation) -> PyMutation:
```

#### Parameters

| Parameter | Type         | Required | Description                                                                            |
| --------- | ------------ | -------- | -------------------------------------------------------------------------------------- |
| `pymut1`  | `PyMutation` | Yes      | First mutation cohort. Its assembly must match that of `pymut2`.                       |
| `pymut2`  | `PyMutation` | Yes      | Second mutation cohort to merge with `pymut1`. Samples and annotations are integrated. |

#### Return value

`PyMutation` – **new** object containing the union of variants, samples, and annotations from both inputs. Originals are left untouched.

#### Exceptions

* `ValueError`: if `pymut1` and `pymut2` have different `assembly` values.

#### Minimal usage example

```python
from pymutation import PyMutation
from pymutation.analysis import combine_pymutations

cohort_a = PyMutation.from_maf("tumours_setA.maf")
cohort_b = PyMutation.from_maf("tumours_setB.maf")

merged = combine_pymutations(cohort_a, cohort_b)
print(f"Merged cohort: {len(merged.samples)} samples, {len(merged.data)} variants")
```
