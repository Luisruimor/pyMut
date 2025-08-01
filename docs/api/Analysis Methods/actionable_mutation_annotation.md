### actionable_mutations_oncokb

#### Short description

Annotates mutation data using the OncoKB API to add actionable mutation information and therapeutic implications.

#### Signature

```python
def actionable_mutations_oncokb(self, token: str, batch_size: int = 5000, timeout: int = 30, max_retries: int = 3, retry_backoff: float = 1.0) -> pd.DataFrame:
```

#### Parameters

| Parameter        | Type    | Required | Description                                                                                    |
| ---------------- | ------- | -------- | ---------------------------------------------------------------------------------------------- |
| `token`          | `str`   | Yes      | OncoKB API authentication token required for accessing the OncoKB service.                    |
| `batch_size`     | `int`   | No       | Maximum number of variants to send in a single API request (default: 5000).                  |
| `timeout`        | `int`   | No       | Timeout for API requests in seconds (default: 30).                                            |
| `max_retries`    | `int`   | No       | Maximum number of retries for failed API requests (default: 3).                               |
| `retry_backoff`  | `float` | No       | Backoff factor for retries, controls delay between retry attempts (default: 1.0).             |

#### Return value

Returns a `pd.DataFrame` containing the original self.data DataFrame with OncoKB annotations added as columns. The annotations include therapeutic implications and actionable mutation information.

#### Exceptions

List only those the user should handle:

* `ValueError`: if the DataFrame doesn't contain the necessary data for export (missing required columns like CHROM, POS, REF, ALT) or if the reference genome is invalid.
* `requests.exceptions.RequestException`: if there's an error with the API request that can't be resolved with retries.

#### Minimal usage example

```python
>>> # Assuming you have a valid OncoKB API token
>>> token = "your_oncokb_api_token"
>>> annotated_data = pymutation_obj.actionable_mutations_oncokb(token)
>>> print(annotated_data.columns)
```