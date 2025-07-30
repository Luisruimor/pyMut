# Tissue Expression Filter

The **filter_by_tissue_expression** method allows filtering mutation data based on gene expression in specific tissues using cancer RNA consensus data.

## What is the Tissue Expression Filter?

This filter allows identifying and filtering mutations in genes that are sufficiently expressed (or not expressed) in specific tissues/cancer types. It uses cancer RNA consensus data to determine if a gene is "active" in a particular tissue according to a defined expression threshold.

## Main Features

- **Multi-tissue filtering**: Allows specifying multiple tissues with independent thresholds
- **Automatic column detection**: Automatically identifies gene symbol columns
- **Smart caching**: Expression data is loaded once and kept in cache
- **Bidirectional filtering**: Can filter expressed or non-expressed genes
- **TCGA codes**: Uses standard TCGA codes for cancer types
- **Robust validation**: Complete validation of input parameters

## Basic Usage

```python
from pyMut.io import read_maf

# Load data
py_mut = read_maf("mutations.maf")

# Filter genes expressed in bladder cancer (threshold 5)
filtered_mut = py_mut.filter_by_tissue_expression([('BLCA', 5)])

# Filter genes NOT expressed in lung adenocarcinoma
not_expressed_mut = py_mut.filter_by_tissue_expression(
    [('LUAD', 4)], 
    keep_expressed=False
)

# Filter genes expressed in multiple tissues with different thresholds
multi_tissue_mut = py_mut.filter_by_tissue_expression([
    ('BLCA', 5),    # Bladder cancer, threshold 5
    ('BRCA', 3),    # Breast cancer, threshold 3
    ('LUAD', 4),    # Lung adenocarcinoma, threshold 4
    ('COAD', 6)     # Colon adenocarcinoma, threshold 6
])

print(f"Original mutations: {len(py_mut.data)}")
print(f"Filtered mutations: {len(filtered_mut.data)}")
```

## Parameters

### tissues (List[Tuple[str, float]], required)
- **Description**: List of tuples with tissue and threshold specifications
- **Format**: `[('tissue_code', threshold), ...]`
- **TCGA Codes**: Uses standard codes like 'BLCA', 'BRCA', 'LUAD', etc.
- **Example**: `[('BLCA', 5), ('BRCA', 3)]`

### keep_expressed (bool, default=True)
- **Description**: Determines which mutations to keep
- **True**: Keeps genes expressed in at least one of the specified tissues
- **False**: Keeps genes NOT expressed in any of the specified tissues

## Return Value

Returns a new **PyMutation** object with data filtered according to the specified tissue expression criteria.

```python
# The returned object maintains the same structure
filtered_mut = py_mut.filter_by_tissue_expression([('BLCA', 5)])
print(type(filtered_mut))  # <class 'pyMut.core.PyMutation'>
```

## Supported TCGA Tissue Codes

The filter uses standard TCGA codes for different cancer types:

```python
# Examples of common TCGA codes
tissue_codes = {
    'BLCA': 'Bladder Urothelial Carcinoma',
    'BRCA': 'Breast Invasive Carcinoma', 
    'COAD': 'Colon Adenocarcinoma',
    'LUAD': 'Lung Adenocarcinoma',
    'LUSC': 'Lung Squamous Cell Carcinoma',
    'PRAD': 'Prostate Adenocarcinoma',
    'THCA': 'Thyroid Carcinoma',
    'LIHC': 'Liver Hepatocellular Carcinoma',
    'STAD': 'Stomach Adenocarcinoma',
    'SKCM': 'Skin Cutaneous Melanoma'
    # ... and many more
}
```

## Helper Function: tissue_expression

In addition to the filtering method, the `tissue_expression` helper function is available for individual checks:

```python
from pyMut.filters.tissue_expression import tissue_expression
import pandas as pd

# Check expression using gene symbol directly
is_expressed = tissue_expression("TSPAN6", ["BLCA", 5])
print(f"TSPAN6 expressed in BLCA (threshold 5): {is_expressed}")

# Check expression using a data row
row = pd.Series({'Hugo_Symbol': 'TSPAN6', 'Chromosome': 'X'})
is_expressed = tissue_expression(row, ["BRCA", 10])
print(f"TSPAN6 expressed in BRCA (threshold 10): {is_expressed}")
```

### tissue_expression Parameters

- **data**: `Union[str, pd.Series]` - Gene symbol or data row
- **tissue**: `List[Union[str, float]]` - `[tissue_code, threshold]`

## Advanced Examples

### Conditional Filtering by Cancer Type

```python
# Filter mutations relevant to breast cancer
breast_cancer_mut = py_mut.filter_by_tissue_expression([
    ('BRCA', 5)  # Genes expressed in breast cancer
])

# Filter silenced genes in lung cancer
lung_silenced_mut = py_mut.filter_by_tissue_expression([
    ('LUAD', 3),
    ('LUSC', 3)
], keep_expressed=False)
```

### Multi-Tissue Analysis

```python
# Genes expressed in at least one of multiple digestive cancers
digestive_cancers_mut = py_mut.filter_by_tissue_expression([
    ('COAD', 4),    # Colon
    ('STAD', 4),    # Stomach  
    ('LIHC', 5),    # Liver
    ('PAAD', 6)     # Pancreas
])

print(f"Mutations in genes expressed in digestive cancers: {len(digestive_cancers_mut.data)}")
```

### Combination with Other Filters

```python
# Combine tissue expression filter with other filters
filtered_mut = (py_mut
    .filter_by_tissue_expression([('BRCA', 5)])
    .filter_by_pass()
    .filter_by_chromosome(['1', '2', '3']))

print(f"Mutations after combined filters: {len(filtered_mut.data)}")
```

## Error Handling

The method includes robust validation and error handling:

```python
try:
    # Error: empty list
    py_mut.filter_by_tissue_expression([])
except ValueError as e:
    print(f"Error: {e}")

try:
    # Error: incorrect tuple format
    py_mut.filter_by_tissue_expression([('BLCA', 5, 'extra')])
except ValueError as e:
    print(f"Error: {e}")

try:
    # Error: tissue code is not string
    py_mut.filter_by_tissue_expression([(123, 5)])
except ValueError as e:
    print(f"Error: {e}")
```

## Performance Considerations

- **Caching**: Expression data is loaded once and kept in cache
- **Memory**: Filtering creates a copy of the data, does not modify the original object
- **Scalability**: Efficient for large datasets thanks to smart caching

## Common Use Cases

1. **Cancer-specific analysis**: Filter mutations relevant to a specific cancer type
2. **Silenced genes**: Identify mutations in non-expressed genes (potential tumor suppressors)
3. **Comparative analysis**: Compare mutation patterns between different cancer types
4. **Variant prioritization**: Focus analysis on genes active in relevant tissues

## Technical Notes

- Expression data comes from `rna_cancer_consensus.json`
- Automatically detects columns like 'Hugo_Symbol', 'Gene_Symbol', etc.
- Filtering preserves all metadata from the original PyMutation object
- Compatible with MAF format data and converted VCF