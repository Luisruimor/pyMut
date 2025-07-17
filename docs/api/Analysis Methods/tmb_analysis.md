# TMB Analysis - Tumor Mutation Burden Analysis

The **calculate_tmb_analysis** method allows calculating the Tumor Mutation Burden (TMB) for each sample in a PyMutation object.

## What is TMB Analysis?

TMB is a metric that quantifies the total number of mutations in a tumor sample, normalized by the size of the interrogated genome. It is an important biomarker in oncology for predicting response to immunotherapy.

## Main Features

- **Automatic calculation per sample**: Analyzes each sample individually
- **Automatic column detection**: Automatically identifies variant classification columns
- **Non-synonymous mutations**: Distinguishes between total and non-synonymous mutations
- **Genome size normalization**: Calculates TMB per million bases
- **Global statistics**: Generates descriptive statistics of the dataset
- **Automatic export**: Saves results to TSV files

## Basic Usage

```python
from pyMut.io import read_maf

# Load data
py_mut = read_maf("mutations.maf")

# Basic TMB analysis
tmb_results = py_mut.calculate_tmb_analysis()

# Access results
analysis_df = tmb_results['analysis']
statistics_df = tmb_results['statistics']

print(f"Samples analyzed: {len(analysis_df)}")
print(f"Average TMB: {analysis_df['TMB_Total_Normalized'].mean():.2f} mut/Mb")
```

## Parameters

### variant_classification_column (str, optional)
- **Description**: Name of the column with variant classification
- **Automatic detection**: If None, automatically detects columns like 'Variant_Classification'
- **Example**: `"Variant_Classification"` or `"gencode_19_variant_classification"`

### genome_size_bp (int, default=60456963)
- **Description**: Size of the interrogated region in base pairs
- **WES (default)**: 60,456,963 bp (complete exome)
- **WGS**: ~3,000,000,000 bp (complete genome)
- **Panels**: Specific size of the panel used

### output_dir (str, default=".")
- **Description**: Directory where to save result files
- **Generated files**: `TMB_analysis.tsv` and `TMB_statistics.tsv`

### save_files (bool, default=True)
- **Description**: Whether to save results to TSV files
- **False**: Only returns DataFrames without saving files

## Return Value

Returns a dictionary with two DataFrames:

```python
{
    'analysis': DataFrame,      # Analysis per sample
    'statistics': DataFrame     # Global statistics
}
```

### 'analysis' DataFrame
```
Sample | Total_Mutations | Non_Synonymous_Mutations | TMB_Total_Normalized | TMB_Non_Synonymous_Normalized
SAMPLE_001 | 45 | 32 | 0.744 | 0.529
SAMPLE_002 | 123 | 89 | 2.034 | 1.472
```

### 'statistics' DataFrame
```
Metric | Count | Mean | Median | Min | Max | Q1 | Q3 | Std
Total_Mutations | 100 | 67.5 | 45.0 | 12 | 234 | 32.0 | 89.0 | 45.2
TMB_Total_Normalized | 100 | 1.117 | 0.744 | 0.198 | 3.871 | 0.529 | 1.472 | 0.748
```

## Non-Synonymous Mutation Types

The analysis automatically identifies mutations with biological impact:

```python
# Types considered non-synonymous
non_synonymous_types = {
    'MISSENSE_MUTATION',        # Amino acid change
    'NONSENSE_MUTATION',        # Premature stop codon
    'FRAME_SHIFT_DEL',          # Frame-changing deletion
    'FRAME_SHIFT_INS',          # Frame-changing insertion
    'NONSTOP_MUTATION',         # Loss of stop codon
    'TRANSLATION_START_SITE',   # Affects start site
    'SPLICE_SITE',              # Affects splicing site
    'IN_FRAME_DEL',             # In-frame deletion
    'IN_FRAME_INS',             # In-frame insertion
    'START_CODON_SNP',          # SNP in start codon
    'START_CODON_DEL',          # Deletion in start codon
    'START_CODON_INS',          # Insertion in start codon
    'STOP_CODON_DEL',           # Deletion in stop codon
    'STOP_CODON_INS'            # Insertion in stop codon
}
```

## Result Files

The method automatically generates two TSV files with analysis results:

### TMB_analysis.tsv - Analysis per Sample
Contains detailed results for each sample with columns:
- **Sample**: Sample identifier
- **Total_Mutations**: Total number of mutations
- **Non_Synonymous_Mutations**: Number of non-synonymous mutations
- **TMB_Total_Normalized**: TMB normalized by genome size (mut/Mb)
- **TMB_Non_Synonymous_Normalized**: Non-synonymous TMB normalized (mut/Mb)

### TMB_statistics.tsv - Global Statistics
Contains descriptive statistics for the entire dataset:
- **Count**: Number of samples
- **Mean**: Average value
- **Median**: Median value
- **Min/Max**: Minimum and maximum values
- **Q1/Q3**: First and third quartiles
- **Std**: Standard deviation

## Complete Example

```python
from pyMut.io import read_maf
import os

# Load TCGA data
py_mut = read_maf("src/pyMut/data/examples/tcga_laml.maf.gz")
print(f"Total mutations: {len(py_mut.data)}")
print(f"Samples: {len(py_mut.samples)}")

# Create output directory
os.makedirs("tmb_results", exist_ok=True)

# Complete TMB analysis
tmb_results = py_mut.calculate_tmb_analysis(
    variant_classification_column="Variant_Classification",
    genome_size_bp=60456963,  # WES standard
    output_dir="tmb_results",
    save_files=True
)

# Access results
analysis_df = tmb_results['analysis']
statistics_df = tmb_results['statistics']

# Display summary
print(f"\n=== TMB Analysis Summary ===")
print(f"Samples analyzed: {len(analysis_df)}")
print(f"Average TMB: {analysis_df['TMB_Total_Normalized'].mean():.3f} mut/Mb")
print(f"Median TMB: {analysis_df['TMB_Total_Normalized'].median():.3f} mut/Mb")
print(f"TMB range: {analysis_df['TMB_Total_Normalized'].min():.3f} - {analysis_df['TMB_Total_Normalized'].max():.3f} mut/Mb")

# Identify high TMB samples
high_tmb_threshold = 10  # mut/Mb
high_tmb_samples = analysis_df[analysis_df['TMB_Total_Normalized'] > high_tmb_threshold]
print(f"\nSamples with high TMB (>{high_tmb_threshold} mut/Mb): {len(high_tmb_samples)}")

if len(high_tmb_samples) > 0:
    print("Top 5 samples with highest TMB:")
    top_samples = high_tmb_samples.nlargest(5, 'TMB_Total_Normalized')
    for _, row in top_samples.iterrows():
        print(f"  - {row['Sample']}: {row['TMB_Total_Normalized']:.3f} mut/Mb")

# Non-synonymous vs total mutations
print(f"\n=== Mutation Type Analysis ===")
total_avg = analysis_df['Total_Mutations'].mean()
nonsyn_avg = analysis_df['Non_Synonymous_Mutations'].mean()
nonsyn_ratio = nonsyn_avg / total_avg if total_avg > 0 else 0

print(f"Average total mutations per sample: {total_avg:.1f}")
print(f"Average non-synonymous mutations per sample: {nonsyn_avg:.1f}")
print(f"Non-synonymous ratio: {nonsyn_ratio:.2%}")
```

## Advanced Usage

### Custom Genome Size

```python
# For whole genome sequencing
wgs_results = py_mut.calculate_tmb_analysis(
    genome_size_bp=3000000000,  # ~3 Gb for WGS
    output_dir="tmb_wgs_results"
)

# For targeted panel (example: 1.2 Mb panel)
panel_results = py_mut.calculate_tmb_analysis(
    genome_size_bp=1200000,  # 1.2 Mb panel
    output_dir="tmb_panel_results"
)
```

### Custom Variant Classification Column

```python
# If using different column name
custom_results = py_mut.calculate_tmb_analysis(
    variant_classification_column="gencode_19_variant_classification",
    output_dir="tmb_custom_results"
)
```

### Memory-efficient Analysis (no file saving)

```python
# For large datasets, avoid file I/O
tmb_results = py_mut.calculate_tmb_analysis(
    save_files=False  # Only return DataFrames
)

# Process results in memory
analysis_df = tmb_results['analysis']
# ... perform analysis ...
```

## TMB Interpretation

### Clinical Thresholds
```python
# Common TMB thresholds in clinical practice
def classify_tmb(tmb_value):
    if tmb_value >= 20:
        return "Very High"
    elif tmb_value >= 10:
        return "High"
    elif tmb_value >= 6:
        return "Intermediate"
    else:
        return "Low"

# Apply classification
analysis_df['TMB_Category'] = analysis_df['TMB_Total_Normalized'].apply(classify_tmb)
category_counts = analysis_df['TMB_Category'].value_counts()
print("TMB Distribution:")
for category, count in category_counts.items():
    print(f"  {category}: {count} samples")
```

### Statistical Analysis
```python
import numpy as np
from scipy import stats

# TMB distribution analysis
tmb_values = analysis_df['TMB_Total_Normalized']

# Test for normality
shapiro_stat, shapiro_p = stats.shapiro(tmb_values)
print(f"Shapiro-Wilk test p-value: {shapiro_p:.4f}")

# Log transformation if needed
if shapiro_p < 0.05:
    log_tmb = np.log1p(tmb_values)  # log(1+x) to handle zeros
    print("TMB distribution is not normal, consider log transformation")

# Percentile analysis
percentiles = [10, 25, 50, 75, 90, 95, 99]
tmb_percentiles = np.percentile(tmb_values, percentiles)
print("\nTMB Percentiles:")
for p, value in zip(percentiles, tmb_percentiles):
    print(f"  {p}th percentile: {value:.3f} mut/Mb")
```


## Performance Considerations

- **Memory usage**: TMB analysis is memory-efficient and processes samples iteratively
- **Genome size**: Ensure correct genome size for accurate normalization
- **Column detection**: Automatic detection works for standard MAF formats

## Common Use Cases

1. **Immunotherapy prediction**: Identify high TMB samples likely to respond to checkpoint inhibitors
2. **Cohort characterization**: Describe TMB distribution in a patient cohort
3. **Quality control**: Identify samples with unusual mutation patterns
4. **Comparative analysis**: Compare TMB between different cancer types or treatments
5. **Biomarker development**: Use TMB as a continuous or categorical biomarker

## Technical Notes

- TMB calculation follows standard clinical guidelines
- Non-synonymous mutations are identified using standard variant classifications
- Results are automatically saved with timestamps for reproducibility
- The method handles missing or malformed variant classifications gracefully
- Compatible with both MAF and VCF-derived data formats