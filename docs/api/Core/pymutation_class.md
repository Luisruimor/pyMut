# PyMutation – Main Class for Mutation Analysis

The **`PyMutation`** class is the central object in **pyMut**, encapsulating mutation data and providing methods for both analysis and visualisation.

## What is PyMutation?

`PyMutation` represents an entire mutation dataset: it stores the data in a structured table, keeps metadata about its origin, and exposes analysis and plotting helpers.

## Object Structure

### Main Attributes

```python
class PyMutation:
    def __init__(self, data: pd.DataFrame, metadata: MutationMetadata, samples: List[str]):
        self.data = data           # DataFrame of mutations (VCF-like columns)
        self.samples = samples     # List of sample IDs
        self.metadata = metadata   # Provenance and configuration metadata
```

### `data` Attribute (`pd.DataFrame`)

VCF-style mutation table with standard columns:

```
CHROM | POS | ID | REF | ALT | QUAL | FILTER | SAMPLE_001 | SAMPLE_002 | ANNOTATIONS | ...
chr1  | 100 | .  | A   | G   | .    | .      | A|G        | A|A        | ...         | ...
chr2  | 200 | .  | C   | T   | .    | .      | C|C        | C|T        | ...         | ...
```

### `samples` Attribute (`List[str]`)

Sample identifiers:

```python
['SAMPLE_001', 'SAMPLE_002', 'SAMPLE_003', ...]
```

### `metadata` Attribute (`MutationMetadata`)

Provenance and settings:

```python
metadata.source_format   # "MAF" or "VCF"
metadata.file_path       # Path to the original file
metadata.loaded_at       # Timestamp when loaded
metadata.filters         # Applied filters
metadata.fasta           # Reference FASTA
metadata.notes           # Original file comments
```

## Creating `PyMutation` Objects

### From a MAF file

```python
from pyMut.io import read_maf

# Recommended entry-point
py_mut = read_maf("mutations.maf")
```

### Manual creation (advanced)

```python
import pandas as pd
from pyMut.core import PyMutation, MutationMetadata

# Build a DataFrame with the required columns
data = pd.DataFrame({
    'CHROM': ['chr1', 'chr2'],
    'POS': [100, 200],
    'REF': ['A', 'C'],
    'ALT': ['G', 'T'],
    'SAMPLE_001': ['A|G', 'C|C'],
    'SAMPLE_002': ['A|A', 'C|T'],
    'Hugo_Symbol': ['GENE1', 'GENE2'],
    'Variant_Classification': ['Missense_Mutation', 'Nonsense_Mutation']
})

# Create metadata
metadata = MutationMetadata(
    source_format="Manual",
    file_path="manual_creation",
    filters=["."],
    fasta="",
    notes="Manually created"
)

# Instantiate the object
samples = ['SAMPLE_001', 'SAMPLE_002']
py_mut = PyMutation(data, metadata, samples)
```

## Visualisation Methods

### Summary Plot – Complete Analysis

```python
# Produce a 6-panel overview figure
fig = py_mut.summary_plot(
    title="My Mutation Analysis",
    figsize=(16, 12),
    max_samples=100,
    top_genes_count=15,
)
fig.savefig("summary_analysis.png")
```

### Individual Plots

#### Variant Classification

```python
fig = py_mut.variant_classification_plot(
    title="Mutation Type Distribution",
    figsize=(10, 6)
)
```

#### Variant Types

```python
fig = py_mut.variant_type_plot(
    title="Variant Type Distribution",
    figsize=(10, 6)
)
```

#### SNV Classes

```python
fig = py_mut.snv_class_plot(
    title="Nucleotide Change Distribution",
    figsize=(10, 6)
)
```

#### Variants per Sample (TMB)

```python
fig = py_mut.variants_per_sample_plot(
    title="Tumour Mutation Burden per Sample",
    max_samples=50,
    figsize=(12, 6)
)
```

#### Per-Sample Classification Summary

```python
fig = py_mut.variant_classification_summary_plot(
    title="Per-Sample Classification Summary",
    figsize=(10, 6)
)
```

#### Most-Mutated Genes

```python
fig = py_mut.top_mutated_genes_plot(
    title="Top Mutated Genes",
    count=20,
    figsize=(10, 8)
)
```

## Analysis Methods

### Tumour Mutation Burden (TMB)

```python
# Full TMB analysis
tmb_results = py_mut.calculate_tmb_analysis(
    variant_classification_column="Variant_Classification",
    genome_size_bp=60_456_963,  # WES
    output_dir="results",
    save_files=True
)

# Access results
analysis_df = tmb_results['analysis']
statistics_df = tmb_results['statistics']
```

## Utility Methods

### High-Quality Plot Configuration

```python
# Global matplotlib settings for publication-ready output
PyMutation.configure_high_quality_plots()

# All subsequent figures are rendered at high DPI
fig = py_mut.summary_plot()
fig.savefig("high_quality_plot.png")  # Automatically DPI=300
```

### Centralised Figure Saving

```python
fig = py_mut.summary_plot()

# Basic save
py_mut.save_figure(fig, "analysis.png")

# Custom save
py_mut.save_figure(fig, "analysis.pdf", dpi=600, bbox_inches="tight")
```

## End-to-End Example

```python
from pyMut.io import read_maf
from pyMut import PyMutation
import matplotlib.pyplot as plt

# 1. LOAD DATA
print("📂 Loading MAF data...")
py_mut = read_maf("src/pyMut/data/examples/tcga_laml.maf.gz")

# 2. EXPLORE STRUCTURE
print("✅ Data successfully loaded:")
print(f"   • Samples: {len(py_mut.samples)}")
print(f"   • Mutations: {len(py_mut.data)}")
print(f"   • Source format: {py_mut.metadata.source_format}")
print(f"   • File: {py_mut.metadata.file_path}")

# 3. HIGH-QUALITY SETTINGS
PyMutation.configure_high_quality_plots()

# 4. COMPREHENSIVE VISUAL ANALYSIS
print("\n📊 Generating visual analysis...")
summary_fig = py_mut.summary_plot(
    title="TCGA-LAML: Comprehensive Mutation Analysis",
    figsize=(18, 14),
    max_samples=150,
    top_genes_count=20
)

# Save in multiple formats
py_mut.save_figure(summary_fig, "tcga_laml_summary.png")
py_mut.save_figure(summary_fig, "tcga_laml_summary.pdf", dpi=300)

# 5. TMB ANALYSIS
print("\n🧬 Calculating TMB...")
tmb_results = py_mut.calculate_tmb_analysis(
    genome_size_bp=60_456_963,  # Standard WES
    output_dir="results/tmb_analysis",
    save_files=True
)

# Explore TMB output
analysis_df = tmb_results['analysis']
statistics_df = tmb_results['statistics']

print(f"   • Analysed samples: {len(analysis_df)}")
print(f"   • Mean TMB: {analysis_df['TMB_Total_Normalized'].mean():.3f} mut/Mb")
print(f"   • Median TMB: {analysis_df['TMB_Total_Normalized'].median():.3f} mut/Mb")

# 6. INDIVIDUAL VISUALISATIONS
print("\n📈 Generating specific plots...")

# TMB per sample
tmb_fig = py_mut.variants_per_sample_plot(
    title="Tumour Mutation Burden (TMB)",
    max_samples=100
)
py_mut.save_figure(tmb_fig, "tmb_per_sample.png")

# Top mutated genes
genes_fig = py_mut.top_mutated_genes_plot(
    title="Top 25 Mutated Genes",
    count=25
)
py_mut.save_figure(genes_fig, "top_mutated_genes.png")

# Mutation types
classification_fig = py_mut.variant_classification_plot(
    title="Mutation Type Distribution"
)
py_mut.save_figure(classification_fig, "variant_classification.png")

# 7. CUSTOM ANALYSIS
print("\n🔍 Custom analysis...")

# Identify high-TMB samples
high_tmb_threshold = 10  # mut/Mb
high_tmb_samples = analysis_df[
    analysis_df['TMB_Total_Normalized'] > high_tmb_threshold
]

print(f"   • High-TMB samples (>{high_tmb_threshold} mut/Mb): {len(high_tmb_samples)}")

if len(high_tmb_samples) > 0:
    print("   • Top 5 samples by TMB:")
    top_samples = high_tmb_samples.nlargest(5, 'TMB_Total_Normalized')
    for _, row in top_samples.iterrows():
        print(f"     - {row['Sample']}: {row['TMB_Total_Normalized']:.3f} mut/Mb")

# Most frequently mutated genes
gene_counts = py_mut.data['Hugo_Symbol'].value_counts().head(10)
print("\n   • Top 10 mutated genes:")
for gene, count in gene_counts.items():
    print(f"     - {gene}: {count} mutations")

print("\n✅ Complete analysis finished!")
print("📁 Generated files:")
print("   • tcga_laml_summary.png/pdf – full overview")
print("   • tmb_per_sample.png – TMB per sample")
print("   • top_mutated_genes.png – most mutated genes")
print("   • variant_classification.png – mutation types")
print("   • results/tmb_analysis/ – detailed TMB results")
```

## Data Access

### Basic Exploration

```python
# Basic info
print(f"Shape: {py_mut.data.shape}")
print(f"Columns: {list(py_mut.data.columns)}")
print(f"Samples: {py_mut.samples}")

# Preview rows
print(py_mut.data.head())

# Column info
print(py_mut.data.info())
```

### Filtering

```python
# By gene
tp53_data = py_mut.data[py_mut.data['Hugo_Symbol'] == 'TP53']

# By mutation type
missense_data = py_mut.data[
    py_mut.data['Variant_Classification'] == 'Missense_Mutation'
]

# By chromosome
chr1_data = py_mut.data[py_mut.data['CHROM'] == 'chr1']
```

### Descriptive Statistics

```python
# Count by mutation type
mutation_counts = py_mut.data['Variant_Classification'].value_counts()
print(mutation_counts)

# Unique genes
unique_genes = py_mut.data['Hugo_Symbol'].nunique()
print(f"Unique genes: {unique_genes}")

# Mutations per sample
mutations_per_sample = {}
for sample in py_mut.samples:
    count = 0
    for _, row in py_mut.data.iterrows():
        genotype = row[sample]
        ref = row['REF']
        if genotype != f"{ref}|{ref}":  # not REF|REF
            count += 1
    mutations_per_sample[sample] = count

print("Mutations per sample:")
for sample, count in sorted(mutations_per_sample.items(),
                            key=lambda x: x[1], reverse=True)[:10]:
    print(f"  {sample}: {count}")
```

## Integration with Other Analyses

### Export Data for External Tools

```python
# Full export
py_mut.data.to_csv("mutations_export.tsv", sep="\t", index=False)

# Export mutations for a single sample
sample = py_mut.samples[0]
sample_mut = py_mut.data[
    py_mut.data[sample] != f"{py_mut.data['REF']}|{py_mut.data['REF']}"
]
sample_mut.to_csv(f"{sample}_mutations.tsv", sep="\t", index=False)

# Export metadata
metadata_info = {
    "source_format": py_mut.metadata.source_format,
    "file_path": py_mut.metadata.file_path,
    "loaded_at": str(py_mut.metadata.loaded_at),
    "total_samples": len(py_mut.samples),
    "total_mutations": len(py_mut.data)
}

import json
with open("metadata.json", "w") as f:
    json.dump(metadata_info, f, indent=2)
```

### Combine with `pandas` for Advanced Analyses

```python
import pandas as pd
import numpy as np

def analyse_gene_cooccurrence(py_mut, genes):
    """Analyse mutation co-occurrence in specific genes."""
    rows = []
    for sample in py_mut.samples:
        mutated = []
        for gene in genes:
            gene_rows = py_mut.data[py_mut.data['Hugo_Symbol'] == gene]
            if gene_rows.empty:
                continue
            for _, row in gene_rows.iterrows():
                genotype = row[sample]
                ref = row['REF']
                if genotype != f"{ref}|{ref}":
                    mutated.append(gene)
                    break
        rows.append({
            "Sample": sample,
            "Mutated_Genes": mutated,
            "Gene_Count": len(mutated)
        })
    return pd.DataFrame(rows)

# Example usage
oncogenes = ["TP53", "KRAS", "PIK3CA", "APC", "EGFR"]
cooccurrence_df = analyse_gene_cooccurrence(py_mut, oncogenes)
print(cooccurrence_df.head())
```

## Best Practices

### 1. Recommended Initial Setup

```python
from pyMut.io import read_maf
from pyMut import PyMutation

# Always enable high-quality plotting early
PyMutation.configure_high_quality_plots()

# Load data
py_mut = read_maf("data.maf")
```

### 2. Data Validation

```python
assert len(py_mut.samples) > 0, "No samples found"
assert len(py_mut.data) > 0, "No mutations found"
assert "Hugo_Symbol" in py_mut.data.columns, "Missing gene column"
```

### 3. Systematic Workflow

```python
# 1. Visual overview
summary_fig = py_mut.summary_plot()

# 2. TMB analysis
tmb_results = py_mut.calculate_tmb_analysis()

# 3. Specific plots as needed
# 4. Custom data analysis
```

### 4. File Management

```python
import os

# Create results directory
os.makedirs("results", exist_ok=True)

# Save all figures centrally
py_mut.save_figure(summary_fig, "results/summary.png")
```