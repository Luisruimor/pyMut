# pyMut

A Python library for visualizing genetic mutations from TSV files.

## Description

pyMut is a visualization tool for genetic mutation data, inspired by tools like Maftools and Mutscape. It allows generating summary visualizations to understand genetic mutations in a TSV-format dataset.

## Features

- Mutation summary visualizations:
  - **Variant Classification**: Distribution of variant classifications
  - **Variant Type**: Distribution of variant types (SNP, INS, DEL, etc.)
  - **SNV Class**: Distribution of SNV classes (nucleotide changes like A>G, C>T, etc.)
  - **Variants per Sample**: Distribution of variants per sample and median (TMB)
  - **Variant Classification Summary**: Box and whisker plot showing the distribution of each variant type across samples, allowing identification of which classifications show greater variability between patients
  - **Top Mutated Genes**: Horizontal bar chart showing the most mutated genes with two modes:
    - "variants" mode: Shows the total number of variants per gene
    - "samples" mode: Shows the percentage of affected samples per gene

## Installation

```bash
pip install pyMut
```

## Basic Usage

```python
from pyMut import PyMutation
import pandas as pd

# Load mutation data from a TSV file
data = pd.read_csv("mutations.tsv", sep="\t")

# Create a PyMutation object
pyMut = PyMutation(data)

# Generate a complete summary plot
fig = pyMut.summary_plot(title="Mutation Summary")
fig.savefig("summary.png")

# Generate individual visualizations
tmb_fig = pyMut.variants_per_sample_plot(title="Tumor Mutation Burden per Sample")
tmb_fig.savefig("variants_per_sample.png")

boxplot_fig = pyMut.variant_classification_summary_plot(title="Variant Distribution per Sample")
boxplot_fig.savefig("variant_classification_summary.png")

# Top mutated genes visualization (by number of variants)
top_genes_fig = pyMut.top_mutated_genes_plot(
    title="Top 10 Most Mutated Genes",
    mode="variants",  # Count total number of variants
    count=10  # Show top 10 genes
)
top_genes_fig.savefig("top_mutated_genes_variants.png")

# Top mutated genes visualization (by sample prevalence)
top_genes_samples_fig = pyMut.top_mutated_genes_plot(
    title="Top 10 Most Prevalent Genes",
    mode="samples",  # Count percentage of affected samples
    count=10  # Show top 10 genes
)
top_genes_samples_fig.savefig("top_mutated_genes_samples.png")
```

## Supported Data Formats

pyMut can work with data in two main formats:

- **Long Format**: Each row represents a mutation, with columns like `Variant_Classification` and `Tumor_Sample_Barcode`.
- **Wide Format**: Samples are represented as columns (e.g., columns with names like `TCGA-XX-YYYY`).

The library automatically detects the format and adapts visualizations accordingly. For wide format, pyMut can interpret data in different notations:
- Pipe-separated genotypes (`|`): such as "A|B" where A and B are different alleles
- Slash-separated genotypes (`/`): such as "A/B"
- Other formats: numerical or textual values indicating variant presence

## Requirements

- Python 3.7+
- pandas
- matplotlib
- numpy
- seaborn

## Documentation

For more information, see the [complete documentation](https://pymut.readthedocs.io/).

## Contributing

Contributions are welcome. Please open an issue first to discuss what you would like to change.

## License

[MIT](LICENSE)