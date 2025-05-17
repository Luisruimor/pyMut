# pyMut

A Python library for visualizing genetic mutations from TSV files.

## Description

pyMut is a visualization tool for genetic mutation data. It allows generating summary visualizations to better understand patterns and distributions of mutations in the data.

## Features

- Statistical summary visualizations:
  - Variant Classification: Distribution of variant classifications
  - Variant Type: Distribution of variant types (SNP, INS, DEL, etc.)
  - SNV Class: Distribution of SNV classes (nucleotide changes like A>G, C>T, etc.)

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

# Generate a summary plot
fig = pyMut.summary_plot()
fig.savefig("summary.png")
```

## Future Development

In the future, pyMut will be expanded with more visualizations and functionalities for a more comprehensive analysis of genetic mutations, while maintaining simplicity and ease of use.
