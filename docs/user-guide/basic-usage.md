# Basic Usage

## Import the library

To start using pyMut, first import the main class, `PyMutation`:

```python
from pyMut import PyMutation
import pandas as pd
```

## Load data

pyMut works with mutation data in pandas DataFrame format. You can load your data from a TSV file:

```python
# Load data from a TSV file
data = pd.read_csv("path/to/mutations.tsv", sep="\t")

# Create a PyMutation object
pyMut = PyMutation(data)
```

## Generate visualizations

### Summary Plot

A summary plot shows general mutation statistics. This visualization includes:

- Variant Classification: Distribution of variant classifications
- Variant Type: Distribution of variant types (SNP, INS, DEL, etc.)
- SNV Class: Distribution of SNV classes (nucleotide changes like A>G, C>T, etc.)

To generate a summary plot:

```python
# Generate a summary plot
fig = pyMut.summary_plot(title="Mutation Summary")

# Save the plot as a PNG file
fig.savefig("summary.png")
```

### Customization

You can customize the summary plot by changing the parameters:

```python
# Generate a customized summary plot
fig = pyMut.summary_plot(
    figsize=(15, 12),     # Figure size
    title="Mutation Analysis in Cancer Patients"  # Custom title
)

# Save the plot as a high-resolution PNG file
fig.savefig("custom_summary.png", dpi=300)
```

### Interactive Visualization

You can also display the plot interactively:

```python
# Generate and display the plot interactively
fig = pyMut.summary_plot()
import matplotlib.pyplot as plt
plt.show()
``` 