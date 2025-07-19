# pyMut 🧬

[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://img.shields.io/badge/tests-17%20passed-green.svg)](tests/)

A professional Python library for visualizing genetic mutations from mutation data files, inspired by tools like **Maftools** and **Mutscape**.

## 🎯 Features

pyMut provides comprehensive mutation visualization capabilities:

### 📊 **Complete Summary Visualizations**
- **Variant Classification**: Distribution of mutation types (Missense, Nonsense, etc.)
- **Variant Type**: Distribution of variant types (SNP, INS, DEL, etc.)
- **SNV Class**: Distribution of nucleotide changes (A>G, C>T, etc.)
- **Variants per Sample (TMB)**: Tumor mutation burden analysis
- **Variant Classification Summary**: Boxplot analysis across samples
- **Top Mutated Genes**: Most frequently mutated genes with two analysis modes

### 🎨 **Professional Visualization Features**
- **High-quality graphics** with publication-ready output (DPI 300+)
- **Consistent color schemes** across all visualizations
- **Interactive mode** for data exploration
- **Automatic format detection** (wide vs long format)
- **Flexible customization** options

### 🔧 **Advanced Capabilities**
- **Automatic data preprocessing** from FUNCOTATION fields
- **Multi-format support** (pipe-separated, slash-separated genotypes)
- **Sample detection** (TCGA format and custom identifiers)
- **Memory-efficient** processing of large datasets
- **Comprehensive error handling** and validation

## 🚀 Quick Start

### Installation

```bash
pip install pyMut
```

### Basic Usage

```python
from pyMut import PyMutation
import pandas as pd

# Load your mutation data
data = pd.read_csv("mutations.tsv", sep="\t")

# Create PyMutation object
py_mut = PyMutation(data)

# Configure high-quality output (recommended)
PyMutation.configure_high_quality_plots()

# Generate complete summary plot
summary_fig = py_mut.summary_plot(
    title="Mutation Analysis Summary",
    figsize=(16, 12),
    max_samples=200,
    top_genes_count=10
)
summary_fig.savefig("mutation_summary.png")  # Automatically high quality!
```

### Individual Visualizations

```python
# Tumor Mutation Burden (TMB) analysis
tmb_fig = py_mut.variants_per_sample_plot(
    title="Tumor Mutation Burden per Sample",
    max_samples=100
)
tmb_fig.savefig("tmb_analysis.png")

# Top mutated genes (by variant count)
genes_fig = py_mut.top_mutated_genes_plot(
    mode="variants",  # Count total variants
    count=15,
    title="Top 15 Most Mutated Genes"
)
genes_fig.savefig("top_genes_variants.png")

# Top mutated genes (by sample prevalence)
prevalence_fig = py_mut.top_mutated_genes_plot(
    mode="samples",   # Count affected samples percentage
    count=15,
    title="Top 15 Genes by Sample Prevalence"
)
prevalence_fig.savefig("top_genes_prevalence.png")

# Variant classification boxplot
boxplot_fig = py_mut.variant_classification_summary_plot(
    title="Variant Classification Distribution Across Samples"
)
boxplot_fig.savefig("variant_boxplot.png")
```

## 📁 Supported Data Formats

### Long Format (Recommended)
Each row represents one mutation:
```
Hugo_Symbol | Variant_Classification | Tumor_Sample_Barcode | REF | ALT
GENE1      | Missense_Mutation      | SAMPLE_001          | A   | G
GENE2      | Nonsense_Mutation      | SAMPLE_001          | C   | T
```

### Wide Format
Samples as columns with genotype information:
```
Hugo_Symbol | Variant_Classification | SAMPLE_001 | SAMPLE_002 | SAMPLE_003
GENE1      | Missense_Mutation      | A|G        | A|A        | A|G
GENE2      | Nonsense_Mutation      | C|T        | C|C        | C|C
```

### Automatic Detection
pyMut automatically detects your data format and handles:
- **TCGA sample identifiers** (TCGA-XX-YYYY format)
- **Custom sample naming** conventions
- **Multiple genotype formats**: `A|G`, `A/G`, or custom notation
- **FUNCOTATION field parsing** for variant extraction

## 🎨 High-Quality Output

### Automatic High-Quality Configuration
```python
# Configure once for all figures
PyMutation.configure_high_quality_plots()

# All subsequent saves will be high quality automatically
fig.savefig("my_plot.png")  # DPI 300, optimized margins, PNG compression
```

### Manual Quality Control
```python
# Save with custom quality settings
py_mut.save_figure(fig, "publication_plot.png", dpi=600)

# Multiple formats
fig.savefig("plot.pdf", dpi=300, bbox_inches='tight')  # PDF for publications
fig.savefig("plot.svg", bbox_inches='tight')           # SVG for editing
```

## 🔬 Advanced Features

### Custom Parameters
```python
# Highly customizable visualizations
py_mut.variants_per_sample_plot(
    figsize=(14, 8),
    max_samples=50,
    variant_column="Custom_Variant_Col",
    sample_column="Custom_Sample_Col"
)
```

### Data Preprocessing
```python
# Automatic extraction from FUNCOTATION fields
# No manual preprocessing required - pyMut handles it automatically
```

## 📋 Requirements

- **Python**: 3.7+ (tested on 3.7, 3.8, 3.9, 3.10, 3.11)
- **Core dependencies**:
  - `pandas` >= 1.2.0
  - `matplotlib` >= 3.3.0
  - `numpy` >= 1.19.0
- **Optional dependencies**:
  - `seaborn` >= 0.11.0 (enhanced styling)

## 🧪 Testing

Run the comprehensive test suite:

```bash
# Clean output (recommended)
./run_clean_tests.sh

# Standard pytest
pytest

# Detailed output
python -m pytest tests/ -v
```

**Test Coverage**: 17 tests covering all major functionality, validation, and edge cases.

## 📚 Documentation

- **[Complete Documentation](docs/)** - Comprehensive guides and API reference
- **[Installation Guide](docs/user-guide/installation.md)** - Detailed installation instructions
- **[User Guide](docs/user-guide/basic-usage.md)** - Step-by-step tutorials
- **[API Reference](docs/api/)** - Complete API documentation
- **[Examples](docs/examples/)** - Real-world usage examples

## 🤝 Contributing

Contributions are welcome! Please:

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feature/amazing-feature`)
3. **Run tests** to ensure everything works (`./run_clean_tests.sh`)
4. **Commit** your changes (`git commit -m 'Add amazing feature'`)
5. **Push** to the branch (`git push origin feature/amazing-feature`)
6. **Open** a Pull Request

### Development Setup
```bash
git clone https://github.com/your-username/pyMut.git
cd pyMut
pip install -e .
./run_clean_tests.sh  # Verify everything works
```

## 📄 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

## 🎯 Comparison with Other Tools

| Feature | pyMut | Maftools (R) | Mutscape (R) |
|---------|-------|--------------|--------------|
| Language | Python 🐍 | R | R |
| Installation | `pip install` | Complex R setup | Complex R setup |
| Data Format | Auto-detection | Manual preparation | Manual preparation |
| Interactive | ✅ | ❌ | ❌ |
| High-Quality Output | Auto-configuration | Manual setup | Manual setup |
| Testing | 17 comprehensive tests | Limited | Limited |
| Documentation | Complete | Partial | Limited |

## 📈 Roadmap

- [ ] **Pathway analysis** integration
- [ ] **Survival analysis** plots
- [ ] **Mutation signatures** analysis
- [ ] **Copy number** visualization
- [ ] **Multi-sample** comparison tools
- [ ] **Export to** common formats (VCF, MAF)

## 🙏 Acknowledgments

Inspired by the excellent work of:
- **Maftools** (R package for mutation analysis)
- **Mutscape** (R package for mutation landscape)
- **TCGA** consortium for standardized data formats

---

**⭐ If pyMut helps your research, please consider giving it a star!**

*Made with ❤️ for the genomics community*