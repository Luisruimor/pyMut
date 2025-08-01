# pyMut 🧬

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/pymut-bio.svg)](https://badge.fury.io/py/pymut-bio)

A Python library designed for the preprocessing, analysis and visualization of somatic genetic variants in standard formats such as VCF and MAF.

## 🚀 Quick Start

### Installation

#### Option 1: Basic Installation (pip)

```bash
pip install pymut-bio
```

**Note**: The pip installation provides core functionality for mutation data visualization, but some advanced features may be limited as certain bioinformatics tools are not available through PyPI.

#### Option 2: Full Installation (Recommended - Conda)

For complete functionality including all bioinformatics tools, use the conda environment:

```bash
# 0) Descargar el environment.yml (elige curl o wget)
curl -fsSL https://raw.githubusercontent.com/Luisruimor/pyMut/main/environment.yml -o environment.yml
# ó:
# wget -O environment.yml https://raw.githubusercontent.com/Luisruimor/pyMut/main/environment.yml

# 1) crear el entorno (añade tus binarios al environment.yml)
conda env create -f environment.yml

# 2) activar el entorno
conda activate NOMBRE-DEL-ENTORNO

# 3) instalar tu librería desde PyPI en ese entorno
pip install pymut-bio
```

The conda environment includes essential bioinformatics tools:
- **bcftools**: VCF/BCF file manipulation
- **ensembl-vep**: Variant Effect Predictor
- **htslib**: High-throughput sequencing data processing
- **tabix**: Generic indexer for TAB-delimited genome position files

These tools enable advanced genomic data processing capabilities that are not available with pip-only installation.

## 📚 Documentation

- **[Complete Documentation](https://luisruimor.github.io/pyMut/)** - Comprehensive guides and API reference
- **[Installation Guide](https://luisruimor.github.io/pyMut/installation/)** - Detailed installation instructions
- **[API Reference](https://luisruimor.github.io/pyMut/api/Core/pymutation_class/)** - Complete API documentation
- **[Examples](https://luisruimor.github.io/pyMut/examples/data/input_read_maf/#example-loading-tcga-laml-maf-file)** - Real-world usage examples


## 📋 Requirements

| Librería                  | Dependencias inmediatas                                                                                                                                                                                                                                                                     |
|---------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **duckdb** 1.3.2          | – Ninguna                                                                                                                                                                                                                                                                                   |
| **fastparquet** 2024.11.0 | – cramjam ≥ 2.3<br>– fsspec<br>– numpy<br>– packaging<br>– pandas ≥ 1.5.0                                                                                                                                                                                                                   |
| **matplotlib** 3.10.3     | – contourpy ≥ 1.0.1<br>– cycler ≥ 0.10<br>– fonttools ≥ 4.22.0<br>– kiwisolver ≥ 1.3.1<br>– numpy ≥ 1.23<br>– packaging ≥ 20.0<br>– pillow ≥ 8<br>– pyparsing ≥ 2.3.1<br>– python-dateutil ≥ 2.7                                                                                            |
| **mkdocs** 1.6.1          | – click ≥ 7.0<br>– colorama ≥ 0.4<br>– ghp-import ≥ 1.0<br>– jinja2 ≥ 2.11.1<br>– markdown ≥ 3.3.6<br>– markupsafe ≥ 2.0.1<br>– mergedeep ≥ 1.3.4<br>– mkdocs-get-deps ≥ 0.2.0<br>– packaging ≥ 20.5<br>– pathspec ≥ 0.11.1<br>– pyyaml ≥ 5.1<br>– pyyaml-env-tag ≥ 0.1<br>– watchdog ≥ 2.0 |
| **numpy** 1.26.4          | – Ninguna                                                                                                                                                                                                                                                                                   |
| **pandas** 2.3.1          | – numpy ≥ 1.22.4<br>– python-dateutil ≥ 2.8.2<br>– pytz ≥ 2020.1<br>– tzdata ≥ 2022.7                                                                                                                                                                                                       |
| **pyarrow** 14.0.2        | – numpy ≥ 1.16.6                                                                                                                                                                                                                                                                            |
| **pyensembl** 2.3.13      | – datacache ≥ 1.4.0,<2.0.0<br>– gtfparse ≥ 2.5.0,<3.0.0<br>– memoized-property ≥ 1.0.2<br>– pylint ≥ 2.17.2,<3.0.0<br>– serializable ≥ 0.2.1,<1.0.0<br>– tinytimer ≥ 0.0.0,<1.0.0<br>– typechecks ≥ 0.0.2,<1.0.0                                                                            |
| **pyfaidx** 0.8.1.4       | – packaging                                                                                                                                                                                                                                                                                 |
| **requests** 2.32.4       | – certifi ≥ 2017.4.17<br>– charset-normalizer ≥ 2,<4<br>– idna ≥ 2.5,<4<br>– urllib3 ≥ 1.21.1,<3                                                                                                                                                                                            |
| **scikit-learn** 1.7.1    | – joblib ≥ 1.2.0<br>– numpy ≥ 1.22.0<br>– scipy ≥ 1.8.0<br>– threadpoolctl ≥ 3.1.0                                                                                                                                                                                                          |
| **scipy** 1.11. 4         | – numpy ≥ 1.21.6,<1.28.0                                                                                                                                                                                                                                                                    |
| **seaborn** 0.13.2        | – matplotlib ≥ 3.4,<3.6.1 or >3.6.1<br>– numpy ≥ 1.20,<1.24.0 or >1.24.0<br>– pandas ≥ 1.2                                                                                                                                                                                                  |
| **urllib3** 2.5.0         | – Ninguna                                                                                                                                                                                                                                                                                   |


## 📄 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

## 🎯 Comparison with Other Tools

| FUNCTIONAL CRITERIA                         | PYMUT (PROPOSAL)   | MUTSCAPE              | MAFTOOLS              |
|---------------------------------------------|--------------------|-----------------------|-----------------------|
| Input formats                               | VCF & MAF (native) | MAF                   | MAF                   |
| VEP annotation                              | ✓                  |                       |                       |
| Genomic range filtering                     | ✓                  | ✓                     | ✓                     |
| PASS category variant filtering             | ✓                  | ✓                     |                       |
| Sample filtering                            | ✓                  |                       | ✓                     |
| Tissue expression filtering                 | ✓                  | ✓                     |                       |
| File format transformation                  | ✓                  | ✓ *(VCF to MAF only)* | ✓ *(VCF to MAF only)* |
| File combination                            | ✓                  | ✓                     |                       |
| Significantly mutated genes (SMG) detection |                    | ✓                     |                       |
| Cancer-related gene annotation              | ✓                  | ✓                     |                       |
| Tumor mutational burden (TMB) calculation   | ✓                  | ✓                     |                       |
| Mutational signature identification         | ✓                  |                       |                       |
| Medical implications mutation annotation    | ✓                  | ✓                     |                       |
| PFAM annotation support                     | ✓                  |                       | ✓                     |
