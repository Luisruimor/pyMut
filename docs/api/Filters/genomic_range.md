# Genomic Range Filters - region and gen_region

The **region** and **gen_region** methods allow filtering PyMutation data by specific genomic location, either by chromosomal coordinates or by gene name.

## What are Genomic Range Filters?

These are methods that allow extracting mutation subsets based on their location in the genome, facilitating analysis of specific regions of interest.

## Main Features

- **Coordinate filtering**: Specify chromosome, start and end positions
- **Gene filtering**: Automatically search for coordinates of a specific gene
- **pyarrow optimization**: Uses pyarrow for fast queries when available
- **Metadata preservation**: Maintains original information and records applied filters
- **Automatic validation**: Verifies that coordinates are valid
- **Detailed logging**: Provides information about the filtering process

## region Method - Coordinate-based Filtering

### Basic Usage

```python
from pyMut.input import read_maf

# Load data
py_mut = read_maf("mutations.maf")

# Filter by specific region
# Chromosome 17, positions 7,500,000 to 7,600,000
tp53_region = py_mut.region("chr17", 7500000, 7600000)

print(f"Mutations in region: {len(tp53_region.data)}")
```

### region Parameters

#### chrom (str) [required]
- **Description**: Chromosome to filter
- **Accepted formats**: `"chr17"`, `"17"`, `"X"`, `"Y"`, `"chrX"`, `"chrY"`
- **Example**: `"chr17"`

#### start (int) [required]
- **Description**: Start position of the region (inclusive)
- **Coordinates**: 1-based (genomic standard)
- **Example**: `7500000`

#### end (int) [required]
- **Description**: End position of the region (inclusive)
- **Coordinates**: 1-based (genomic standard)
- **Example**: `7600000`

### region Examples

```python
# TP53 gene region (chromosome 17)
tp53_mutations = py_mut.region("chr17", 7571720, 7590868)

# BRCA1 gene region (chromosome 17)
brca1_mutations = py_mut.region("17", 43044295, 43125483)

# Complete chromosome X (large region example)
chrx_mutations = py_mut.region("X", 1, 156040895)

# Specific region on chromosome Y
chry_region = py_mut.region("chrY", 2781479, 2781479)  # Specific position
```

## gen_region Method - Gene Name Filtering

### Basic Usage

```python
# Filter by gene name (automatic)
tp53_mutations = py_mut.gen_region("TP53")
brca1_mutations = py_mut.gen_region("BRCA1")

print(f"Mutations in TP53: {len(tp53_mutations.data)}")
print(f"Mutations in BRCA1: {len(brca1_mutations.data)}")
```

### gen_region Parameters

#### gen_name (str) [required]
- **Description**: Name of the gene to search
- **Format**: Official gene symbol (HUGO)
- **Examples**: `"TP53"`, `"BRCA1"`, `"EGFR"`, `"KRAS"`

### Gene Database

The method uses an internal database with genomic information:

```python
# Supported genes include:
available_genes = [
    "TP53",      # Chromosome 17: 7,571,720-7,590,868
    "BRCA1",     # Chromosome 17: 43,044,295-43,125,483
    "BRCA2",     # Chromosome 13: 32,315,086-32,400,268
    "EGFR",      # Chromosome 7: 55,019,017-55,211,628
    "KRAS",      # Chromosome 12: 25,205,246-25,250,929
    "PIK3CA",    # Chromosome 3: 179,148,114-179,240,093
    # ... and many more
]
```

## Complete Example

```python
from pyMut.input import read_maf
import logging

# Configure logging to see details
logging.basicConfig(level=logging.INFO)

# Load TCGA data
py_mut = read_maf("src/pyMut/data/examples/tcga_laml.maf.gz")
print(f"Total mutations: {len(py_mut.data)}")

# Analysis of specific genes
genes_of_interest = ["TP53", "KRAS", "PIK3CA", "EGFR"]

for gene in genes_of_interest:
    try:
        # Filter by gene
        gene_mutations = py_mut.gen_region(gene)
        print(f"\n=== Analysis of {gene} ===")
        print(f"Mutations found: {len(gene_mutations.data)}")
        
        if len(gene_mutations.data) > 0:
            # Most common mutation types
            types = gene_mutations.data['Variant_Classification'].value_counts()
            print(f"Mutation types:")
            for mutation_type, count in types.head(3).items():
                print(f"  - {mutation_type}: {count}")
            
            # Affected samples
            samples = gene_mutations.data['Tumor_Sample_Barcode'].nunique()
            print(f"Affected samples: {samples}")
            
    except Exception as e:
        print(f"❌ Error processing {gene}: {e}")

# Specific region analysis
print(f"\n=== Chromosomal Region Analysis ===")
# Cancer gene-rich region on chromosome 17
chr17_region = py_mut.region("chr17", 7000000, 8000000)
print(f"Mutations in chr17:7M-8M: {len(chr17_region.data)}")

# Genes in this region
if len(chr17_region.data) > 0:
    genes_region = chr17_region.data['Hugo_Symbol'].value_counts()
    print("Most mutated genes in the region:")
    for gene, count in genes_region.head(5).items():
        print(f"  - {gene}: {count} mutations")
```

## Coordinate Validation

```python
# Automatic validations performed:
try:
    # Valid coordinates
    valid_region = py_mut.region("chr1", 1000000, 2000000)
    
    # Error: start > end
    invalid_region = py_mut.region("chr1", 2000000, 1000000)
    
except ValueError as e:
    print(f"❌ Invalid coordinates: {e}")

# Error: non-existent chromosome
try:
    invalid_chr = py_mut.region("chr99", 1000000, 2000000)
except KeyError as e:
    print(f"❌ Chromosome not found: {e}")
```

## Metadata Handling

```python
# Filters are automatically registered
original = py_mut
filtered = py_mut.region("chr17", 7500000, 7600000)

print("Applied filters:")
for filter_info in filtered.metadata.filters:
    print(f"  - {filter_info}")

# Example output:
# - region:chr17:7500000-7600000
```

## Common Use Cases

### Candidate Gene Analysis

```python
# List of oncological genes of interest
cancer_genes = ["TP53", "KRAS", "PIK3CA", "EGFR", "BRAF", "APC"]

# Analyze each gene individually
gene_analysis = {}
for gene in cancer_genes:
    mutations = py_mut.gen_region(gene)
    gene_analysis[gene] = {
        'mutations': len(mutations.data),
        'samples': mutations.data['Tumor_Sample_Barcode'].nunique() if len(mutations.data) > 0 else 0
    }

# Summary
for gene, stats in gene_analysis.items():
    print(f"{gene}: {stats['mutations']} mutations in {stats['samples']} samples")
```

### Chromosomal Region Analysis

```python
# Chromosomal arm analysis
# Chromosome 17p (short arm)
chr17p = py_mut.region("chr17", 1, 22300000)

# Chromosome 17q (long arm)  
chr17q = py_mut.region("chr17", 22300001, 83257441)

print(f"Mutations in 17p: {len(chr17p.data)}")
print(f"Mutations in 17q: {len(chr17q.data)}")
```

### Hotspot Analysis

```python
# Known hotspot regions
hotspots = {
    "TP53_DBD": ("chr17", 7571720, 7590868),      # DNA binding domain
    "KRAS_G12": ("chr12", 25245274, 25245276),    # Codon 12
    "PIK3CA_E545": ("chr3", 179218303, 179218305) # Codon 545
}

for name, (chrom, start, end) in hotspots.items():
    hotspot_muts = py_mut.region(chrom, start, end)
    print(f"{name}: {len(hotspot_muts.data)} mutations")
```

## Combination with Other Filters

```python
# Combine genomic filters with sample filters
# 1. Filter by region
tp53_region = py_mut.region("chr17", 7571720, 7590868)

# 2. Filter by specific samples
specific_samples = ["TCGA-AB-2802", "TCGA-AB-2803"]
tp53_samples = tp53_region.filter_by_chrom_sample(sample=specific_samples)

print(f"TP53 mutations in specific samples: {len(tp53_samples.data)}")
```