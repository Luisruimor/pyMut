# filter_by_chrom_sample - Chromosome and Sample Filter

The **filter_by_chrom_sample** method allows filtering PyMutation data by chromosome and/or sample, providing granular control over which data to include in the analysis.

## What is filter_by_chrom_sample?

It is a versatile method that allows filtering data by chromosome, sample, or both criteria simultaneously. In addition to filtering rows, it also handles column filtering when samples are specified, maintaining data format integrity.

## Main Features

- **Dual filtering**: By chromosome and/or sample in a single operation
- **Row and column filtering**: Removes both irrelevant rows and columns
- **Multiple value support**: Accepts lists of chromosomes and samples
- **MAF/VCF compatibility**: Handles both formats automatically
- **Metadata preservation**: Records all applied filters
- **Automatic validation**: Verifies the existence of chromosomes and samples
- **Detailed logging**: Provides information about the filtering process

## Basic Usage

```python
from pyMut.input import read_maf

# Load data
py_mut = read_maf("mutations.maf")

# Filter by chromosome only
chr17_data = py_mut.filter_by_chrom_sample(chrom="chr17")

# Filter by sample only
sample_data = py_mut.filter_by_chrom_sample(sample="TCGA-AB-2802")

# Filter by both criteria
chr17_sample = py_mut.filter_by_chrom_sample(
    chrom="chr17", 
    sample="TCGA-AB-2802"
)

print(f"Mutations in chr17: {len(chr17_data.data)}")
print(f"Mutations in sample: {len(sample_data.data)}")
print(f"Mutations chr17 + sample: {len(chr17_sample.data)}")
```

## Parameters

### chrom (str, list, optional)
- **Description**: Chromosome(s) to filter
- **Accepted formats**: `"chr17"`, `"17"`, `["chr1", "chr17"]`, `["X", "Y"]`
- **Normalization**: Automatically normalized to standard format
- **Example**: `"chr17"` or `["chr1", "chr2", "chrX"]`

### sample (str, list, optional)
- **Description**: Sample(s) to filter
- **Format**: Sample identifiers as they appear in the data
- **Example**: `"TCGA-AB-2802"` or `["TCGA-AB-2802", "TCGA-AB-2803"]`

### sample_column (str, optional)
- **Description**: Name of the column containing sample information
- **Default**: `"Tumor_Sample_Barcode"` (MAF standard)
- **Usage**: For data with non-standard column names

## Filtering Behavior

### Chromosome-Only Filtering
```python
# Keeps all columns, filters only rows
chr_filtered = py_mut.filter_by_chrom_sample(chrom="chr17")

# Result: Only mutations in chr17, all samples preserved
print(f"Unique chromosomes: {chr_filtered.data['CHROM'].unique()}")
print(f"Preserved samples: {len(chr_filtered.samples)}")
```

### Sample-Only Filtering
```python
# Filters both rows and columns
sample_filtered = py_mut.filter_by_chrom_sample(sample="TCGA-AB-2802")

# Result: Only mutations from the sample, only relevant columns
print(f"Samples in data: {sample_filtered.samples}")
print(f"Sample columns: {[col for col in sample_filtered.data.columns if 'TCGA' in col]}")
```

### Combined Filtering
```python
# Applies both filters simultaneously
combined = py_mut.filter_by_chrom_sample(
    chrom=["chr17", "chrX"], 
    sample=["TCGA-AB-2802", "TCGA-AB-2803"]
)

print(f"Chromosomes: {combined.data['CHROM'].unique()}")
print(f"Samples: {combined.samples}")
```

## Detailed Examples

### Specific Chromosome Analysis

```python
from pyMut.input import read_maf
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)

# Load TCGA data
py_mut = read_maf("src/pyMut/data/examples/tcga_laml.maf.gz")
print(f"Original data: {len(py_mut.data)} mutations, {len(py_mut.samples)} samples")

# Analysis of oncologically interesting chromosomes
oncology_chromosomes = ["chr17", "chr13", "chr3", "chr7", "chr12"]

print(f"\n=== Chromosome Analysis ===")
for chrom in oncology_chromosomes:
    try:
        # Filter by chromosome
        chrom_data = py_mut.filter_by_chrom_sample(chrom=chrom)
        
        print(f"\n{chrom}:")
        print(f"  • Mutations: {len(chrom_data.data)}")
        
        if len(chrom_data.data) > 0:
            # Most mutated genes in this chromosome
            top_genes = chrom_data.data['Hugo_Symbol'].value_counts().head(3)
            print(f"  • Top genes:")
            for gene, count in top_genes.items():
                print(f"    - {gene}: {count} mutations")
            
            # Affected samples
            affected_samples = chrom_data.data['Tumor_Sample_Barcode'].nunique()
            print(f"  • Affected samples: {affected_samples}")
            
    except Exception as e:
        print(f"❌ Error processing {chrom}: {e}")
```

### Sample-Specific Analysis

```python
# Select samples of interest
samples_of_interest = py_mut.samples[:5]  # First 5 samples

print(f"\n=== Sample Analysis ===")
for sample in samples_of_interest:
    try:
        # Filter by sample
        sample_data = py_mut.filter_by_chrom_sample(sample=sample)
        
        print(f"\n{sample}:")
        print(f"  • Total mutations: {len(sample_data.data)}")
        
        if len(sample_data.data) > 0:
            # Mutation types
            mutation_types = sample_data.data['Variant_Classification'].value_counts()
            print(f"  • Main mutation types:")
            for mut_type, count in mutation_types.head(3).items():
                print(f"    - {mut_type}: {count}")
            
            # Chromosomal distribution
            chrom_dist = sample_data.data['CHROM'].value_counts()
            print(f"  • Most affected chromosomes:")
            for chrom, count in chrom_dist.head(3).items():
                print(f"    - {chrom}: {count} mutations")
                
    except Exception as e:
        print(f"❌ Error processing {sample}: {e}")
```

### Multi-Sample Comparative Analysis

```python
# Compare multiple samples
sample_groups = {
    "Group_A": ["TCGA-AB-2802", "TCGA-AB-2803", "TCGA-AB-2804"],
    "Group_B": ["TCGA-AB-2805", "TCGA-AB-2806", "TCGA-AB-2807"]
}

print(f"\n=== Multi-Sample Comparative Analysis ===")
for group_name, samples in sample_groups.items():
    try:
        # Filter by sample group
        group_data = py_mut.filter_by_chrom_sample(sample=samples)
        
        print(f"\n{group_name} ({len(samples)} samples):")
        print(f"  • Total mutations: {len(group_data.data)}")
        print(f"  • Average mutations per sample: {len(group_data.data)/len(samples):.1f}")
        
        if len(group_data.data) > 0:
            # Most mutated genes in the group
            top_genes = group_data.data['Hugo_Symbol'].value_counts().head(5)
            print(f"  • Most mutated genes:")
            for gene, count in top_genes.items():
                print(f"    - {gene}: {count} mutations")
                
    except Exception as e:
        print(f"❌ Error processing {group_name}: {e}")
```

### Combined Filtering Analysis

```python
# Analysis combining chromosome and sample filters
print(f"\n=== Combined Analysis: chr17 + Specific Samples ===")

# Select samples with high mutation load
high_mutation_samples = []
for sample in py_mut.samples[:10]:  # Check first 10 samples
    sample_data = py_mut.filter_by_chrom_sample(sample=sample)
    if len(sample_data.data) > 50:  # Samples with >50 mutations
        high_mutation_samples.append(sample)

print(f"Samples with high mutation load: {len(high_mutation_samples)}")

if high_mutation_samples:
    # Filter chr17 in high-mutation samples
    chr17_high_mut = py_mut.filter_by_chrom_sample(
        chrom="chr17",
        sample=high_mutation_samples
    )
    
    print(f"chr17 mutations in high-mutation samples: {len(chr17_high_mut.data)}")
    
    if len(chr17_high_mut.data) > 0:
        # Genes most affected in chr17
        chr17_genes = chr17_high_mut.data['Hugo_Symbol'].value_counts()
        print(f"Most mutated genes in chr17:")
        for gene, count in chr17_genes.head(5).items():
            print(f"  - {gene}: {count} mutations")
```

## Advanced Usage

### Filtering by Chromosome Lists

```python
# Filter multiple chromosomes simultaneously
sex_chromosomes = py_mut.filter_by_chrom_sample(chrom=["chrX", "chrY"])
autosomes_1_5 = py_mut.filter_by_chrom_sample(chrom=["chr1", "chr2", "chr3", "chr4", "chr5"])

print(f"Mutations in sex chromosomes: {len(sex_chromosomes.data)}")
print(f"Mutations in chromosomes 1-5: {len(autosomes_1_5.data)}")
```

### Filtering by Sample List

```python
# Define a list of samples of interest
samples_of_interest = [
    "TCGA-AB-2802",
    "TCGA-AB-2803", 
    "TCGA-AB-2804",
    "TCGA-AB-2805",
    "TCGA-AB-2806"
]

# Filter data to include only specified samples
filtered_data = py_mut.filter_by_chrom_sample(sample=samples_of_interest)
print(f"Mutations in selected samples: {len(filtered_data.data)}")
print(f"Number of samples in filtered data: {len(filtered_data.samples)}")

# You can also combine with chromosome filtering
chr17_selected_samples = py_mut.filter_by_chrom_sample(
    chrom="chr17",
    sample=samples_of_interest
)
print(f"chr17 mutations in selected samples: {len(chr17_selected_samples.data)}")

# Filter by sample subsets for comparative analysis
group_a = ["TCGA-AB-2802", "TCGA-AB-2803"]
group_b = ["TCGA-AB-2804", "TCGA-AB-2805"]

group_a_data = py_mut.filter_by_chrom_sample(sample=group_a)
group_b_data = py_mut.filter_by_chrom_sample(sample=group_b)

print(f"Group A mutations: {len(group_a_data.data)}")
print(f"Group B mutations: {len(group_b_data.data)}")
```

## Error Handling and Validation

```python
# The method includes robust validation
try:
    # Valid filtering
    valid_filter = py_mut.filter_by_chrom_sample(chrom="chr17", sample="TCGA-AB-2802")
    print("✅ Valid filtering successful")
    
except ValueError as e:
    print(f"❌ Validation error: {e}")

# Handle non-existent chromosomes
try:
    invalid_chrom = py_mut.filter_by_chrom_sample(chrom="chr99")
except KeyError as e:
    print(f"❌ Chromosome not found: {e}")

# Handle non-existent samples
try:
    invalid_sample = py_mut.filter_by_chrom_sample(sample="NON_EXISTENT_SAMPLE")
except KeyError as e:
    print(f"❌ Sample not found: {e}")
```

## Integration with Other Methods

### Chaining with Other Filters

```python
# Chain multiple filters
filtered_data = (py_mut
    .filter_by_chrom_sample(chrom="chr17")
    .filter_by_pass()  # If available
    .filter_by_tissue_expression([('BRCA', 5)]))  # If available

print(f"Data after chained filters: {len(filtered_data.data)}")
```

### Combination with Analysis Methods

```python
# Filter and then analyze
chr17_data = py_mut.filter_by_chrom_sample(chrom="chr17")

# Perform TMB analysis on filtered data
if hasattr(chr17_data, 'calculate_tmb_analysis'):
    tmb_results = chr17_data.calculate_tmb_analysis()
    print(f"TMB analysis on chr17: {len(tmb_results['analysis'])} samples")

# Generate visualizations
if hasattr(chr17_data, 'summary_plot'):
    fig = chr17_data.summary_plot(title="chr17 Mutation Analysis")
```

## Metadata and Tracking

```python
# Filters are automatically recorded in metadata
original = py_mut
filtered = py_mut.filter_by_chrom_sample(chrom="chr17", sample="TCGA-AB-2802")

print("Applied filters:")
for filter_info in filtered.metadata.filters:
    print(f"  - {filter_info}")

# Example output:
# - filter_by_chrom_sample:chrom=chr17,sample=TCGA-AB-2802
```

## Common Use Cases

1. **Chromosome-specific analysis**: Focus on specific chromosomes of oncological interest
2. **Sample subset analysis**: Analyze specific patient cohorts
3. **Quality control**: Remove problematic samples or chromosomes
4. **Comparative studies**: Compare mutation patterns between sample groups
5. **Performance optimization**: Reduce dataset size for faster processing

## Technical Notes

- The method preserves all original metadata and sample information
- Chromosome normalization handles formats with and without "chr" prefix
- Sample filtering also removes corresponding columns from the data
- Compatible with both MAF and VCF-derived data formats
- Thread-safe and suitable for parallel processing
- Maintains data integrity and format consistency