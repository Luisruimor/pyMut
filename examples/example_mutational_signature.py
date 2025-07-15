#!/usr/bin/env python3
"""
Example script demonstrating mutational signature analysis in pyMut.

This script shows how to extract and visualize mutational signatures
from mutation data using Non-negative Matrix Factorization (NMF).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyMut import PyMutation
import os

# Configure matplotlib for high quality plots
PyMutation.configure_high_quality_plots()

def main():
    """Main function to demonstrate mutational signature analysis."""
    
    print("=== pyMut Mutational Signature Analysis Example ===\n")
    
    # Option 1: Load real mutation data
    # Uncomment the following lines if you have a MAF file:
    # data = pd.read_csv("mutations.maf", sep="\t")
    
    # Option 2: Generate synthetic mutation data for demonstration
    print("Generating synthetic mutation data for demonstration...")
    data = generate_synthetic_mutation_data(n_samples=50, n_mutations=5000)
    
    # Create PyMutation object with metadata
    from pyMut.core import MutationMetadata
    metadata = MutationMetadata(
        source_format="synthetic",
        file_path="synthetic_data.tsv",
        filters=[],
        fasta="synthetic.fasta"
    )
    py_mut = PyMutation(data, metadata)
    print(f"\nLoaded {len(data)} mutations from {data['Tumor_Sample_Barcode'].nunique()} samples")
    
    # Example 1: Basic mutational signature analysis
    print("\n1. Basic mutational signature analysis (3 signatures)...")
    fig1 = py_mut.mutational_signature_plot(
        n_signatures=3,
        title="Mutational Signature Analysis - 3 Signatures"
    )
    fig1.savefig("mutational_signatures_basic.png")
    print("   Saved: mutational_signatures_basic.png")
    
    # Example 2: Analysis with more signatures
    print("\n2. Analysis with 5 signatures...")
    fig2 = py_mut.mutational_signature_plot(
        n_signatures=5,
        figsize=(22, 28),
        title="Mutational Signature Analysis - 5 Signatures"
    )
    fig2.savefig("mutational_signatures_5sig.png")
    print("   Saved: mutational_signatures_5sig.png")
    
    # Example 3: With custom context column
    print("\n3. Analysis with specific context column...")
    # Add a context column if not present
    if 'Trinucleotide_Context' not in data.columns:
        data['Trinucleotide_Context'] = generate_trinucleotide_contexts(data)
        py_mut = PyMutation(data, metadata)  # Recreate with new data
    
    fig3 = py_mut.mutational_signature_plot(
        n_signatures=3,
        context_column='Trinucleotide_Context',
        title="Mutational Signatures with Trinucleotide Context"
    )
    fig3.savefig("mutational_signatures_context.png")
    print("   Saved: mutational_signatures_context.png")
    
    # Example 4: Load and use COSMIC signatures (if available)
    cosmic_file = "COSMIC_signatures.tsv"
    if os.path.exists(cosmic_file):
        print(f"\n4. Analysis with COSMIC signatures from {cosmic_file}...")
        cosmic_df = pd.read_csv(cosmic_file, sep='\t', index_col=0)
        
        fig4 = py_mut.mutational_signature_plot(
            n_signatures=4,
            cosmic_signatures=cosmic_df,
            figsize=(24, 28),
            title="Mutational Signatures with COSMIC Comparison"
        )
        fig4.savefig("mutational_signatures_cosmic.png")
        print("   Saved: mutational_signatures_cosmic.png")
    else:
        print(f"\n4. COSMIC signatures file not found ({cosmic_file})")
        print("   The analysis will use synthetic COSMIC signatures for demonstration")
    
    # Example 5: Interactive visualization
    print("\n5. Interactive visualization (will open in new window)...")
    fig5 = py_mut.mutational_signature_plot(
        n_signatures=3,
        figsize=(18, 20),
        show_interactive=True,
        title="Interactive Mutational Signature Analysis"
    )
    
    print("\n✅ All examples completed successfully!")
    print("\nGenerated files:")
    print("  - mutational_signatures_basic.png")
    print("  - mutational_signatures_5sig.png")
    print("  - mutational_signatures_context.png")
    if os.path.exists(cosmic_file):
        print("  - mutational_signatures_cosmic.png")
    
    return py_mut, data


def generate_synthetic_mutation_data(n_samples=50, n_mutations=5000):
    """
    Generate synthetic mutation data for demonstration purposes.
    
    This creates realistic-looking mutation data with:
    - SNV mutations (single nucleotide variants)
    - Trinucleotide contexts
    - Multiple samples
    - Various genes
    """
    np.random.seed(42)  # For reproducibility
    
    # Sample names
    samples = [f"Sample_{i:03d}" for i in range(1, n_samples + 1)]
    
    # Common cancer genes
    genes = ['TP53', 'KRAS', 'EGFR', 'BRAF', 'PIK3CA', 'PTEN', 'APC', 
             'BRCA1', 'BRCA2', 'MYC', 'RB1', 'VHL', 'MLH1', 'MSH2', 
             'CDKN2A', 'NOTCH1', 'JAK2', 'FLT3', 'IDH1', 'IDH2']
    
    # Nucleotides
    nucleotides = ['A', 'C', 'G', 'T']
    
    # Generate mutations
    mutations = []
    
    for _ in range(n_mutations):
        # Random sample
        sample = np.random.choice(samples)
        
        # Random gene
        gene = np.random.choice(genes)
        
        # Random reference and alt alleles (ensuring they're different)
        ref = np.random.choice(nucleotides)
        alt = np.random.choice([n for n in nucleotides if n != ref])
        
        # Generate trinucleotide context
        upstream = np.random.choice(nucleotides)
        downstream = np.random.choice(nucleotides)
        context = f"{upstream}{ref}{downstream}"
        
        # Random position
        position = np.random.randint(1000000, 100000000)
        
        # Variant classification
        variant_classes = ['Missense_Mutation', 'Nonsense_Mutation', 
                          'Silent', 'Frame_Shift_Del', 'Frame_Shift_Ins']
        variant_class = np.random.choice(variant_classes, p=[0.6, 0.15, 0.15, 0.05, 0.05])
        
        # Variant type
        variant_type = 'SNP'  # For this example, we're focusing on SNPs
        
        mutations.append({
            'Hugo_Symbol': gene,
            'Chromosome': np.random.choice([str(i) for i in range(1, 23)] + ['X', 'Y']),
            'Start_Position': position,
            'End_Position': position,
            'Variant_Classification': variant_class,
            'Variant_Type': variant_type,
            'Reference_Allele': ref,
            'Tumor_Seq_Allele2': alt,
            'REF': ref,
            'ALT': alt,
            'Tumor_Sample_Barcode': sample,
            'HGVSp': f'p.{ref}{position}{alt}',
            'Transcript_ID': f'ENST{np.random.randint(100000, 999999)}',
            'RefSeq': f'NM_{np.random.randint(100000, 999999)}',
            'Protein_position': f'{position // 3}',
            'Codons': f'{ref}NN',
            'Reference_Context': context,
            'Strand': np.random.choice(['+', '-'])
        })
    
    return pd.DataFrame(mutations)


def generate_trinucleotide_contexts(data):
    """
    Generate trinucleotide contexts for mutations.
    
    In real data, this would come from the reference genome.
    For demonstration, we generate random contexts.
    """
    nucleotides = ['A', 'C', 'G', 'T']
    contexts = []
    
    for _, row in data.iterrows():
        ref = row.get('REF', row.get('Reference_Allele', 'N'))
        if ref in nucleotides:
            upstream = np.random.choice(nucleotides)
            downstream = np.random.choice(nucleotides)
            context = f"{upstream}{ref}{downstream}"
        else:
            context = "NNN"
        contexts.append(context)
    
    return contexts


if __name__ == "__main__":
    # Run the examples
    py_mut, data = main()
    
    # Additional analysis can be done here
    print("\n📊 Additional Information:")
    print(f"Total mutations analyzed: {len(data)}")
    print(f"Total samples: {data['Tumor_Sample_Barcode'].nunique()}")
    print(f"Most mutated genes: {data['Hugo_Symbol'].value_counts().head(5).to_dict()}")
    print(f"Variant types: {data['Variant_Type'].value_counts().to_dict()}") 