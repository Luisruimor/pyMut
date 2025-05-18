"""
Basic example for mutation visualization with pyMut.

This example shows how to generate genetic mutation data visualizations using
pyMut, showing both the complete summary and individual visualizations in the simplest
possible way. Each plot method is shown with default parameters first, followed by
a commented version with all available parameters.
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

# Add src directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import PyMutation
from src.pyMut import PyMutation

def main():
    """
    Main script for generating genetic mutation visualizations.
    """
    # Get the absolute path to the project directory
    project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    
    # Path to the example file
    tsv_file = os.path.join(project_dir, "src", "pyMut", "data", "examples", "tcga_laml_converted.tsv")
    
    # Check if the file exists
    if not os.path.exists(tsv_file):
        print(f"Error: The file {tsv_file} does not exist.")
        exit(1)
    
    # Read the data
    mutation_data = pd.read_csv(tsv_file, sep='\t')
    
    # Creating PyMutation object
    print("Creating PyMutation object...")
    py_mut = PyMutation(mutation_data)
    
    # 1. Generate the complete summary plot
    print("\n1. Generating complete summary plot...")
    summary_fig = py_mut.summary_plot()
    # To view this visualization interactively:
    # py_mut.show_interactive_plots([summary_fig])
    
    # Commented version with all available parameters:
    # summary_fig = py_mut.summary_plot(
    #     title="Example summary plot",     # Custom title for the plot
    #     figsize=(16, 12),                 # Figure size in inches (width, height)
    #     max_samples=50,                   # Maximum number of samples to display
    #     top_genes_count=5,                # Number of genes to show in Top Mutated Genes (default is 10)
    #     show_interactive=True             # Whether to show the plot interactively
    # )
    
    # Save the figure
    summary_fig.savefig('examples/summary_plot.png')
    # Additional options for saving:
    # summary_fig.savefig(output_path, dpi=300, bbox_inches='tight')
    
    #########################################################################################################
    
    # 2. Variant Classification Plot
    print("\n2. Generating variant classification plot...")
    vc_fig = py_mut.variant_classification_plot()
    # To view this visualization interactively, uncomment:
    # py_mut.show_interactive_plots([vc_fig])
    
    # Commented version with all available parameters:
    # vc_fig = py_mut.variant_classification_plot(
    #     figsize=(12, 6),                  # Figure size in inches (width, height)
    #     title="Variant Classification",   # Custom title for the plot
    #     show_interactive=True             # Whether to show the plot interactively
    # )

    # Save the figure    
    vc_fig.savefig('examples/variant_classification.png')
    
    #########################################################################################################
    
    # 3. Variant Type Plot
    print("\n3. Generating variant types plot...")
    vt_fig = py_mut.variant_type_plot()
    # To view this visualization interactively, uncomment:
    # py_mut.show_interactive_plots([vt_fig])
    
    # Commented version with all available parameters:
    # vt_fig = py_mut.variant_type_plot(
    #     figsize=(12, 6),                  # Figure size in inches (width, height)
    #     title="Variant Type",             # Custom title for the plot
    #     show_interactive=True             # Whether to show the plot interactively
    # )
    
    # Save the figure
    vt_fig.savefig('examples/variant_type.png')
    
    #########################################################################################################
    
    # 4. SNV Class Plot
    print("\n4. Generating SNV classes plot...")
    snv_fig = py_mut.snv_class_plot()
    # To view this visualization interactively, uncomment:
    # py_mut.show_interactive_plots([snv_fig])
    
    # Commented version with all available parameters:
    # snv_fig = py_mut.snv_class_plot(
    #     figsize=(12, 6),                  # Figure size in inches (width, height)
    #     title="SNV Class",                # Custom title for the plot
    #     ref_column="REF",                 # Column containing reference allele
    #     alt_column="ALT",                 # Column containing alternate allele
    #     show_interactive=True             # Whether to show the plot interactively
    # )
    
    # Save the figure
    snv_fig.savefig('examples/snv_class.png')
    
    #########################################################################################################
    
    # 5. Variants per Sample Plot (TMB)
    print("\n5. Generating variants per sample plot (TMB)...")
    tmb_fig = py_mut.variants_per_sample_plot()
    # To view this visualization interactively, uncomment:
    # py_mut.show_interactive_plots([tmb_fig])
    
    # Commented version with all available parameters:
    # tmb_fig = py_mut.variants_per_sample_plot(
    #     figsize=(12, 6),                      # Figure size in inches (width, height)
    #     max_samples=50,                       # Maximum number of samples to display
    #     title="Variants per Sample",          # Custom title for the plot
    #     variant_column="Variant_Classification", # Column with variant classification
    #     sample_column="Tumor_Sample_Barcode", # Column with sample ID
    #     show_interactive=True                 # Whether to show the plot interactively
    # )
    
    # Save the figure
    tmb_fig.savefig('examples/variants_per_sample.png')
    
    #########################################################################################################
    
    # 6. Variant Classification Summary Plot (Boxplot)
    print("\n6. Generating variant classification summary plot (Boxplot)...")
    vcs_fig = py_mut.variant_classification_summary_plot()
    # To view this visualization interactively, uncomment:
    # py_mut.show_interactive_plots([vcs_fig])
    
    # Commented version with all available parameters:
    # vcs_fig = py_mut.variant_classification_summary_plot(
    #     figsize=(12, 6),                      # Figure size in inches (width, height)
    #     title="Variant Classification Summary", # Custom title for the plot
    #     variant_column="Variant_Classification", # Column with variant classification
    #     sample_column="Tumor_Sample_Barcode", # Column with sample ID
    #     show_interactive=True                 # Whether to show the plot interactively
    # )
    
    # Save the figure
    vcs_fig.savefig('examples/variant_classification_summary.png')
    
    #########################################################################################################
    
    # 7a. Top Mutated Genes Plot (variants mode)
    print("\n7a. Generating top mutated genes plot (variants mode)...")
    tmg_variants_fig = py_mut.top_mutated_genes_plot(mode="variants")
    # To view this visualization interactively, uncomment:
    # py_mut.show_interactive_plots([tmg_variants_fig])
    
    # Commented version with all available parameters:
    # tmg_variants_fig = py_mut.top_mutated_genes_plot(
    #     mode="variants",                      # Counting mode (required): "variants" or "samples"
    #     figsize=(12, 6),                      # Figure size in inches (width, height)
    #     title="Top Mutated Genes",            # Custom title for the plot
    #     variant_column="Variant_Classification", # Column with variant classification
    #     gene_column="Hugo_Symbol",            # Column with gene symbol
    #     sample_column="Tumor_Sample_Barcode", # Column with sample ID
    #     count=10,                             # Number of genes to show (default 10, shows all if less)
    #     show_interactive=True                 # Whether to show the plot interactively
    # )
    
    # Save the figure
    tmg_variants_fig.savefig('examples/top_mutated_genes_variants.png')
    
    #########################################################################################################
    
    # 7b. Top Mutated Genes Plot (samples mode)
    print("\n7b. Generating top mutated genes plot (samples mode)...")
    tmg_samples_fig = py_mut.top_mutated_genes_plot(mode="samples")
    # To view this visualization interactively, uncomment:
    # py_mut.show_interactive_plots([tmg_samples_fig])
    
    # Commented version with all available parameters:
    # tmg_samples_fig = py_mut.top_mutated_genes_plot(
    #     mode="samples",                       # Counting mode (required): "variants" or "samples"
    #     figsize=(12, 6),                      # Figure size in inches (width, height)
    #     title="Top Mutated Genes",            # Custom title for the plot
    #     variant_column="Variant_Classification", # Column with variant classification
    #     gene_column="Hugo_Symbol",            # Column with gene symbol
    #     sample_column="Tumor_Sample_Barcode", # Column with sample ID
    #     count=15,                             # Number of genes to show
    #     show_interactive=True                 # Whether to show the plot interactively
    # )
    
    # Save the figure
    tmg_samples_fig.savefig('examples/top_mutated_genes_samples.png')


if __name__ == "__main__":
    main()
