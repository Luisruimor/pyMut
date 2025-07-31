"""
Basic example for mutation visualization with pyMut.

This example shows how to generate genetic mutation data visualizations using
pyMut, showing both the complete summary and individual visualizations in the simplest
possible way. Each plot method is shown with default parameters first, followed by
a commented version with all available parameters.

SAVING OPTIONS REFERENCE:
========================

This file demonstrates multiple ways to save figures with different quality settings:

1. AUTOMATIC HIGH QUALITY (RECOMMENDED):
   PyMutation.configure_high_quality_plots()  # Configure once at the beginning
   fig.savefig('plot.png')                    # Automatically DPI=300, bbox_inches='tight'

2. MANUAL HIGH QUALITY:
   fig.savefig('plot.png', dpi=300, bbox_inches='tight')

3. CENTRALIZED METHOD:
   py_mut.save_figure(fig, 'plot.png')        # Uses pyMut's centralized method

4. CUSTOM QUALITY:
   fig.savefig('plot.png', dpi=600, bbox_inches='tight', facecolor='white')

5. DIFFERENT FORMATS:
   fig.savefig('plot.pdf', dpi=300, bbox_inches='tight')  # PDF
   fig.savefig('plot.svg', bbox_inches='tight')           # SVG (vector)
   fig.savefig('plot.jpg', dpi=300, bbox_inches='tight')  # JPEG

DPI QUALITY LEVELS:
- 72-100: Screen quality (small files, pixelated when printed)
- 150-200: Good quality (medium files)
- 300: Publication quality (large files, crisp printing) ⭐ RECOMMENDED
- 600+: Ultra-high quality (very large files, professional printing)

BBOX_INCHES OPTIONS:
- None: Fixed margins (may waste space)
- 'tight': Automatic margins (removes unnecessary whitespace) ⭐ RECOMMENDED
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

# Import PyMutation directly (now installed in development mode)
from pyMut import PyMutation

def main():
    """
    Main script for generating genetic mutation visualizations.
    """
    # ========================================================================
    # STEP 0: Configure high quality once (RECOMMENDED)
    # ========================================================================
    print("⚙️  Configuring automatic high quality...")
    PyMutation.configure_high_quality_plots()
    print("✅ All figures will now be saved automatically in high quality!")
    
    # Path to the example file (using relative path from project root)
    tsv_file = "src/pyMut/data/examples/tcga_laml_converted.tsv"
    
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

    # Commented version with all available parameters:
    # summary_fig = py_mut.summary_plot(
    #     title="Example summary plot",     # Custom title for the plot
    #     figsize=(16, 12),                 # Figure size in inches (width, height)
    #     max_samples=50,                   # Maximum number of samples to display
    #     top_genes_count=5,                # Number of genes to show in Top Mutated Genes (default is 10)
    # )
    
    # Save the figure (automatically high quality!)
    summary_fig.savefig('examples/summary_plot.png')
    
    # Alternative saving methods (commented examples):
    # Method 1: Manual high quality parameters
    # summary_fig.savefig('examples/summary_plot.png', dpi=300, bbox_inches='tight')
    # 
    # Method 2: Using pyMut's centralized save_figure() method
    # py_mut.save_figure(summary_fig, 'examples/summary_plot.png')
    # 
    # Method 3: Custom quality settings
    # summary_fig.savefig('examples/summary_plot.png', dpi=600, bbox_inches='tight', 
    #                     facecolor='white', edgecolor='none', transparent=False)
    # 
    # Method 4: Different formats
    # summary_fig.savefig('examples/summary_plot.pdf', dpi=300, bbox_inches='tight')  # PDF
    # summary_fig.savefig('examples/summary_plot.svg', bbox_inches='tight')           # SVG (vector)
    # summary_fig.savefig('examples/summary_plot.jpg', dpi=300, bbox_inches='tight')  # JPEG
    
    #########################################################################################################
    
    # 2. Variant Classification Plot
    print("\n2. Generating variant classification plot...")
    vc_fig = py_mut.variant_classification_plot()

    # Commented version with all available parameters:
    # vc_fig = py_mut.variant_classification_plot(
    #     figsize=(12, 6),                  # Figure size in inches (width, height)
    #     title="Variant Classification",   # Custom title for the plot
    # )

    # Save the figure (automatically high quality!)
    vc_fig.savefig('examples/variant_classification.png')
    
    # Alternative saving methods (commented examples):
    # vc_fig.savefig('examples/variant_classification.png', dpi=300, bbox_inches='tight')
    # py_mut.save_figure(vc_fig, 'examples/variant_classification.png')
    # vc_fig.savefig('examples/variant_classification.pdf', dpi=300, bbox_inches='tight')
    
    #########################################################################################################
    
    # 3. Variant Type Plot
    print("\n3. Generating variant types plot...")
    vt_fig = py_mut.variant_type_plot()

    # Commented version with all available parameters:
    # vt_fig = py_mut.variant_type_plot(
    #     figsize=(12, 6),                  # Figure size in inches (width, height)
    #     title="Variant Type",             # Custom title for the plot
    # )
    
    # Save the figure (automatically high quality!)
    vt_fig.savefig('examples/variant_type.png')
    
    # Alternative saving methods (commented examples):
    # vt_fig.savefig('examples/variant_type.png', dpi=300, bbox_inches='tight')
    # py_mut.save_figure(vt_fig, 'examples/variant_type.png')
    # vt_fig.savefig('examples/variant_type.pdf', dpi=300, bbox_inches='tight')
    
    #########################################################################################################
    
    # 4. SNV Class Plot
    print("\n4. Generating SNV classes plot...")
    snv_fig = py_mut.snv_class_plot()

    # Commented version with all available parameters:
    # snv_fig = py_mut.snv_class_plot(
    #     figsize=(12, 6),                  # Figure size in inches (width, height)
    #     title="SNV Class",                # Custom title for the plot
    #     ref_column="REF",                 # Column containing reference allele
    #     alt_column="ALT",                 # Column containing alternate allele
    # )
    
    # Save the figure (automatically high quality!)
    snv_fig.savefig('examples/snv_class.png')
    
    # Alternative saving methods (commented examples):
    # snv_fig.savefig('examples/snv_class.png', dpi=300, bbox_inches='tight')
    # py_mut.save_figure(snv_fig, 'examples/snv_class.png')
    # snv_fig.savefig('examples/snv_class.pdf', dpi=300, bbox_inches='tight')
    
    #########################################################################################################
    
    # 5. Variants per Sample Plot (TMB)
    print("\n5. Generating variants per sample plot (TMB)...")
    tmb_fig = py_mut.variants_per_sample_plot()

    # Commented version with all available parameters:
    # tmb_fig = py_mut.variants_per_sample_plot(
    #     figsize=(12, 6),                      # Figure size in inches (width, height)
    #     max_samples=50,                       # Maximum number of samples to display
    #     title="Variants per Sample",          # Custom title for the plot
    #     variant_column="Variant_Classification", # Column with variant classification
    #     sample_column="Tumor_Sample_Barcode", # Column with sample ID
    # )
    
    # Save the figure (automatically high quality!)
    tmb_fig.savefig('examples/variants_per_sample.png')
    
    # Alternative saving methods (commented examples):
    # tmb_fig.savefig('examples/variants_per_sample.png', dpi=300, bbox_inches='tight')
    # py_mut.save_figure(tmb_fig, 'examples/variants_per_sample.png')
    # tmb_fig.savefig('examples/variants_per_sample.pdf', dpi=300, bbox_inches='tight')
    
    #########################################################################################################
    
    # 6. Variant Classification Summary Plot (Boxplot)
    print("\n6. Generating variant classification summary plot (Boxplot)...")
    vcs_fig = py_mut.variant_classification_summary_plot()

    # Commented version with all available parameters:
    # vcs_fig = py_mut.variant_classification_summary_plot(
    #     figsize=(12, 6),                      # Figure size in inches (width, height)
    #     title="Variant Classification Summary", # Custom title for the plot
    #     variant_column="Variant_Classification", # Column with variant classification
    #     sample_column="Tumor_Sample_Barcode", # Column with sample ID
    # )
    
    # Save the figure (automatically high quality!)
    vcs_fig.savefig('examples/variant_classification_summary.png')
    
    # Alternative saving methods (commented examples):
    # vcs_fig.savefig('examples/variant_classification_summary.png', dpi=300, bbox_inches='tight')
    # py_mut.save_figure(vcs_fig, 'examples/variant_classification_summary.png')
    # vcs_fig.savefig('examples/variant_classification_summary.pdf', dpi=300, bbox_inches='tight')
    
    #########################################################################################################
    
    # 7a. Top Mutated Genes Plot (variants mode)
    print("\n7a. Generating top mutated genes plot (variants mode)...")
    tmg_variants_fig = py_mut.top_mutated_genes_plot(mode="variants")
    
    # Commented version with all available parameters:
    # tmg_variants_fig = py_mut.top_mutated_genes_plot(
    #     mode="variants",                      # Counting mode (required): "variants" or "samples"
    #     figsize=(12, 6),                      # Figure size in inches (width, height)
    #     title="Top Mutated Genes",            # Custom title for the plot
    #     variant_column="Variant_Classification", # Column with variant classification
    #     gene_column="Hugo_Symbol",            # Column with gene symbol
    #     sample_column="Tumor_Sample_Barcode", # Column with sample ID
    #     count=10,                             # Number of genes to show (default 10, shows all if less)
    # )
    
    # Save the figure (automatically high quality!)
    tmg_variants_fig.savefig('examples/top_mutated_genes_variants.png')
    
    # Alternative saving methods (commented examples):
    # tmg_variants_fig.savefig('examples/top_mutated_genes_variants.png', dpi=300, bbox_inches='tight')
    # py_mut.save_figure(tmg_variants_fig, 'examples/top_mutated_genes_variants.png')
    # tmg_variants_fig.savefig('examples/top_mutated_genes_variants.pdf', dpi=300, bbox_inches='tight')
    
    #########################################################################################################
    
    # 7b. Top Mutated Genes Plot (samples mode)
    print("\n7b. Generating top mutated genes plot (samples mode)...")
    tmg_samples_fig = py_mut.top_mutated_genes_plot(mode="samples")

    # Commented version with all available parameters:
    # tmg_samples_fig = py_mut.top_mutated_genes_plot(
    #     mode="samples",                       # Counting mode (required): "variants" or "samples"
    #     figsize=(12, 6),                      # Figure size in inches (width, height)
    #     title="Top Mutated Genes",            # Custom title for the plot
    #     variant_column="Variant_Classification", # Column with variant classification
    #     gene_column="Hugo_Symbol",            # Column with gene symbol
    #     sample_column="Tumor_Sample_Barcode", # Column with sample ID
    #     count=15,                             # Number of genes to show
    # )
    
    # Save the figure (automatically high quality!)
    tmg_samples_fig.savefig('examples/top_mutated_genes_samples.png')
    
    # Alternative saving methods (commented examples):
    # tmg_samples_fig.savefig('examples/top_mutated_genes_samples.png', dpi=300, bbox_inches='tight')
    # py_mut.save_figure(tmg_samples_fig, 'examples/top_mutated_genes_samples.png')
    # tmg_samples_fig.savefig('examples/top_mutated_genes_samples.pdf', dpi=300, bbox_inches='tight')
    
    # Advanced saving examples:
    # Ultra-high quality for publications
    # py_mut.save_figure(tmg_samples_fig, 'examples/top_mutated_genes_samples_ultra.png', dpi=600)
    # 
    # Multiple formats at once
    # for fmt in ['png', 'pdf', 'svg']:
    #     tmg_samples_fig.savefig(f'examples/top_mutated_genes_samples.{fmt}', 
    #                            dpi=300 if fmt != 'svg' else None, bbox_inches='tight')


if __name__ == "__main__":
    main()