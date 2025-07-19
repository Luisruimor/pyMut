import sys
import os
import pandas as pd

# Add src to path to import pyMut
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from pyMut import PyMutation


def main():
    """Generate oncoplots from example data."""
    print("=" * 60)
    print("ONCOPLOT EXAMPLE - pyMut")
    print("=" * 60)
    
    # Load example data
    data_path = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyMut', 'data', 'examples', 'tcga_laml_converted.tsv')
    print(f"Loading data from: {data_path}")
    data = pd.read_csv(data_path, sep='\t')
    print(f"Data loaded: {data.shape[0]} rows, {data.shape[1]} columns")
    
    # Configure high quality plots
    PyMutation.configure_high_quality_plots()
    
    # Create PyMutation object
    py_mut = PyMutation(data)
    
    # Create output directory
    output_dir = os.path.join(os.path.dirname(__file__), 'oncoplot')
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")
    
    # Generate oncoplot with 50 samples limit
    print("\nGenerating oncoplot (50 samples)...")
    fig1 = py_mut.oncoplot(
        title="TCGA LAML Oncoplot (Top 50 samples)",
        top_genes_count=15,
        max_samples=50,
        show_interactive=False
    )
    
    # Save limited oncoplot
    fig1.savefig(os.path.join(output_dir, "oncoplot_50_samples.png"))
    fig1.savefig(os.path.join(output_dir, "oncoplot_50_samples.pdf"))
    print("✓ Oncoplot (50 samples) saved")
    
    # Generate natural oncoplot without sample limit
    print("\nGenerating natural oncoplot (all samples)...")
    fig2 = py_mut.oncoplot(
        title="TCGA LAML Oncoplot (All samples)",
        top_genes_count=15,
        max_samples=None,
        show_interactive=False
    )
    
    # Save natural oncoplot
    fig2.savefig(os.path.join(output_dir, "oncoplot_natural.png"))
    fig2.savefig(os.path.join(output_dir, "oncoplot_natural.pdf"))
    print("✓ Natural oncoplot saved")
    
    print("\n" + "=" * 60)
    print("FILES GENERATED:")
    print("=" * 60)
    for file in sorted(os.listdir(output_dir)):
        if file.endswith(('.png', '.pdf')):
            filepath = os.path.join(output_dir, file)
            size_mb = os.path.getsize(filepath) / (1024 * 1024)
            print(f"• {file:<25} ({size_mb:.1f} MB)")
    print("=" * 60)


if __name__ == "__main__":
    main() 