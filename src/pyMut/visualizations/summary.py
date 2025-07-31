"""
Module for generating summary charts.

This module contains functions for creating summary visualizations
that show different statistics from mutation data.
"""

import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, Optional, Tuple
from ..core import PyMutation


def _create_variant_classification_plot(py_mut: PyMutation,
                                     variant_column: str = "Variant_Classification",
                                     ax: Optional[plt.Axes] = None,
                                     color_map: Optional[Dict] = None,
                                     set_title: bool = True) -> plt.Axes:
    """
    Create a horizontal bar chart showing the count for each type of variant classification.
    
    Args:
        py_mut: PyMutation object with mutation data.
        variant_column: Name of the column containing the variant classification.
        ax: Matplotlib axis to draw on. If None, a new one is created.
        color_map: Optional dictionary mapping variant classifications to colors.
        set_title: Whether to set the title on the plot.
        
    Returns:
        Matplotlib axis with the visualization.
    """
    data = py_mut.data

    # Count variants by classification
    variant_counts = data[variant_column].value_counts().to_dict()
    
    # Sort counts
    variant_counts = dict(sorted(variant_counts.items(), key=lambda item: item[1]))
    
    # Create axis if none was provided
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))
    
    # Color each bar using a color palette
    if color_map:
        colors = [color_map.get(variant, plt.colormaps['tab20'](i % 20)) for i, variant in enumerate(variant_counts.keys())]
    else:
        # Use the 'tab20' colormap
        cmap = plt.colormaps['tab20']  # Instead of cm.get_cmap('tab20')
        colors = [cmap(i % 20) for i in range(len(variant_counts))]
    
    bars = ax.barh(list(variant_counts.keys()), list(variant_counts.values()), color=colors)
    
    # Adjust title and labels, only if set_title is True
    if set_title:
        ax.set_title("Variant Classification", fontsize=14, fontweight='bold')
    ax.set_ylabel("")  # Remove Y-axis title
    
    # Add labels with count on each bar
    for bar in bars:
        ax.text(bar.get_width() + 10,
                 bar.get_y() + bar.get_height()/2,
                 f'{int(bar.get_width())}',
                 va='center', fontsize=10)
    
    # Adjust margins and remove unnecessary borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    return ax


def _create_variant_type_plot(py_mut: PyMutation,
                           variant_column: str = "Variant_Type",
                           ax: Optional[plt.Axes] = None,
                           set_title: bool = True) -> plt.Axes:
    """
    Create a horizontal bar chart showing the count for each variant type.
    
    Args:
        py_mut: PyMutation object with mutation data.
        variant_column: Name of the column containing the variant type.
        ax: Matplotlib axis to draw on. If None, a new one is created.
        set_title: Whether to set the title on the plot.
        
    Returns:
        Matplotlib axis with the visualization.
    """

    data = py_mut.data

    # Count variants by type
    variant_counts = data[variant_column].value_counts().to_dict()
    
    # Sort counts
    variant_counts = dict(sorted(variant_counts.items(), key=lambda item: item[1]))
    
    # Create axis if none was provided
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))
    
    # Specific colors (if number of categories exceeds the list, a colormap will be used)
    colors = ['#D3C4E7', '#FFFACD', '#87CEEB']
    if len(variant_counts) > len(colors):
        cmap = plt.colormaps['tab20']
        colors = cmap(range(len(variant_counts)))
    else:
        colors = colors[:len(variant_counts)]
    
    # Create horizontal bar chart
    bars = ax.barh(list(variant_counts.keys()), list(variant_counts.values()), color=colors)
    
    # Adjust title and labels, only if set_title is True
    if set_title:
        ax.set_title("Variant Type", fontsize=14, fontweight='bold')
    ax.set_ylabel("")  # Remove Y-axis title
    
    # Y-axis labels in normal style (not italic)
    # Italic style configuration has been removed
    
    # Add labels with count on each bar
    for bar in bars:
        ax.text(bar.get_width() + 10,
                 bar.get_y() + bar.get_height()/2,
                 f'{int(bar.get_width())}',
                 va='center', fontsize=10)
    
    # Adjust margins and remove unnecessary borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    return ax


def _create_snv_class_plot(py_mut: PyMutation,
                        ref_column: str = "REF",
                        alt_column: str = "ALT",
                        ax: Optional[plt.Axes] = None,
                        set_title: bool = True) -> plt.Axes:
    """
    Create a horizontal bar chart showing the count for each SNV class.
    
    Args:
        py_mut: PyMutation object with mutation data.
        ref_column: Name of the column containing the reference allele.
        alt_column: Name of the column containing the alternative (tumor) allele.
        ax: Matplotlib axis to draw on. If None, a new one is created.
        set_title: Whether to set the title on the plot.
        
    Returns:
        Matplotlib axis with the visualization.
    """
    data = py_mut.data

    # Create a copy of the DataFrame to not modify the original
    df_copy = data.copy()
    
    # Check if we can generate SNV Class information
    if ref_column in df_copy.columns and alt_column in df_copy.columns:
        # Generate SNV Class by combining reference and alternative alleles
        df_copy['SNV_Class'] = df_copy[ref_column] + '>' + df_copy[alt_column]
        # Filter to keep only single nucleotide changes (one character on each side)
        df_copy = df_copy[df_copy['SNV_Class'].str.match(r'^[A-Z]>[A-Z]$')]
    else:
        # If we can't generate the information, return an empty axis
        print(f"Cannot generate SNV classes: columns {ref_column} and/or {alt_column} are missing.")
        if ax is None:
            _, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No data available for SNV Class", 
               ha='center', va='center', fontsize=12)
        if set_title:
            ax.set_title("SNV Class", fontsize=14, fontweight='bold')
        ax.axis('off')
        return ax
    
    # Remove rows without SNV Class information
    df_copy = df_copy[df_copy['SNV_Class'].notnull()]
    
    # Count variants by SNV class
    snv_counts = df_copy['SNV_Class'].value_counts().to_dict()
    
    # If there's no data, show a message
    if not snv_counts:
        if ax is None:
            _, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No data available for SNV Class", 
               ha='center', va='center', fontsize=12)
        if set_title:
            ax.set_title("SNV Class", fontsize=14, fontweight='bold')
        ax.axis('off')
        return ax
    
    # Sort counts
    snv_counts = dict(sorted(snv_counts.items(), key=lambda item: item[1]))
    
    # Create axis if none was provided
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))
    
    # Define specific colors for each class
    colors = ['#FF8C00', '#9ACD32', '#FFD700', '#FF4500', '#4169E1', '#1E90FF']
    if len(snv_counts) > len(colors):
        cmap = plt.colormaps['tab20']
        colors = cmap(range(len(snv_counts)))
    else:
        colors = colors[:len(snv_counts)]
    
    # Create horizontal bar chart
    bars = ax.barh(list(snv_counts.keys()), list(snv_counts.values()), color=colors)
    
    # Adjust title and labels, only if set_title is True
    if set_title:
        ax.set_title("SNV Class", fontsize=14, fontweight='bold')
    ax.set_ylabel("")  # Remove Y-axis title
    
    # Add labels with count on each bar
    for bar in bars:
        ax.text(bar.get_width() + 10,
                 bar.get_y() + bar.get_height()/2,
                 f'{int(bar.get_width())}',
                 va='center', fontsize=10)
    
    # Adjust margins and remove unnecessary borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    return ax


def _create_variant_classification_summary_plot(py_mut: PyMutation,
                                             variant_column: str = "Variant_Classification",
                                             sample_column: str = "Tumor_Sample_Barcode",
                                             ax: Optional[plt.Axes] = None,
                                             color_map: Optional[Dict] = None,
                                             show_labels: bool = True,
                                             set_title: bool = True) -> plt.Axes:
    """
    Create a box-and-whiskers plot (boxplot) that summarizes, for each variant classification,
    the distribution (among samples) of the number of detected alternative alleles.

    Args:
        py_mut: PyMutation object with mutation data.
        variant_column: Name of the column containing the variant classification.
        sample_column: Name of the column containing the sample identifier.
                       If it doesn't exist, samples are assumed to be columns (wide format).
        ax: Matplotlib axis to draw on. If None, a new one is created.
        color_map: Optional dictionary mapping variant classifications to colors.
        show_labels: If True, shows the classification labels on the X-axis.
        set_title: Whether to set the title on the plot.
        
    Returns:
        Matplotlib axis with the visualization.
    """
    data = py_mut.data

    # Check if we have the classification column
    if variant_column not in data.columns:
        print(f"Column not found: {variant_column}")
        if ax is None:
            _, ax = plt.subplots(figsize=(12, 6))
        ax.text(0.5, 0.5, f"No data available for Variant Classification Summary\nMissing column: {variant_column}", 
               ha='center', va='center', fontsize=12)
        if set_title:
            ax.set_title("Variant Classification Summary", fontsize=14, fontweight='bold')
        ax.axis('off')
        return ax
    
    # Detect if we have long or wide format
    samples_as_columns = sample_column not in data.columns
    
    if samples_as_columns:
        # Wide format: samples are columns
        # Find columns that could be samples (TCGA IDs, etc.)
        potential_sample_cols = [col for col in data.columns if col.startswith('TCGA-') or 
                                (isinstance(col, str) and col.count('-') >= 2)]
        
        if not potential_sample_cols:
            print("No sample columns found that look like identifiers")
            if ax is None:
                _, ax = plt.subplots(figsize=(12, 6))
            ax.text(0.5, 0.5, "No data available for Variant Classification Summary\nNo sample columns detected", 
                  ha='center', va='center', fontsize=12)
            if set_title:
                ax.set_title("Variant Classification Summary", fontsize=14, fontweight='bold')
            ax.axis('off')
            return ax
            
        print(f"Detected {len(potential_sample_cols)} sample columns")
        
        # Accumulate counts by sample and classification
        sample_variant_counts = {}
        
        # Get unique variant values for grouping
        unique_variants = data[variant_column].unique()
        unique_variants = [v for v in unique_variants if pd.notna(v) and v != "Unknown"]
        
        # First determine the format of values in sample columns
        # Check some rows to see if they are of type "A|B"
        sample_col = potential_sample_cols[0]
        sample_format = "unknown"
        sample_values = data[sample_col].dropna().values[:10]  # Take some samples
        
        if len(sample_values) > 0:
            # Check if there are values in "A|B" or "A/B" format
            if any('|' in str(x) for x in sample_values):
                sample_format = "pipe_separated"
            elif any('/' in str(x) for x in sample_values):
                sample_format = "slash_separated"
            else:
                sample_format = "other"
                
        print(f"Detected sample format: {sample_format}")
        
        # For each variant classification, count how many samples have that variant
        for variant_class in unique_variants:
            # Filter rows with this variant classification
            variant_subset = data[data[variant_column] == variant_class]
            
            # For each sample, count how many variants of this type it has
            for sample in potential_sample_cols:
                # Count according to detected format
                if sample_format == "pipe_separated":
                    # For "A|B" format, count when A and B are different (variant allele)
                    sample_count = variant_subset[sample].apply(
                        lambda x: 1 if (isinstance(x, str) and '|' in x and 
                                       x.split('|')[0] != x.split('|')[1]) else 0
                    ).sum()
                elif sample_format == "slash_separated":
                    # Similar for "A/B" format
                    sample_count = variant_subset[sample].apply(
                        lambda x: 1 if (isinstance(x, str) and '/' in x and 
                                       x.split('/')[0] != x.split('/')[1]) else 0
                    ).sum()
                else:
                    # For other formats, assume values other than "0", "0/0" or "0|0" indicate variant
                    sample_count = variant_subset[sample].apply(
                        lambda x: 1 if (x != 0 and x != '0' and x != '0/0' and 
                                       x != '0|0' and not pd.isnull(x)) else 0
                    ).sum()
                
                if sample not in sample_variant_counts:
                    sample_variant_counts[sample] = {}
                
                if sample_count > 0:
                    sample_variant_counts[sample][variant_class] = sample_count
    else:
        # Long format: there is a specific column for the sample
        # Accumulate counts by sample and classification
        sample_variant_counts = {}
        
        # Process DataFrame to count by sample and variant
        for _, row in data.iterrows():
            sample = row[sample_column]
            variant_class = row[variant_column]
            
            if pd.isnull(variant_class) or variant_class == "Unknown":
                continue
                
            if sample not in sample_variant_counts:
                sample_variant_counts[sample] = {}
                
            if variant_class not in sample_variant_counts[sample]:
                sample_variant_counts[sample][variant_class] = 0
                
            # Increment counter (assuming 1 alternative allele per variant)
            sample_variant_counts[sample][variant_class] += 1
    
    # Convert dictionary to DataFrame (rows: samples, columns: variant classifications)
    df_data = {sample: dict(variants) for sample, variants in sample_variant_counts.items()}
    df = pd.DataFrame.from_dict(df_data, orient='index').fillna(0)
    
    if df.empty:
        print("No data found for analysis after processing.")
        if ax is None:
            _, ax = plt.subplots(figsize=(12, 6))
        ax.text(0.5, 0.5, "No data available for Variant Classification Summary\nNo data found for analysis", 
               ha='center', va='center', fontsize=12)
        if set_title:
            ax.set_title("Variant Classification Summary", fontsize=14, fontweight='bold')
        ax.axis('off')
        return ax
    
    # (Optional) Reorder columns by total sum (highest to lowest)
    col_order = df.sum(axis=0).sort_values(ascending=False).index.tolist()
    df = df[col_order]
    
    # Remove columns with total sum 0
    df = df.loc[:, df.sum() > 0]
    
    # Prepare data for boxplot: each column (variant classification) is a series of counts per sample
    variant_types = df.columns.tolist()
    data_to_plot = [df[vt].values for vt in variant_types]
    
    # Create axis if none was provided
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 6))
    
    # Draw boxplot with improved configuration
    bp = ax.boxplot(data_to_plot, 
                    patch_artist=True,  # Fill boxes with color
                    medianprops=dict(color="red", linewidth=1.5),
                    boxprops=dict(linewidth=1.5),
                    whiskerprops=dict(linewidth=1.5),
                    capprops=dict(linewidth=1.5),
                    flierprops=dict(marker='o', markerfacecolor='gray', markersize=4, alpha=0.5),
                    showfliers=True,  # Show outliers
                    widths=0.7)
    
    # Color each box using a color palette or the provided color map
    if color_map:
        colors = [color_map.get(vt, plt.colormaps['tab20'](i % 20)) for i, vt in enumerate(variant_types)]
    else:
        # Use a high-differentiation colormap for boxes
        cmap = plt.colormaps['tab20']
        colors = [cmap(i % 20) for i in range(len(variant_types))]
    
    # Apply colors to each box
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Automatically detect if we're in summary plot or individual visualization
    # If ax is in a figure with multiple subplots, we're in summary plot
    is_in_summary_plot = hasattr(ax.figure, 'axes') and len(ax.figure.axes) > 1
    
    # Configure X-axis labels according to context
    if is_in_summary_plot:
        # In summary plot: don't show labels to avoid overlap
        ax.set_xticklabels([])
    else:
        # In individual visualization: always show labels rotated 45 degrees
        ax.set_xticklabels(variant_types, rotation=45, ha='right', fontsize=10)
    
    # Adjust Y-axis limits
    ymin = 0
    ymax = max(max(d) if len(d) > 0 else 0 for d in data_to_plot) * 1.1
    if ymax == 0:  # In case there's no positive data
        ymax = 1
    ax.set_ylim(ymin, ymax)
    
    # Improve visual presentation
    # Add more descriptive title, only if set_title is True
    if set_title:
        ax.set_title("Variant Classification Summary", fontsize=14, fontweight='bold')
    
    # Add grid for better readability
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)
    
    # Apply styling according to context
    if is_in_summary_plot:
        # In summary plot: apply modern style where lines don't cross
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_position(('outward', 5))
        ax.spines['bottom'].set_position(('outward', 5))
    else:
        # In individual visualization: use original classic style (axes that cross)
        # Only remove top and right borders, keep normal axis position
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # DO NOT use outward positioning - this allows axes to cross
    
    return ax


def _create_summary_plot(py_mut: PyMutation,
                      figsize: Tuple[int, int] = (16, 12),
                      title: str = "Mutation Summary",
                      max_samples: Optional[int] = 200,
                      top_genes_count: int = 10) -> plt.Figure:
    """
    Creates a summary plot with multiple visualizations of mutation data.
    
    Args:
        py_mut: PyMutation object with mutation data.
        figsize: Figure size.
        title: Plot title.
        max_samples: Maximum number of samples to show in the variants per sample plot.
                    If None, all samples are shown.
        top_genes_count: Number of genes to show in the top mutated genes plot.
                    If there are fewer genes than this number, all will be shown.
        
    Returns:
        Figure with summary visualizations.
    """
    data = py_mut.data

    # Create a figure with multiple subplots, making the charts wider
    fig, axs = plt.subplots(2, 3, figsize=figsize, gridspec_kw={'width_ratios': [1.5, 1.5, 1.5], 'height_ratios': [1, 1]})
    fig.suptitle(title, fontsize=16)
    
    # Detect column names respecting original capitalization
    variant_classification_col = "Variant_Classification"
    sample_column = "Tumor_Sample_Barcode"
    gene_column = "Hugo_Symbol"
    
    # Look for capitalization variants for variant_classification
    if variant_classification_col not in data.columns:
        for col in data.columns:
            if col.lower() == variant_classification_col.lower():
                variant_classification_col = col
                break
    
    # Look for capitalization variants for gene_column
    if gene_column not in data.columns:
        for col in data.columns:
            if col.lower() == gene_column.lower():
                gene_column = col
                break
    
    # Generate a coherent color map for all variant classifications
    unique_variants = data[variant_classification_col].unique()
    # Use a fixed colormap to ensure consistent colors
    cmap = plt.colormaps['tab20']  # Colormap with good variety of colors
    variant_color_map = {variant: cmap(i % 20) for i, variant in enumerate(unique_variants) if pd.notna(variant)}
    
    # Create the variant classification plot using the predefined color map
    var_class_ax = _create_variant_classification_plot(
        py_mut, 
        variant_column=variant_classification_col, 
        ax=axs[0, 0],
        color_map=variant_color_map,  # Pass the color map
        set_title=True  # Set title for summary
    )
    
    # Create the variant type plot
    _create_variant_type_plot(py_mut, ax=axs[0, 1], set_title=True)
    
    # Create the SNV class plot
    _create_snv_class_plot(py_mut, 
                         ref_column="REF",
                         alt_column="ALT",
                         ax=axs[0, 2],
                         set_title=True)
    
    # Create the variants per sample plot (TMB)
    # Passing the same color map used for the classification plot
    variants_ax = _create_variants_per_sample_plot(
        py_mut,
        variant_column=variant_classification_col,
        sample_column=sample_column,
        ax=axs[1, 0],
        color_map=variant_color_map,  # Use the same color map
        set_title=True,
        max_samples=max_samples  # Pass the max_samples parameter
    )
    
    # Create the variant classification summary plot
    var_boxplot_ax = _create_variant_classification_summary_plot(
        py_mut,
        variant_column=variant_classification_col,
        sample_column=sample_column,
        ax=axs[1, 1],
        color_map=variant_color_map,  # Use the same color map
        show_labels=False,  # Don't show labels in the summary plot
        set_title=True
    )
    
    # Create the top mutated genes plot
    top_genes_ax = _create_top_mutated_genes_plot(
        py_mut,
        variant_column=variant_classification_col,
        gene_column=gene_column,
        sample_column=sample_column,
        mode="variants",
        count=top_genes_count,  # Use the configurable parameter
        ax=axs[1, 2],
        color_map=variant_color_map,  # Use the same color map
        set_title=True
    )
    
    # Remove individual legends to avoid duplication
    if var_class_ax.get_legend() is not None:
        var_class_ax.get_legend().remove()
    
    if variants_ax.get_legend() is not None:
        variants_ax.get_legend().remove()
        
    if var_boxplot_ax.get_legend() is not None:
        var_boxplot_ax.get_legend().remove()
        
    if top_genes_ax.get_legend() is not None:
        top_genes_ax.get_legend().remove()
    
    # Create a common legend for the plots and place it at the bottom
    # Create handles and labels manually for the global legend
    handles = []
    labels = []
    
    # Count variants to order by frequency
    variant_counts = data[variant_classification_col].value_counts()
    
    # Sort variants by frequency (descending order)
    ordered_variants = variant_counts.index.tolist()
    
    # Create handles and labels in order of frequency
    for variant in ordered_variants:
        if variant in variant_color_map and not pd.isnull(variant) and variant != "Unknown":
            color = variant_color_map[variant]
            patch = plt.Rectangle((0,0), 1, 1, fc=color)
            handles.append(patch)
            labels.append(variant)
    
    # Add the common legend with variants ordered by frequency
    fig.legend(handles, labels, loc='lower center', ncol=min(len(labels), 5), 
               bbox_to_anchor=(0.5, 0.02))
    
    # Adjust spacing between subplots with padding to increase margins
    plt.tight_layout(pad=2.0)  # Increase pad for larger margin
    
    # Increase separation between first and second row, adjust left and right margins
    # Increase left substantially to prevent text from going off the left margin
    plt.subplots_adjust(top=0.9, bottom=0.15, hspace=0.4, left=0.15, right=0.93)
    
    return fig


def _create_variants_per_sample_plot(py_mut: PyMutation,
                                   variant_column: str = "Variant_Classification",
                                   sample_column: str = "Tumor_Sample_Barcode",
                                   ax: Optional[plt.Axes] = None,
                                   color_map: Optional[Dict] = None,
                                   set_title: bool = True,
                                   max_samples: Optional[int] = 200) -> plt.Axes:
    """
    Create a stacked bar plot showing the number of variants per sample (TMB)
    and their composition by variant type.

    Args:
        py_mut: PyMutation object with mutation data.
        variant_column: Name of the column containing the variant classification.
        sample_column: Name of the column containing the sample identifier,
                      or string used to identify sample columns if samples
                      are stored as columns.
        ax: Matplotlib axis to draw on. If None, a new one is created.
        color_map: Optional dictionary mapping variant classifications to colors.
        set_title: Whether to set the title on the plot.
        max_samples: Maximum number of samples to show. If None, all are shown.
                    If there are more samples than this number, only the first max_samples are shown.
        
    Returns:
        Matplotlib axis with the visualization.
    """
    data = py_mut.data

    # Check if we have the necessary columns
    if variant_column not in data.columns:
        print(f"Column not found: {variant_column}")
        if ax is None:
            _, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, f"No data available for Variants per Sample\nMissing column: {variant_column}", 
               ha='center', va='center', fontsize=12)
        if set_title:
            ax.set_title("Variants per Sample", fontsize=14, fontweight='bold')
        ax.axis('off')
        return ax
    
    # Detect data format
    # If sample_column exists as a column, assume "long" format
    # If not, assume "wide" format where samples are columns
    samples_as_columns = sample_column not in data.columns
    
    if samples_as_columns:
        # Find columns that could be samples (TCGA IDs, etc.)
        potential_sample_cols = [col for col in data.columns if col.startswith('TCGA-') or '|' in str(data[col].iloc[0])]
        
        if not potential_sample_cols:
            print("No sample columns found that start with TCGA- or contain the | character")
            if ax is None:
                _, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, "No data available for Variants per Sample\nNo sample columns detected", 
                  ha='center', va='center', fontsize=12)
            if set_title:
                ax.set_title("Variants per Sample", fontsize=14, fontweight='bold')
            ax.axis('off')
            return ax
        
        # Convert format "wide" to "long" for counting
        # First, count variants by type for each sample
        variant_counts = {}
        
        # Group by Variant_Classification
        for variant_type in data[variant_column].unique():
            variant_subset = data[data[variant_column] == variant_type]
            
            for sample_col in potential_sample_cols:
                # If the sample column contains the genotypes in REF|ALT format
                sample_variants = variant_subset[sample_col].apply(
                    lambda x: 1 if '|' in str(x) and str(x).split('|')[0] != str(x).split('|')[1] else 0
                ).sum()
                
                if sample_col not in variant_counts:
                    variant_counts[sample_col] = {}
                
                if sample_variants > 0:
                    variant_counts[sample_col][variant_type] = sample_variants
        
        # Create DataFrame from dictionary
        samples_df = []
        for sample, variants in variant_counts.items():
            for var_type, count in variants.items():
                samples_df.append({'Sample': sample, 'Variant_Classification': var_type, 'Count': count})
        
        if not samples_df:
            print("Unable to process variants per sample")
            if ax is None:
                _, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, "No data available for Variants per Sample\nUnable to process variants", 
                  ha='center', va='center', fontsize=12)
            if set_title:
                ax.set_title("Variants per Sample", fontsize=14, fontweight='bold')
            ax.axis('off')
            return ax
        
        processed_df = pd.DataFrame(samples_df)
        variant_counts = processed_df.pivot(index='Sample', columns='Variant_Classification', values='Count').fillna(0)
        
    else:
        # "long" format where there is a specific column for the sample
        # Count variants by sample and classification
        variant_counts = data.groupby([sample_column, variant_column]).size().unstack(fill_value=0)
    
    # Calculate total variants per sample and sort from highest to lowest
    variant_counts['total'] = variant_counts.sum(axis=1)
    variant_counts = variant_counts.sort_values('total', ascending=False)
    
    # Limit the number of samples if max_samples is specified
    if max_samples is not None and len(variant_counts) > max_samples:
        variant_counts = variant_counts.iloc[:max_samples]
    
    # Calculate median variants per sample
    median_tmb = variant_counts['total'].median()
    
    # Remove total column for visualization
    if 'total' in variant_counts.columns:
        variant_counts = variant_counts.drop('total', axis=1)
    
    # Create axis if none was provided
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))
    
    # Generate colors for different variant classifications
    if color_map is not None:
        # If specific colors are provided, use them for corresponding variants
        colors = []
        for variant in variant_counts.columns:
            # Find color in color map for this exact variant
            if variant in color_map:
                colors.append(color_map[variant])
            else:
                # If exact color not found, use a default one
                # Use plt.colormaps instead of cm.get_cmap
                colors.append(plt.colormaps['tab20'](len(colors) % 20))
    else:
        # Use default colormap
        # Use plt.colormaps instead of cm.get_cmap
        cmap = plt.colormaps['tab20']
        colors = [cmap(i % cmap.N) for i in range(len(variant_counts.columns))]
    
    # Create stacked bar plot
    variant_counts.plot(kind='bar', stacked=True, ax=ax, color=colors, width=0.8)
    
    # Add horizontal line for median
    ax.axhline(y=median_tmb, color='red', linestyle='--', linewidth=1)
    
    # Configure labels and title, only if set_title is True
    if set_title:
        ax.set_title("Variants per Sample", fontsize=14, fontweight='bold')
    
    ax.set_xlabel("")  # Remove X-axis label "Samples"
    
    # Remove X-axis labels for better clarity when there are many samples
    ax.set_xticklabels([])
    ax.tick_params(axis='x', which='both', bottom=False)
    
    # Automatically detect if we're in summary plot or individual visualization
    # If ax is in a figure with multiple subplots, we're in summary plot
    is_in_summary_plot = hasattr(ax.figure, 'axes') and len(ax.figure.axes) > 1
    
    if is_in_summary_plot:
        # In summary plot: use clean legend and add median below title
        ax.legend(title="Variant Classification", bbox_to_anchor=(1.05, 1), loc='upper left')
        # Add median below title in summary plot
        ax.text(0.5, 0.92, f"Median: {median_tmb:.1f}", transform=ax.transAxes, ha='center', fontsize=12)
    else:
        # In individual visualization: place median inside legend box
        # Include median in legend title for it to appear inside the box
        # Use rich text formatting to put "Variant Classification" and "Median:" in bold
        legend_title = f"$\\mathbf{{Variant Classification}}$\n\n$\\mathbf{{Median:}}$ {median_tmb:.1f}"
        ax.legend(title=legend_title, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Adjust margins and remove unnecessary borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    
    return ax


def _create_top_mutated_genes_plot(py_mut: PyMutation,
                               mode: str = "variants",
                               variant_column: str = "Variant_Classification",
                               gene_column: str = "Hugo_Symbol",
                               sample_column: str = "Tumor_Sample_Barcode",
                               count: int = 10,
                               ax: Optional[plt.Axes] = None,
                               color_map: Optional[Dict] = None,
                               set_title: bool = True) -> plt.Axes:
    """
    Create a horizontal bar plot showing the most mutated genes and distribution
    of variants according to their classification.

    Args:
        py_mut: PyMutation object with mutation data.
        mode: Mutation counting mode: "variants" (counts total number of variants)
              or "samples" (counts number of affected samples).
        variant_column: Name of the column containing the variant classification.
        gene_column: Name of the column containing the gene symbol.
        sample_column: Name of the column containing the sample identifier,
                      or prefix used to identify sample columns if samples
                      are stored as columns.
        count: Number of top genes to show.
        ax: Matplotlib axis to draw on. If None, a new one is created.
        color_map: Optional dictionary mapping variant classifications to colors.
        set_title: Whether to set the title on the plot.
        
    Returns:
        Matplotlib axis with the visualization.
    """
    data = py_mut.data

    # Check if we have the necessary columns
    if gene_column not in data.columns:
        print(f"Column not found: {gene_column}")
        if ax is None:
            _, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, f"No data available for Top Mutated Genes\nMissing column: {gene_column}", 
               ha='center', va='center', fontsize=12)
        if set_title:
            ax.set_title("Top Mutated Genes", fontsize=14, fontweight='bold')
        ax.axis('off')
        return ax
        
    if variant_column not in data.columns:
        print(f"Column not found: {variant_column}")
        if ax is None:
            _, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, f"No data available for Top Mutated Genes\nMissing column: {variant_column}", 
               ha='center', va='center', fontsize=12)
        if set_title:
            ax.set_title("Top Mutated Genes", fontsize=14, fontweight='bold')
        ax.axis('off')
        return ax
    
    # Detect if samples are stored as columns (wide format)
    # Search for columns that could be samples (e.g., start with TCGA-)
    sample_cols = [col for col in data.columns if str(col).startswith("TCGA-")]
    samples_as_columns = len(sample_cols) > 0
    
    if not samples_as_columns and sample_column not in data.columns:
        print(f"Sample column not found: {sample_column}")
        print("No sample columns with TCGA-* format detected either")
        if ax is None:
            _, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, "No data available for Top Mutated Genes\nNo samples detected", 
               ha='center', va='center', fontsize=12)
        if set_title:
            ax.set_title("Top Mutated Genes", fontsize=14, fontweight='bold')
        ax.axis('off')
        return ax
    
    # Create axis if none was provided
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 8))
    
    # Filter rows with missing or "Unknown" values
    data_filtered = data[(data[gene_column].notna()) & (data[variant_column].notna())]
    data_filtered = data_filtered[(data_filtered[gene_column] != "Unknown") & 
                                   (data_filtered[variant_column] != "Unknown")]
    
    if mode == "variants":
        # MODE "variants" - Count total number of variants
        # Count variants by gen and variant type
        gene_variant_counts = data_filtered.groupby([gene_column, variant_column]).size().unstack(fill_value=0)
        
        # Calculate total for each gen and sort
        gene_totals = gene_variant_counts.sum(axis=1).sort_values(ascending=False)
        top_genes = gene_totals.index[:count].tolist()
        
        # Select only top genes
        df_top = gene_variant_counts.loc[top_genes]
        
        # Manually sort them from lowest to highest for visualization
        ordered_genes = sorted(top_genes, key=lambda g: gene_totals[g])
        df_plot = df_top.loc[ordered_genes]
        
        # Assign colors for each variant type
        if color_map is None:
            cmap = plt.colormaps['tab20']  # Instead of cm.get_cmap
            colors = [cmap(i % 20) for i in range(len(df_plot.columns))]
        else:
            colors = [color_map.get(variant, plt.colormaps['tab20'](i % 20)) 
                     for i, variant in enumerate(df_plot.columns)]
        
        # Create horizontal bar plot - use exactly the same width in both modes
        df_plot.plot(kind='barh', stacked=True, ax=ax, color=colors, width=0.65)
        
        # Configure the same aspect and positioning for both modes
        ax.margins(x=0.05, y=0.01)
        
        # Add labels with total count to the right of each bar
        for i, gene in enumerate(df_plot.index):
            total = gene_totals[gene]
            # Adjust offset for numbers (smaller for being closer to bars)
            offset_variants = 0.01 * ax.get_xlim()[1] if ax.get_xlim()[1] > 0 else 0.1
            ax.text(total + offset_variants, i, f'{int(total)}', va='center', fontsize=10)
        
        title_text_variants = f"Top {count} Mutated Genes (variants)" # Use count variable
        # Common title and axis configuration for variants mode
        if set_title:
            ax.set_title(title_text_variants, fontsize=14, fontweight='bold') 
        ax.set_ylabel("")  # Remove Y-axis label
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False) # Hide left axis line
        ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=True) # No ticks on Y axis, but with labels
        
        # Configure the same exact margins in both modes for total consistency
        ax.margins(x=0.05, y=0.01)

        # Improve legend for variants mode
        handles_v, labels_v = ax.get_legend_handles_labels()
        by_label_v = dict(zip(labels_v, handles_v))
        
        if by_label_v:
            num_legend_items_v = len(by_label_v)
            ncol_legend_v = min(num_legend_items_v, 4)
            base_offset_v = -0.20
            row_offset_factor_v = 0.06
            num_legend_rows_v = (num_legend_items_v + ncol_legend_v - 1) // ncol_legend_v
            vertical_offset_v = base_offset_v - (row_offset_factor_v * num_legend_rows_v)
            ax.legend(by_label_v.values(), by_label_v.keys(),
                      title="Variant Classification",
                      loc="lower center",
                      bbox_to_anchor=(0.5, vertical_offset_v),
                      ncol=ncol_legend_v)
        elif ax.get_legend() is not None:
            ax.get_legend().remove()
            
        return ax
    
    else:  # mode == "samples"
        # MODE "samples" - Count number of affected samples by gen and variant type
        
        # gene_variant_sample_counts: Dict[str (gene), Dict[str (variant_type), Set[str (sample_id)]]]
        gene_variant_sample_counts = {}

        # Iterate over each unique gen present in filtered data
        for gene_iter_val in data_filtered[gene_column].unique():
            if pd.isna(gene_iter_val) or gene_iter_val == "Unknown":
                continue
            
            gene_variant_sample_counts[gene_iter_val] = {}
            # Filter DataFrame to get only rows corresponding to current gen
            gene_data_for_current_gene = data_filtered[data_filtered[gene_column] == gene_iter_val]

            if samples_as_columns:  # Wide format: samples are columns TCGA-*
                for sample_col_name in sample_cols:
                    # For each row (specific variant) of current gen
                    for _, row_series in gene_data_for_current_gene.iterrows():
                        sample_genotype_value = str(row_series[sample_col_name]).strip().upper()
                        actual_variant_type = row_series[variant_column]
                        
                        if pd.isna(actual_variant_type) or actual_variant_type == "Unknown":
                            continue

                        is_mutation_present_in_sample = False
                        if '|' in sample_genotype_value:
                            alleles = sample_genotype_value.split('|')
                            if len(alleles) >= 2 and alleles[0] != alleles[1]: 
                                is_mutation_present_in_sample = True
                        elif '/' in sample_genotype_value: 
                            alleles = sample_genotype_value.split('/')
                            if len(alleles) >= 2 and alleles[0] != alleles[1]:
                                is_mutation_present_in_sample = True
                        elif sample_genotype_value not in ["", ".", "0", "0/0", "0|0"] and not pd.isna(row_series[sample_col_name]):
                            is_mutation_present_in_sample = True
                        
                        if is_mutation_present_in_sample:
                            if actual_variant_type not in gene_variant_sample_counts[gene_iter_val]:
                                gene_variant_sample_counts[gene_iter_val][actual_variant_type] = set()
                            gene_variant_sample_counts[gene_iter_val][actual_variant_type].add(sample_col_name)
            
            else:  # "long" format: there is a 'sample_column' column
                for _, row_series in gene_data_for_current_gene.iterrows():
                    actual_sample_id = row_series[sample_column]
                    actual_variant_type = row_series[variant_column]

                    if pd.isna(actual_variant_type) or actual_variant_type == "Unknown" or pd.isna(actual_sample_id):
                        continue
                    
                    if actual_variant_type not in gene_variant_sample_counts[gene_iter_val]:
                        gene_variant_sample_counts[gene_iter_val][actual_variant_type] = set()
                    gene_variant_sample_counts[gene_iter_val][actual_variant_type].add(actual_sample_id)

        plot_data_list = []
        for gene, variant_dict in gene_variant_sample_counts.items():
            row_for_df = {gene_column: gene}
            has_data_for_gene = False
            for variant_type, sample_set in variant_dict.items():
                if sample_set: 
                    row_for_df[variant_type] = len(sample_set)
                    has_data_for_gene = True
            if has_data_for_gene: 
                plot_data_list.append(row_for_df)
        
        if not plot_data_list:
            ax.text(0.5, 0.5, "No data available for analysis (samples mode)", 
                      ha='center', va='center', fontsize=12)
            if set_title:
                ax.set_title(f"Top {count} Mutated Genes (Sample Prevalence)", fontsize=14, fontweight='bold')
            ax.axis('off')
            return ax

        gene_variant_counts_df = pd.DataFrame(plot_data_list).set_index(gene_column).fillna(0)
        
        gene_total_affected_samples = {}
        for gene, variant_dict in gene_variant_sample_counts.items():
            all_samples_for_gene = set()
            for _, sample_set in variant_dict.items(): 
                all_samples_for_gene.update(sample_set) 
            if all_samples_for_gene: 
                 gene_total_affected_samples[gene] = len(all_samples_for_gene)

        if not gene_total_affected_samples: 
            ax.text(0.5, 0.5, "No genes with affected samples to display (samples mode)", 
                      ha='center', va='center', fontsize=12)
            if set_title:
                ax.set_title(f"Top {count} Mutated Genes (Sample Prevalence)", fontsize=14, fontweight='bold')
            ax.axis('off')
            return ax

        gene_totals_series = pd.Series(gene_total_affected_samples).sort_values(ascending=False)
        
        total_samples_in_dataset = len(sample_cols) if samples_as_columns else data_filtered[sample_column].nunique()

        top_genes_list = gene_totals_series.index[:count].tolist()
        
        df_top_plot = gene_variant_counts_df.loc[gene_variant_counts_df.index.isin(top_genes_list)]
        df_top_plot = df_top_plot.loc[:, (df_top_plot != 0).any(axis=0)] 

        valid_top_genes_for_plot = [g for g in top_genes_list if g in df_top_plot.index]
        ordered_genes_for_plot = sorted(valid_top_genes_for_plot, key=lambda g: gene_totals_series[g])
        
        if not ordered_genes_for_plot: 
            ax.text(0.5, 0.5, "No genes in the top selection to display (samples mode)", 
                      ha='center', va='center', fontsize=12)
            if set_title:
                ax.set_title(f"Top {count} Mutated Genes (Sample Prevalence)", fontsize=14, fontweight='bold')
            ax.axis('off')
            return ax

        df_plot_final = df_top_plot.loc[ordered_genes_for_plot]

        # NORMALIZATION: Modify df_plot_final to make total sum of each row (gen) equal 
        # to total number of unique samples affected by that gen
        normalized_df_plot = pd.DataFrame(index=df_plot_final.index, columns=df_plot_final.columns)
        
        for gene in df_plot_final.index:
            # Total number of unique samples affected by this gen (for 100% of the bar)
            total_unique_samples = gene_totals_series[gene]
            
            # Actual values (unnormalized) by variant type for this gen
            current_values = df_plot_final.loc[gene]
            
            # Actual sum of values by variant type
            current_sum = current_values.sum()
            
            if current_sum > 0:  # Avoid division by zero
                # Normalization factor: how much each current unit represents relative to total unique samples
                normalization_factor = total_unique_samples / current_sum
                
                # Calculate new normalized values: current values * normalization factor
                normalized_values = current_values * normalization_factor
                
                # Assign normalized values to row corresponding to current gen
                normalized_df_plot.loc[gene] = normalized_values
        
        # Replace df_plot_final with normalized version
        df_plot_final = normalized_df_plot.fillna(0)

        variant_types_in_plot = df_plot_final.columns.tolist()
        if color_map is None:
            cmap_instance = plt.colormaps.get_cmap('tab20')
            colors_for_plot = [cmap_instance(i % cmap_instance.N) for i in range(len(variant_types_in_plot))]
        else:
            cmap_instance = plt.colormaps.get_cmap('tab20')
            colors_for_plot = [color_map.get(vt, cmap_instance(i % cmap_instance.N)) 
                               for i, vt in enumerate(variant_types_in_plot)]
        
        df_plot_final.plot(kind='barh', stacked=True, ax=ax, color=colors_for_plot, width=0.65)
        
        # Configure exactly the same aspect and positioning as variants mode
        ax.margins(x=0.05, y=0.01)
        
        for i, gene_name_in_plot in enumerate(df_plot_final.index):
            num_unique_samples_affected = gene_totals_series[gene_name_in_plot]
            percentage = (num_unique_samples_affected / total_samples_in_dataset) * 100 if total_samples_in_dataset > 0 else 0
            bar_length = df_plot_final.loc[gene_name_in_plot].sum()
            # Use same offset calculation as variants mode for consistency
            offset = 0.01 * ax.get_xlim()[1] if ax.get_xlim()[1] > 0 else 0.1
            ax.text(bar_length + offset, i, f'{percentage:.1f}%', va='center', fontsize=10)
        
        title_text = f"Top {count} Mutated Genes (Sample Prevalence)"
        if set_title:
            ax.set_title(title_text, fontsize=14, fontweight='bold') 
        ax.set_ylabel("")  # Remove Y-axis label "Genes"
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False) # Hide left axis line
        ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=True) # No ticks on Y axis, but with labels
        
        # Configure the same exact margins in both modes for total consistency
        ax.margins(x=0.05, y=0.01)

        handles, labels = ax.get_legend_handles_labels()
        legend_elements = {label: handle for label, handle in zip(labels, handles)}
        # Ensure valid_handles and valid_labels are based on df_plot_final columns that are actually in legend
        valid_handles = [legend_elements[label] for label in df_plot_final.columns if label in legend_elements]
        valid_labels = [label for label in df_plot_final.columns if label in legend_elements]


        if valid_labels: 
            num_legend_items = len(valid_labels)
            ncol_legend = min(num_legend_items, 4) 
            # More robust adjustment for vertical_offset
            base_offset = -0.20 
            row_offset_factor = 0.06
            num_legend_rows = (num_legend_items + ncol_legend - 1) // ncol_legend
            vertical_offset = base_offset - (row_offset_factor * num_legend_rows)

            ax.legend(valid_handles, valid_labels, 
                    title="Variant Classification", 
                    loc="lower center", 
                    bbox_to_anchor=(0.5, vertical_offset), 
                    ncol=ncol_legend) 
        elif ax.get_legend() is not None: 
            ax.get_legend().remove()
            
        return ax