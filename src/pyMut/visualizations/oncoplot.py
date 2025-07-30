"""
Functions for creating oncoplots (also known as waterfall plots).

This module contains the necessary functions to generate oncoplots, which are
heatmap visualizations showing mutation patterns across samples and genes.
Oncoplots are fundamental in cancer genomics for visualizing mutational landscapes.

The main function is `create_oncoplot_plot()` which can handle:
- Automatic detection of sample columns (TCGA and .GT format)
- Multiple genotype formats (A|G, A/G, etc.)
- Multi_Hit detection for samples with multiple mutation types
- Customizable color schemes for variant classifications
- Configurable parameters for top genes and maximum samples

Main functions:
- is_mutated(): Determines if a genotype represents a mutation
- detect_sample_columns(): Automatically detects sample columns in the DataFrame
- create_variant_color_mapping(): Creates color mapping for variant types
- process_mutation_matrix(): Processes mutation data into matrix format
- create_oncoplot_plot(): Main function to create the oncoplot
"""

from typing import Dict, List, Optional, Set, Tuple, TYPE_CHECKING

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

from ..utils.constants import (
    GENE_COLUMN, VARIANT_CLASSIFICATION_COLUMN, REF_COLUMN, ALT_COLUMN,
    DEFAULT_ONCOPLOT_FIGSIZE
)

if TYPE_CHECKING:
    from ..core import PyMutation

def is_mutated(genotype: str, ref: str, alt: str) -> bool:
    """
    Determina si un genotipo representa una mutación comparando con REF/ALT.
    
    Soporta múltiples formatos de genotipo:
    - Formato pipe: "A|G", "C|T" 
    - Formato slash: "A/G", "C/T"
    - Otros separadores: "A:G", "A;G"
    - Casos especiales: maneja valores faltantes y formatos no estándar
    
    Args:
        genotype: Valor del genotipo de la muestra
        ref: Alelo de referencia
        alt: Alelo alternativo
        
    Returns:
        bool: True si representa una mutación, False en caso contrario
        
    Examples:
        >>> is_mutated("A|G", "A", "G")  # True
        >>> is_mutated("A|A", "A", "G")  # False  
        >>> is_mutated("G|G", "A", "G")  # True (homocigoto alternativo)
        >>> is_mutated("./.", "A", "G")  # False
    """
    # Check for null or empty values
    if pd.isna(genotype) or pd.isna(ref) or pd.isna(alt):
        return False
        
    genotype_str = str(genotype).strip()
    ref_str = str(ref).strip()
    alt_str = str(alt).strip()
    
    # Cases that are clearly NOT mutations
    no_mutation_values = {"", ".", "0", "0/0", "0|0", "./.", ".|.", "NA", "NaN"}
    if genotype_str.upper() in no_mutation_values:
        return False
    
    # Common genotype separators
    separators = ['|', '/', ':', ';', ',']
    alleles = [genotype_str]  # Default: treat as single allele
    
    # Try to split by common separators
    for sep in separators:
        if sep in genotype_str:
            alleles = genotype_str.split(sep)
            break
    
    # Check if any allele matches the alternative
    # This detects both heterozygotes (A|G) and homozygous alternatives (G|G)
    return alt_str in alleles

def detect_sample_columns(data: pd.DataFrame) -> List[str]:
    """
    Automatically detects columns representing samples in the DataFrame.
    
    Searches for common sample naming patterns:
    - TCGA format: columns starting with "TCGA-"
    - GT format: columns ending with ".GT"
    - Other sample patterns according to standard conventions
    
    Args:
        data: DataFrame with mutation data
        
    Returns:
        List[str]: List of column names identified as samples
        
    Raises:
        ValueError: If no sample columns are detected
        
    Examples:
        >>> columns = ["Hugo_Symbol", "TCGA-AB-2988", "TCGA-AB-2869", "Variant_Classification"]
        >>> detect_sample_columns(pd.DataFrame(columns=columns))
        ['TCGA-AB-2988', 'TCGA-AB-2869']
    """
    sample_columns = []
    
    # Detect TCGA format (most common)
    tcga_cols = [col for col in data.columns if str(col).startswith("TCGA-")]
    sample_columns.extend(tcga_cols)
    
    # Detect .GT format (genotype)
    gt_cols = [col for col in data.columns if str(col).endswith(".GT")]
    sample_columns.extend(gt_cols)
    
    # Remove duplicates while maintaining order
    sample_columns = list(dict.fromkeys(sample_columns))
    
    if not sample_columns:
        raise ValueError(
            "Could not automatically detect sample columns. "
            "Please ensure columns follow standard formats like 'TCGA-*' or '*.GT', "
            "or specify columns manually."
        )
    
    return sample_columns

def create_variant_color_mapping(variants: Set[str]) -> Dict[str, np.ndarray]:
    """
    Creates a consistent color mapping for variant types.
    
    Assigns specific colors to known variant types and automatic colors
    for unrecognized variants.
    
    Args:
        variants: Set of unique variant types
        
    Returns:
        Dict[str, np.ndarray]: Dictionary mapping variant -> RGB color
        
    Examples:
        >>> variants = {"Missense_Mutation", "Nonsense_Mutation", "None"}
        >>> mapping = create_variant_color_mapping(variants)
        >>> len(mapping) == 3
        True
    """
    # Define specific colors for known variant types
    predefined_colors = {
        'Missense_Mutation': np.array([0.0, 0.8, 0.0]),        # Green
        'Nonsense_Mutation': np.array([1.0, 0.0, 0.0]),        # Red
        'Frame_Shift_Del': np.array([0.8, 0.0, 0.8]),          # Magenta
        'Frame_Shift_Ins': np.array([0.6, 0.0, 0.6]),          # Dark magenta
        'In_Frame_Del': np.array([0.0, 0.8, 0.8]),             # Cyan
        'In_Frame_Ins': np.array([0.0, 0.6, 0.6]),             # Dark cyan
        'Splice_Site': np.array([1.0, 0.5, 0.0]),              # Orange
        'Translation_Start_Site': np.array([0.8, 0.8, 0.0]),   # Yellow
        'Nonstop_Mutation': np.array([0.5, 0.0, 1.0]),         # Blue violet
        'Silent': np.array([0.7, 0.7, 0.7]),                   # Gray
        'Multi_Hit': np.array([0.0, 0.0, 0.0]),                # Black
        'None': np.array([0.95, 0.95, 0.95])                   # Very light gray
    }
    
    color_mapping = {}
    
    # Use predefined colors when available
    for variant in variants:
        if variant in predefined_colors:
            color_mapping[variant] = predefined_colors[variant]
    
    # Assign automatic colors for unrecognized variants
    unassigned_variants = [v for v in variants if v not in color_mapping]
    
    if unassigned_variants:
        # Generate automatic colors using a palette
        n_colors = len(unassigned_variants)
        auto_colors = plt.cm.tab20(np.linspace(0, 1, n_colors))
        
        for i, variant in enumerate(unassigned_variants):
            color_mapping[variant] = auto_colors[i][:3]  # RGB only, no alpha
    
    return color_mapping

def process_mutation_matrix(data: pd.DataFrame,
                           gene_column: str = GENE_COLUMN,
                           variant_column: str = VARIANT_CLASSIFICATION_COLUMN,
                           ref_column: str = REF_COLUMN,
                           alt_column: str = ALT_COLUMN,
                           sample_columns: Optional[List[str]] = None) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """
    Processes mutation data to create a genes x samples matrix.
    
    Transforms mutation data in long or wide format into an optimized matrix
    for oncoplot visualization. Automatically detects multiple mutations
    (Multi_Hit) and handles different genotype formats.
    
    Args:
        data: DataFrame with mutation data
        gene_column: Name of the gene column
        variant_column: Name of the variant classification column
        ref_column: Name of the reference allele column
        alt_column: Name of the alternative allele column
        sample_columns: Optional list of sample columns (auto-detects if None)
        
    Returns:
        Tuple[pd.DataFrame, Dict[str, int]]: 
            - Mutation matrix (genes x samples)
            - Dictionary with mutation counts per gene
            
    Raises:
        ValueError: If required columns are missing or no valid data
        
    Examples:
        >>> matrix, counts = process_mutation_matrix(mutation_data)
        >>> matrix.shape  # (genes, samples)
        (100, 50)
        >>> "Multi_Hit" in matrix.values
        True
    """
    # Validate required columns
    required_columns = [gene_column, variant_column]
    missing_columns = [col for col in required_columns if col not in data.columns]
    
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")
    
    # Detect sample columns if not provided
    if sample_columns is None:
        try:
            sample_columns = detect_sample_columns(data)
        except ValueError as e:
            # If cannot detect automatically, show warning but continue
            print(f"Warning: {e}")
            print(f"Available columns: {list(data.columns)}")
            # Try to identify columns that could be samples
            excluded_cols = {gene_column, variant_column, ref_column, alt_column, 'Tumor_Sample_Barcode'}
            potential_samples = [col for col in data.columns if col not in excluded_cols]
            if potential_samples:
                print(f"Using potential columns as samples: {potential_samples}")
                sample_columns = potential_samples
            else:
                raise ValueError("No valid sample columns detected")
        
    if not sample_columns:
        raise ValueError("No valid sample columns detected")
    
    print(f"Processing {len(sample_columns)} samples and data from {len(data)} variants...")
    
    # Create base matrix: genes x samples, initialized with 'None'
    unique_genes = data[gene_column].dropna().unique()
    unique_genes = [gene for gene in unique_genes if str(gene) != 'nan' and str(gene).strip() != '']
    
    if not unique_genes:
        raise ValueError("No valid genes found in the data")
    
    # Create DataFrame with 'None' as default value
    mutation_matrix = pd.DataFrame(
        'None',
        index=unique_genes,
        columns=sample_columns
    )
    
    # Check if REF and ALT columns exist if specified
    has_ref_alt = (ref_column in data.columns and alt_column in data.columns)
    
    # Process each mutation data row
    mutation_counts = {gene: 0 for gene in unique_genes}
    
    for _, row in data.iterrows():
        gene = row[gene_column]
        variant_type = row[variant_column]
        
        # Skip rows with missing values
        if pd.isna(gene) or pd.isna(variant_type):
            continue
            
        gene_str = str(gene).strip()
        variant_str = str(variant_type).strip()
        
        if gene_str == '' or variant_str == '' or gene_str not in unique_genes:
            continue
        
        # Get reference and alternative alleles if available
        ref_allele = row[ref_column] if has_ref_alt else None
        alt_allele = row[alt_column] if has_ref_alt else None
        
        # Process each sample for this variant
        for sample_col in sample_columns:
            if sample_col not in data.columns:
                continue
                
            genotype = row[sample_col]
            
            # Determine if there's a mutation in this sample
            has_mutation = False
            
            if has_ref_alt and not pd.isna(ref_allele) and not pd.isna(alt_allele):
                # Use REF/ALT information to determine mutation
                has_mutation = is_mutated(genotype, ref_allele, alt_allele)
            else:
                # Use heuristic based on genotype only
                # For cases without REF/ALT, consider non-empty/non-missing values as mutations
                if not pd.isna(genotype):
                    genotype_str = str(genotype).strip()
                    no_mutation_values = {"", ".", "0", "0/0", "0|0", "./.", ".|.", "NA", "NaN", "None"}
                    has_mutation = genotype_str not in no_mutation_values
            
            if has_mutation:
                current_value = mutation_matrix.loc[gene_str, sample_col]
                
                if current_value == 'None':
                    # First mutation in this cell
                    mutation_matrix.loc[gene_str, sample_col] = variant_str
                    mutation_counts[gene_str] += 1
                elif current_value != variant_str:
                    # Multiple different mutations -> Multi_Hit
                    mutation_matrix.loc[gene_str, sample_col] = 'Multi_Hit'
                # If it's the same variant, do nothing (avoid duplicates)
    
    print(f"Processed matrix: {mutation_matrix.shape[0]} genes x {mutation_matrix.shape[1]} samples")
    print(f"Genes with mutations: {sum(1 for count in mutation_counts.values() if count > 0)}")
    
    return mutation_matrix, mutation_counts

def _create_oncoplot_plot(py_mut: 'PyMutation',
                         gene_column: str = GENE_COLUMN,
                         variant_column: str = VARIANT_CLASSIFICATION_COLUMN,
                         ref_column: str = REF_COLUMN,
                         alt_column: str = ALT_COLUMN,
                         top_genes_count: int = 30,
                         max_samples: int = 180,
                         figsize: Optional[Tuple[int, int]] = None,
                         title: str = "Oncoplot",
                         ax: Optional[plt.Axes] = None) -> plt.Figure:
    """
    Creates a complete oncoplot with upper TMB panel and side gene panel.
    
    Generates a comprehensive oncoplot visualization that includes:
    - Upper panel: TMB (Tumor Mutation Burden) per sample
    - Main panel: Mutation heatmap (genes x samples)
    - Right side panel: Top mutated genes in "samples" mode
    - Legend: Variant classification types
    - Waterfall/cascade sorting algorithm for visual effect
    
    Args:
        py_mut: PyMutation object with mutation data.
        gene_column: Name of the gene column (default 'Hugo_Symbol')
        variant_column: Name of the variant classification column
        ref_column: Name of the reference allele column
        alt_column: Name of the alternative allele column
        top_genes_count: Number of top mutated genes to show (default 30)
        max_samples: Maximum number of samples to show (default 180)
        figsize: Figure size (width, height) in inches
        title: Title for the plot
        ax: Existing matplotlib axes (optional, for simple mode without panels)
        
    Returns:
        plt.Figure: matplotlib Figure object with the complete oncoplot
        
    Raises:
        ValueError: If required columns are missing or no data to visualize
        
    Examples:
        >>> fig = create_oncoplot_plot(mutation_data, top_genes_count=20)
        >>> fig.savefig('complete_oncoplot.png', dpi=300, bbox_inches='tight')
    """
    data = py_mut.data
    try:
        # Process data into mutation matrix
        mutation_matrix, mutation_counts = process_mutation_matrix(
            data, gene_column, variant_column, ref_column, alt_column
        )
    
        if mutation_matrix.empty:
            raise ValueError("No mutation data to visualize")
        
        # CALCULATE TMB USING ALL GENES (BEFORE FILTERING)
        # This must match the calculation in variants_per_sample
        all_genes_tmb = (mutation_matrix != 'None').sum(axis=0)
        median_tmb_all = np.median(all_genes_tmb)
        
        # Select top mutated genes
        top_genes = sorted(mutation_counts.items(), key=lambda x: x[1], reverse=True)
        top_genes = [gene for gene, count in top_genes[:top_genes_count] if count > 0]
        
        if not top_genes:
            raise ValueError("No genes with mutations found")
        
        # Filter matrix to top genes
        plot_matrix = mutation_matrix.loc[top_genes].copy()
        
        # ========================================================================
        # STANDARD WATERFALL/CASCADE ALGORITHM (MAFTOOLS/COMPLEXHEATMAP)
        # ========================================================================
        # This algorithm reproduces the standard behavior of oncoplots:
        # 1. Sort genes by mutation frequency (most mutated at top)
        # 2. Sort samples lexicographically by mutation pattern
        # 3. Creates the characteristic cascade visual effect
        
        print("Applying standard cascade algorithm (maftools)...")
        
        # STEP 1: Convert categorical matrix to numeric for the algorithm
        unique_variants = set()
        for row in plot_matrix.values.flatten():
            unique_variants.add(row)
        
        unique_values = sorted(unique_variants)
        value_to_num = {value: i for i, value in enumerate(unique_values)}
        
        # Convert to numeric matrix (0 = 'None', >0 = mutation)
        numeric_matrix = plot_matrix.replace(value_to_num)
        waterfall_matrix = numeric_matrix.copy()
        
        # STEP 2: Sort genes by frequency (as in maftools)
        # Count zeros per gene (samples NOT mutated)
        # Check that 'None' exists in the mapping
        if 'None' in value_to_num:
            gene_zero_counts = (waterfall_matrix == value_to_num['None']).sum(axis=1)
        else:
            # If no 'None', use minimum value (representing non-mutation)
            min_value = min(value_to_num.values())
            gene_zero_counts = (waterfall_matrix == min_value).sum(axis=1)
        
        # Sort genes: fewer zeros first (most frequent at top)
        sorted_genes = gene_zero_counts.sort_values(ascending=True).index.tolist()
        waterfall_matrix = waterfall_matrix.loc[sorted_genes]
        
        # STEP 3: Implement the cascade/waterfall algorithm as in maftools
        # This is the standard algorithm used by maftools, ComplexHeatmap, etc.
        # Sorts samples lexicographically by mutation pattern
        
        # Create binary matrix for sorting (0 = not mutated, 1 = mutated)
        binary_matrix = plot_matrix.copy()
        for gene in sorted_genes:
            for sample in plot_matrix.columns:
                if plot_matrix.loc[gene, sample] == 'None':
                    binary_matrix.loc[gene, sample] = 0
                else:
                    binary_matrix.loc[gene, sample] = 1
        
        # STANDARD ALGORITHM: Sort samples using all genes as criteria
        # This reproduces the behavior of maftools and ComplexHeatmap
        # The sorting is lexicographic: first by most frequent gene,
        # then by second most frequent gene, etc.
        sorted_samples = binary_matrix.T.sort_values(
            by=sorted_genes,  # Use genes already sorted by frequency
            ascending=False   # Descending so samples with more mutations appear first
        ).index.tolist()
        
        # STEP 4: Limit number of samples if necessary
        if max_samples and len(sorted_samples) > max_samples:
            waterfall_matrix = waterfall_matrix.iloc[:, :max_samples]
            sorted_samples = sorted_samples[:max_samples]
        
        # STEP 5: Apply cascade ordering to original matrix for plotting
        plot_matrix = plot_matrix.loc[sorted_genes, sorted_samples]
        
        # Show cascade algorithm information
        print("Cascade applied:")
        print(f"  - Genes sorted by frequency: {len(sorted_genes)}")
        print(f"  - Samples sorted by cascade algorithm: {len(sorted_samples)}")
        print(f"  - Final samples shown: {plot_matrix.shape[1]}")
        
        # 6. Get TMB for visualization
        # TMB for heatmap: sorted according to cascade effect
        sorted_tmb_heatmap = [all_genes_tmb[sample] for sample in sorted_samples]
        median_tmb = np.median(sorted_tmb_heatmap)
        
        # TMB for upper panel: MUST use same cascade order as heatmap
        # This ensures TMB shows the characteristic descending pattern
        tmb_display_order = sorted_samples  # Same order as heatmap
        sorted_tmb_display = sorted_tmb_heatmap  # Already in correct order
        
        # Create color mapping
        unique_variants = set()
        for row in plot_matrix.values.flatten():
            unique_variants.add(row)
        
        color_mapping = create_variant_color_mapping(unique_variants)
        
        # Convert categorical values to numeric for heatmap
        unique_values = sorted(unique_variants)
        value_to_num = {value: i for i, value in enumerate(unique_values)}
        
        numeric_matrix = plot_matrix.applymap(lambda x: value_to_num[x])
        
        # Configure figure and subplots
        if figsize is None:
            figsize = DEFAULT_ONCOPLOT_FIGSIZE
            
        if ax is None:
            # Create complete figure with grid: TMB above, main oncoplot, genes on right side
            fig = plt.figure(figsize=figsize)
            
            # Define grid layout: 3 rows x 2 columns
            # height_ratios: [2.5, 10, 1.5] = taller TMB above, large heatmap, legend with space below
            # width_ratios: [10, 1] = wide heatmap, narrow gene panel
            gs = fig.add_gridspec(3, 2, 
                                height_ratios=[2.5, 10, 1.5], 
                                width_ratios=[10, 1],
                                hspace=0.05, 
                                wspace=0.02)
            
            # Create subplots according to grid
            ax_tmb = fig.add_subplot(gs[0, 0])      # TMB panel (top left)
            ax_main = fig.add_subplot(gs[1, 0])     # Main heatmap panel (center left)
            ax_genes = fig.add_subplot(gs[1, 1])    # Gene panel (center right)
            ax_legend = fig.add_subplot(gs[2, :])   # Legend panel (bottom, full span)
            
            # === TMB PANEL (UPPER) ===
            # Calculate TMB using same logic as variants_per_sample_plot
            # Use same sample columns detected for heatmap
            detected_sample_columns = mutation_matrix.columns.tolist()
            
            # Use same logic as variants_per_sample_plot
            variant_counts_tmb = {}
            
            # Group by Variant_Classification (same as in summary.py)
            for variant_type in data[variant_column].unique():
                variant_subset = data[data[variant_column] == variant_type]
                
                for sample_col in detected_sample_columns:
                    # Count mutations using same logic as summary.py
                    sample_variants = variant_subset[sample_col].apply(
                        lambda x: 1 if '|' in str(x) and str(x).split('|')[0] != str(x).split('|')[1] else 0
                    ).sum()
                    
                    if sample_col not in variant_counts_tmb:
                        variant_counts_tmb[sample_col] = {}
                    
                    if sample_variants > 0:
                        variant_counts_tmb[sample_col][variant_type] = sample_variants
            
            # Create DataFrame for TMB (same as in summary.py)
            samples_df_tmb = []
            for sample, variants in variant_counts_tmb.items():
                for var_type, count in variants.items():
                    samples_df_tmb.append({'Sample': sample, 'Variant_Classification': var_type, 'Count': count})
            
            if samples_df_tmb:
                processed_df_tmb = pd.DataFrame(samples_df_tmb)
                tmb_df = processed_df_tmb.pivot(index='Sample', columns='Variant_Classification', values='Count').fillna(0)
                
                # Filter to selected samples and order as in heatmap
                tmb_df = tmb_df.loc[tmb_df.index.intersection(sorted_samples)]
                tmb_df = tmb_df.reindex(sorted_samples, fill_value=0)
                
                # Prepare colors for TMB (same mapping as heatmap)
                tmb_colors = [color_mapping.get(vt, [0.7, 0.7, 0.7]) for vt in tmb_df.columns]
                
                # Create stacked bars for TMB
                tmb_df.plot(kind='bar', stacked=True, ax=ax_tmb, color=tmb_colors, width=0.8, legend=False)
                
                # Configure TMB panel
                ax_tmb.set_ylabel('TMB', fontsize=10)
                ax_tmb.set_xlim(-0.5, len(sorted_samples) - 0.5)
                
                # Y limits adjust automatically based on data
                max_tmb = tmb_df.sum(axis=1).max() if not tmb_df.empty else 10
                ax_tmb.set_ylim(0, max_tmb * 1.1)
                
                ax_tmb.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
                # REMOVE "Sample" label from X axis of TMB panel
                ax_tmb.set_xlabel('')
                ax_tmb.spines['top'].set_visible(False)
                ax_tmb.spines['right'].set_visible(False)
                ax_tmb.spines['bottom'].set_visible(False)
                ax_tmb.grid(axis='y', alpha=0.3)
            else:
                # If no TMB data, show message
                ax_tmb.text(0.5, 0.5, 'TMB not available', ha='center', va='center', fontsize=10)
                ax_tmb.set_xlim(0, 1)
                ax_tmb.set_ylim(0, 1)
                ax_tmb.axis('off')
            
            # === GENE PANEL (RIGHT SIDE) ===
            # Create custom side panel with stacked bars and synchronized colors
            try:
                # Calculate variant counts per gene using processed matrix
                gene_variant_counts = {}
                
                # For each gene, count how many samples have each variant type
                for gene in sorted_genes:  # Use cascade order
                    gene_variant_counts[gene] = {}
                    gene_row = plot_matrix.loc[gene]
                    
                    # Count each variant type (excluding 'None')
                    for variant_type in gene_row.value_counts().index:
                        if variant_type != 'None':
                            count = gene_row.value_counts()[variant_type]
                            gene_variant_counts[gene][variant_type] = count
                
                # Create DataFrame for stacked bars
                variant_types = set()
                for gene_variants in gene_variant_counts.values():
                    variant_types.update(gene_variants.keys())
                
                variant_types = sorted(list(variant_types))
                
                # Build DataFrame with genes as index and variant types as columns
                stacked_data = []
                for gene in sorted_genes:  # Use cascade order
                    row_data = {'gene': gene}
                    for variant_type in variant_types:
                        row_data[variant_type] = gene_variant_counts[gene].get(variant_type, 0)
                    stacked_data.append(row_data)
                
                df_stacked = pd.DataFrame(stacked_data).set_index('gene')
                
                # REVERSE order to match heatmap (top to bottom)
                # Heatmap shows genes in sorted_genes order (top to bottom)
                # Horizontal bars plot from bottom to top by default
                # Therefore, we reverse the DataFrame to match visually
                df_stacked = df_stacked.iloc[::-1]  # Reverse row order
                
                # Ensure numeric columns are float to avoid FutureWarning
                for col in df_stacked.columns:
                    if pd.api.types.is_numeric_dtype(df_stacked[col]):
                        df_stacked[col] = df_stacked[col].astype(float)
                
                # Calculate totals and percentages for each gene
                gene_totals = df_stacked.sum(axis=1)
                total_samples = len(plot_matrix.columns)
                
                # Normalize to show sample prevalence (as percentage of total samples)
                df_percentage = df_stacked.copy()
                for gene in df_percentage.index:
                    total_affected = gene_totals[gene]
                    percentage = (total_affected / total_samples) * 100
                    # Scale each variant proportionally so sum equals total percentage
                    if total_affected > 0: # Avoid division by zero
                        scaling_factor = percentage / total_affected
                        # Apply scaling directly to the row
                        df_percentage.loc[gene] = df_percentage.loc[gene] * scaling_factor
                
                # Create color mapping synchronized with main heatmap
                stacked_colors = []
                for variant_type in variant_types:
                    if variant_type in color_mapping:
                        stacked_colors.append(color_mapping[variant_type])
                    else:
                        # Default color if not in mapping
                        stacked_colors.append([0.7, 0.7, 0.7])
                
                # Create horizontal stacked bars
                df_percentage.plot(kind='barh', stacked=True, ax=ax_genes, 
                                 color=stacked_colors, width=0.65)
                
                # Configure gene panel - NO gene labels and NO X label
                ax_genes.set_xlabel('')  # REMOVE "Samples (%)" label from X axis
                ax_genes.set_ylabel('')  # No Y label
                
                # Calculate X limit based on data
                max_percentage = df_percentage.sum(axis=1).max() if not df_percentage.empty else 50
                ax_genes.set_xlim(0, max_percentage * 1.1)
                ax_genes.set_ylim(-0.5, len(sorted_genes) - 0.5)
                
                # REMOVE gene labels from Y axis
                ax_genes.set_yticklabels([''] * len(sorted_genes))  # Empty labels
                
                # FIX PERCENTAGES: use exact same logic as summary.py
                # Detect sample columns same as in summary.py
                sample_cols_for_percentage = [col for col in data.columns if str(col).startswith("TCGA-")]
                total_samples_in_dataset = len(sample_cols_for_percentage)  # Same as in summary.py
                
                for i, gene in enumerate(df_percentage.index):
                    # Reverse the inverted order to get correct gene
                    gene_in_original_order = sorted_genes[len(sorted_genes) - 1 - i]  
                    
                    # USE SAME LOGIC AS SUMMARY.PY:
                    # Count unique samples affected by this gene using original data
                    gene_data_filtered = data[data[gene_column] == gene_in_original_order]
                    
                    # Count unique affected samples (same as in summary.py)
                    affected_samples_set = set()
                    
                    for sample_col_name in sample_cols_for_percentage:
                        for _, row_series in gene_data_filtered.iterrows():
                            sample_genotype_value = str(row_series[sample_col_name]).strip().upper()
                            
                            is_mutation_present = False
                            if '|' in sample_genotype_value:
                                alleles = sample_genotype_value.split('|')
                                if len(alleles) >= 2 and alleles[0] != alleles[1]: 
                                    is_mutation_present = True
                            elif '/' in sample_genotype_value: 
                                alleles = sample_genotype_value.split('/')
                                if len(alleles) >= 2 and alleles[0] != alleles[1]:
                                    is_mutation_present = True
                            elif sample_genotype_value not in ["", ".", "0", "0/0", "0|0"] and not pd.isna(row_series[sample_col_name]):
                                is_mutation_present = True
                            
                            if is_mutation_present:
                                affected_samples_set.add(sample_col_name)
                    
                    # Calculate percentage exactly same as summary.py
                    num_unique_samples_affected = len(affected_samples_set)
                    real_percentage = (num_unique_samples_affected / total_samples_in_dataset) * 100 if total_samples_in_dataset > 0 else 0
                    
                    if real_percentage > 0:
                        # Calculate offset for text
                        bar_length = df_percentage.loc[gene].sum()
                        offset = max_percentage * 0.02
                        ax_genes.text(bar_length + offset, i, f'{real_percentage:.1f}%', 
                                    va='center', fontsize=9)
                
                # Remove legend from side panel (shown below)
                legend = ax_genes.get_legend()
                if legend is not None:
                    legend.remove()
                
                # Configure panel style
                ax_genes.spines['top'].set_visible(False)
                ax_genes.spines['right'].set_visible(False)
                ax_genes.spines['left'].set_visible(False)
                ax_genes.tick_params(axis='y', which='both', left=False, right=False)
                ax_genes.tick_params(axis='x', labelsize=9)
                ax_genes.grid(axis='x', alpha=0.3)
                ax_genes.margins(x=0.05, y=0.01)
                
            except Exception as e:
                print(f"Warning: Error creating custom gene panel: {e}")
                import traceback
                traceback.print_exc()
                # If fails, create empty panel
                ax_genes.text(0.5, 0.5, 'Gene panel\nnot available', 
                             ha='center', va='center', fontsize=10)
                ax_genes.set_xlim(0, 1)
                ax_genes.set_ylim(0, 1)
                ax_genes.axis('off')
            
            # === LEGEND PANEL (BOTTOM) ===
            # Create legend elements for variant types
            legend_elements = []
            for variant in sorted(unique_variants):
                if variant != 'None':  # Don't show 'None' in legend
                    color = color_mapping[variant]
                    label = variant.replace('_', ' ')
                    legend_elements.append(
                        plt.Rectangle((0, 0), 1, 1, facecolor=color, label=label)
                    )
            
            if legend_elements:
                ax_legend.legend(
                    handles=legend_elements,
                    title='Variant Classification',
                    loc='center',
                    ncol=min(len(legend_elements), 6),
                    fontsize=10,
                    title_fontsize=11,
                    frameon=False
                )
            
            # Add sample count centered above legend
            ax_legend.text(0.5, 0.85, f'Samples (n={len(sorted_samples)})', 
                          ha='center', va='center', transform=ax_legend.transAxes,
                          fontsize=11, fontweight='bold')
            
            ax_legend.axis('off')
            
            # Main figure title
            fig.suptitle(title, fontsize=16, fontweight='bold', y=0.95)
            
        else:
            # Simple mode: just heatmap without additional panels
            fig = ax.get_figure()
            ax_main = ax
        
        # === MAIN HEATMAP ===
        # Create custom colormap
        colors = [color_mapping[value] for value in unique_values]
        custom_cmap = mcolors.ListedColormap(colors)
        
        # Create main heatmap
        im = ax_main.imshow(
            numeric_matrix.values,
            cmap=custom_cmap, 
            aspect='auto',
            interpolation='nearest'
        )
        
        # Configure main heatmap axes
        if ax is not None:  # Only set title if no TMB panel
            ax_main.set_title(title, fontsize=16, fontweight='bold', pad=20)
        
        # REMOVE X axis label to avoid overlap with TMB
        # ax_main.set_xlabel(f'Samples (n={len(sorted_samples)})', fontsize=12)  # Commented
        ax_main.set_ylabel('Genes', fontsize=12)
        
        # Configure ticks and labels
        ax_main.set_xticks(range(len(sorted_samples)))
        ax_main.set_yticks(range(len(sorted_genes)))
        ax_main.set_yticklabels(sorted_genes, fontsize=10)
        
        # COMPLETELY REMOVE X axis marks and labels from waterfall plot
        ax_main.set_xticklabels([''] * len(sorted_samples))  # Empty labels
        ax_main.tick_params(axis='x', which='both', bottom=False, top=False)  # No X marks
        
        # Add subtle grid
        ax_main.set_xticks(np.arange(-0.5, len(sorted_samples), 1), minor=True)
        ax_main.set_yticks(np.arange(-0.5, len(sorted_genes), 1), minor=True)
        ax_main.grid(which='minor', color='white', linestyle='-', linewidth=0.5)
        
        # If simple mode, add legend to heatmap
        if ax is not None:
            legend_elements = []
            for variant in sorted(unique_variants):
                if variant != 'None':
                    color = color_mapping[variant]
                    label = variant.replace('_', ' ')
                    legend_elements.append(
                        plt.Rectangle((0, 0), 1, 1, facecolor=color, label=label)
                    )
    
            if legend_elements:
                ax_main.legend(
                    handles=legend_elements, 
                    title='Variant Classification',
                    bbox_to_anchor=(1.05, 1), 
                    loc='upper left',
                    fontsize=9,
                    title_fontsize=10
                )
    
        print("Oncoplot created successfully:")
        print(f"  - {plot_matrix.shape[0]} genes")
        print(f"  - {plot_matrix.shape[1]} samples")
        print(f"  - {len(color_mapping)} variant types")
        print("  - Standard cascade algorithm applied (maftools)")
        
        return fig 
        
    except Exception as e:
        print(f"Error creating oncoplot: {e}")
        raise 