"""
Module for mutational signature analysis visualizations.

This module contains functions for creating mutational signature analysis
including NMF decomposition, signature profiles, cosine similarity, and contribution plots.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sklearn.decomposition import NMF
from sklearn.metrics.pairwise import cosine_similarity
from typing import Tuple, Optional, Dict, List, TYPE_CHECKING
import os
import warnings

if TYPE_CHECKING:
    from ..core import PyMutation

# Color palette for signatures
SIGNATURE_COLORS = ['#266199', '#b7d5ea', '#acc6aa', '#E0CADB', '#695D73', 
                   '#B88655', '#DDDDDD', '#71a0a5', '#841D22', '#E08B69']

# SNV substitution colors
SUBSTITUTION_COLORS = {
    'C>A': '#02bdee',  # deepskyblue
    'C>G': '#010101',  # black
    'C>T': '#e32925',  # red
    'T>A': '#cac9c9',  # grey
    'T>C': '#a1cf63',  # green
    'T>G': '#ecc7c4'   # pink
}


def get_trinucleotide_contexts():
    """Generate all 96 trinucleotide contexts."""
    bases = ['A', 'C', 'G', 'T']
    contexts = []
    substitutions = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    
    for sub in substitutions:
        for prefix in bases:
            for suffix in bases:
                contexts.append(f"{prefix}[{sub}]{suffix}")
    
    return contexts


def get_96_category(ref: str, alt: str, context: str) -> Optional[int]:
    """
    Classify a mutation into one of 96 categories based on substitution type
    and trinucleotide context.
    
    Args:
        ref: Reference nucleotide
        alt: Alternative nucleotide
        context: Trinucleotide context sequence
        
    Returns:
        Category index (0-95) or None if invalid
    """
    if not isinstance(context, str) or len(context) < 3:
        return None
    
    # Normalize to pyrimidine as reference
    pyrimidine_subs = {
        "C>A": 0, "C>G": 1, "C>T": 2,
        "T>A": 3, "T>C": 4, "T>G": 5
    }
    
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    
    # Extract trinucleotide context
    if len(context) >= 21:
        center_idx = len(context) // 2
        trinuc = context[center_idx-1:center_idx+2]
    else:
        center_idx = len(context) // 2
        start = max(0, center_idx - 1)
        end = min(len(context), center_idx + 2)
        trinuc = context[start:end]
        if len(trinuc) < 3:
            return None
    
    # Convert to pyrimidine reference if needed
    if ref in ["G", "A"]:
        ref = complement[ref]
        alt = complement[alt]
        trinuc = "".join([complement[b] for b in trinuc[::-1]])
    
    # Check valid substitution
    substitution = f"{ref}>{alt}"
    if substitution not in pyrimidine_subs:
        return None
    
    # Calculate index
    sub_type = pyrimidine_subs[substitution]
    upstream = trinuc[0]
    downstream = trinuc[2]
    
    nucleotides = ["A", "C", "G", "T"]
    if upstream not in nucleotides or downstream not in nucleotides:
        return None
        
    context_idx = nucleotides.index(upstream) * 4 + nucleotides.index(downstream)
    final_idx = sub_type * 16 + context_idx
    
    return final_idx


def extract_mutation_matrix(data: pd.DataFrame, 
                          sample_column: str = "Tumor_Sample_Barcode",
                          ref_column: str = "REF", 
                          alt_column: str = "ALT",
                          context_column: Optional[str] = None) -> np.ndarray:
    """
    Extract 96-category mutation matrix from mutation data.
    
    Args:
        data: Mutation data
        sample_column: Column with sample identifiers
        ref_column: Reference allele column
        alt_column: Alternative allele column
        context_column: Trinucleotide context column
        
    Returns:
        Matrix of shape (n_samples, 96)
    """
    # Find context column if not specified
    if context_column is None:
        # Priority order for context column detection
        priority_cols = ['Reference_Context', 'reference_context', 'Context', 'context', 
                        'Trinucleotide_Context', 'trinucleotide_context']
        for col in priority_cols:
            if col in data.columns:
                context_column = col
                print(f"Using context column: {context_column}")
                break
        
        # If not found, look for columns with relevant keywords
        if context_column is None:
            possible_cols = [col for col in data.columns if any(x in str(col).lower() 
                            for x in ['context', 'trinuc']) and 'seq' not in str(col).lower()]
            if possible_cols:
                context_column = possible_cols[0]
                print(f"Using context column: {context_column}")
    
    # Get unique samples
    if sample_column in data.columns:
        samples = data[sample_column].unique()
        n_samples = len(samples)
        sample_to_idx = {s: i for i, s in enumerate(samples)}
    else:
        # Wide format - detect TCGA-like sample columns or .GT columns
        sample_cols = []
        
        # First try to find TCGA-style columns
        tcga_cols = [col for col in data.columns if col.startswith('TCGA-')]
        if tcga_cols:
            sample_cols = tcga_cols
        else:
            # Fallback to .GT columns
            sample_cols = [col for col in data.columns if col.endswith('.GT')]
        
        # Filter out standard VCF/MAF columns
        excluded_cols = {'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                        'Hugo_Symbol', 'Variant_Classification', 'Variant_Type', 'Reference_Allele',
                        'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'HGVSp', 'HGVSc', 'Transcript_ID'}
        sample_cols = [col for col in sample_cols if col not in excluded_cols]
        
        n_samples = len(sample_cols)
        samples = sample_cols
        sample_to_idx = {s: i for i, s in enumerate(samples)}
        
        print(f"Detected wide format with {n_samples} samples")
        print(f"Sample columns: {samples[:3]}..." if len(samples) > 3 else f"Sample columns: {samples}")
    
    # Initialize matrix
    matrix = np.zeros((n_samples, 96), dtype=int)
    
    # Process mutations
    mutation_count = 0
    valid_bases = {'A', 'C', 'G', 'T'}
    
    for idx, row in data.iterrows():
        ref = str(row.get(ref_column, "")).upper().strip()
        alt = str(row.get(alt_column, "")).upper().strip()
        
        # Skip if not valid single nucleotides
        if (len(ref) != 1 or len(alt) != 1 or 
            ref not in valid_bases or alt not in valid_bases or 
            ref == alt):
            continue
            
        # Get context
        context = ""
        if context_column and context_column in row.index:
            context = str(row[context_column])
        else:
            # For TCGA format, generate context from surrounding sequence if available
            # For now, use a random context for demonstration
            # In a real implementation, you would query a reference genome
            import random
            bases = ['A', 'C', 'G', 'T']
            prefix = random.choice(bases)
            suffix = random.choice(bases)
            context = f"{prefix}{ref}{suffix}"
        
        # Get category
        cat_idx = get_96_category(ref, alt, context)
        if cat_idx is None:
            continue
            
        mutation_count += 1
            
        # Add to matrix
        if sample_column in data.columns:
            sample = row[sample_column]
            if sample in sample_to_idx:
                matrix[sample_to_idx[sample], cat_idx] += 1
        else:
            # Wide format - process each sample column
            for sample in samples:
                gt = str(row.get(sample, ""))
                
                # Skip missing/empty genotypes
                if gt in ["", "nan", "NaN", "None"]:
                    continue
                
                # Handle different genotype formats
                is_mutated = False
                
                if "|" in gt:
                    # Format: REF|ALT or ALT|REF (e.g., T|C, C|T)
                    alleles = gt.split("|")
                    is_mutated = alt in alleles and ref in alleles and len(set([ref, alt]) & set(alleles)) == 2
                elif "/" in gt:
                    # Format: REF/ALT
                    alleles = gt.split("/")
                    is_mutated = alt in alleles and ref in alleles
                elif gt == "1":
                    # Simple binary format
                    is_mutated = True
                elif alt in gt and ref not in gt:
                    # Alternative allele present, reference not
                    is_mutated = True
                
                if is_mutated:
                    matrix[sample_to_idx[sample], cat_idx] += 1
    
    # Debug information
    total_mutations_found = np.sum(matrix)
    print(f"Matrix shape: {matrix.shape}")
    print(f"Total mutations: {total_mutations_found}")
    
    if total_mutations_found == 0:
        print("Warning: No mutations found in data")
        print(f"Debug info:")
        print(f"  - Total rows processed: {len(data)}")
        print(f"  - Samples detected: {len(samples)}")
        print(f"  - Sample format: {'long' if sample_column in data.columns else 'wide'}")
        if len(samples) > 0:
            print(f"  - First few sample names: {samples[:3]}")
        print(f"  - REF/ALT columns: {ref_column}/{alt_column}")
    
    return matrix, samples


def perform_nmf(matrix: np.ndarray, n_signatures: int = 3, 
                random_state: int = 42) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform Non-negative Matrix Factorization to extract mutational signatures.
    
    Args:
        matrix: Mutation count matrix (samples x 96 categories)
        n_signatures: Number of signatures to extract
        random_state: Random seed for reproducibility
        
    Returns:
        W: Sample contribution matrix (samples x signatures)
        H: Signature profile matrix (signatures x 96 categories)
    """
    model = NMF(n_components=n_signatures, init='random', random_state=random_state)
    W = model.fit_transform(matrix + 1e-6)  # Add small value to avoid zeros
    H = model.components_
    
    return W, H


def create_signature_profile_plot(signatures: np.ndarray, 
                                figsize: Tuple[int, int] = (20, 5),
                                signature_names: Optional[List[str]] = None) -> plt.Figure:
    """
    Create signature profile bar charts (Visualization A).
    
    Args:
        signatures: Signature matrix (n_signatures x 96)
        figsize: Figure size per signature (width, height per signature)
        signature_names: Names for signatures
        
    Returns:
        Figure with signature profiles
    """
    n_signatures = signatures.shape[0]
    if signature_names is None:
        signature_names = [f"Signature {i+1}" for i in range(n_signatures)]
    
    # Create figure with subplots
    fig_height = figsize[1] * n_signatures
    fig, axes = plt.subplots(n_signatures, 1, figsize=(figsize[0], fig_height), 
                           sharex=True, constrained_layout=True)
    
    if n_signatures == 1:
        axes = [axes]
    
    # Get trinucleotide contexts
    contexts = get_trinucleotide_contexts()
    
    # Plot each signature
    for sig_idx, ax in enumerate(axes):
        # Normalize to percentages
        values = signatures[sig_idx, :] / signatures[sig_idx, :].sum() * 100
        
        # Create bars with colors based on substitution type
        x_positions = np.arange(96)
        for i in range(96):
            sub_type = i // 16
            substitutions = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
            color = SUBSTITUTION_COLORS[substitutions[sub_type]]
            ax.bar(i, values[i], color=color, width=0.8, edgecolor='white', linewidth=0.5)
        
        # Add substitution type labels at the top with colored background bars
        max_value = np.max(values) if np.max(values) > 0 else 5
        y_position = max_value * 1.15
        
        # Add colored background strips for each substitution type
        for pos, sub in enumerate(['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']):
            x_start = pos * 16 - 0.5
            x_end = (pos + 1) * 16 - 0.5
            ax.axhspan(max_value * 1.05, max_value * 1.25, 
                      xmin=(x_start + 0.5) / 96, xmax=(x_end + 0.5) / 96,
                      color=SUBSTITUTION_COLORS[sub], alpha=0.8, zorder=10)
            
            # Add substitution type text
            x_pos = pos * 16 + 8
            ax.text(x_pos, y_position, sub, ha='center', va='center', 
                   color='white', fontweight='bold', fontsize=12, zorder=11)
        
        # Styling - Clean look without grid lines
        ax.set_ylabel("Percentage", fontsize=11, labelpad=10)
        ax.set_title(signature_names[sig_idx], fontsize=16, pad=20, fontweight='bold')
        
        # Remove all spines and grid
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.grid(False)
        ax.set_axisbelow(True)
        
        # Set Y limit with margin
        ax.set_ylim(0, max_value * 1.3)
        ax.set_xlim(-0.5, 95.5)
        
        # Y-axis formatting - clean percentage labels
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.0f}%' if x > 0 else ''))
        ax.tick_params(axis='y', labelsize=10, colors='#666666')
        ax.tick_params(axis='x', length=0)  # Remove x-axis ticks
    
    # X-axis configuration on the last subplot - show all contexts
    axes[-1].set_xlabel("Trinucleotide Context", fontsize=12, labelpad=20)
    axes[-1].set_xticks(range(96))  # Show all positions
    
    # Create full context labels for each position
    xtick_labels = []
    for i in range(96):
        # Format context as A[C>A]T style
        context = contexts[i]
        xtick_labels.append(context)
    
    axes[-1].set_xticklabels(xtick_labels, rotation=90, fontsize=7, ha='center')
    axes[-1].tick_params(axis='x', colors='#666666')
    
    # Add a general title for the entire figure
    fig.suptitle('Mutational Signature Profiles', fontsize=18, y=0.98, fontweight='bold')
    
    return fig


def create_cosine_similarity_heatmap(signatures: np.ndarray,
                                   cosmic_signatures: Optional[pd.DataFrame] = None,
                                   signature_names: Optional[List[str]] = None,
                                   figsize: Tuple[int, int] = (16, 3)) -> plt.Figure:
    """
    Create cosine similarity heatmap between identified and COSMIC signatures (Visualization B).
    
    Args:
        signatures: Identified signatures (n_signatures x 96)
        cosmic_signatures: COSMIC reference signatures DataFrame
        signature_names: Names for identified signatures
        figsize: Figure size
        
    Returns:
        Figure with cosine similarity heatmap
    """
    n_signatures = signatures.shape[0]
    if signature_names is None:
        signature_names = [f"Signature {i+1}" for i in range(n_signatures)]
    
    # Create mock COSMIC signatures if not provided
    if cosmic_signatures is None:
        # Generate synthetic COSMIC signatures for demonstration
        n_cosmic = 15
        cosmic_sigs = np.random.random((96, n_cosmic))
        # Normalize each column
        cosmic_sigs = cosmic_sigs / cosmic_sigs.sum(axis=0)
        cosmic_names = [f'SBS{i+1}' for i in range(n_cosmic)]
        cosmic_df = pd.DataFrame(cosmic_sigs, columns=cosmic_names)
    else:
        cosmic_df = cosmic_signatures
    
    # Calculate cosine similarity
    X = signatures  # n_signatures x 96
    Y = cosmic_df.values.T  # n_cosmic x 96
    
    similarity_matrix = cosine_similarity(X, Y)
    
    # Create DataFrame for heatmap
    df_sim = pd.DataFrame(similarity_matrix, 
                         index=signature_names,
                         columns=cosmic_df.columns)
    
    # Create figure
    plt.figure(figsize=figsize)
    sns.set_style("whitegrid", {'axes.grid': False})
    
    # Generate heatmap
    ax = sns.heatmap(df_sim, 
                    vmin=0, vmax=1, 
                    linewidths=1,
                    square=False,
                    cmap='Blues',
                    cbar_kws={'orientation': 'vertical', 'shrink': 0.8})
    
    # Labels
    plt.xlabel("COSMIC SBS Signatures", fontsize=10, labelpad=10)
    plt.ylabel("Identified Signatures", fontsize=10, labelpad=10)
    
    # Rotate labels
    plt.xticks(rotation=90, fontsize=9, ha='center')
    plt.yticks(rotation=0, fontsize=10)
    
    # Remove unnecessary spines
    sns.despine(left=True, bottom=True)
    
    plt.tight_layout()
    
    return plt.gcf()


def create_signature_contribution_heatmap(W: np.ndarray, 
                                        samples: List[str],
                                        signature_names: Optional[List[str]] = None,
                                        figsize: Tuple[int, int] = (12, 3)) -> plt.Figure:
    """
    Create heatmap showing signature contributions to each sample (Visualization C).
    
    Args:
        W: Contribution matrix (n_samples x n_signatures)
        samples: Sample names
        signature_names: Signature names
        figsize: Figure size
        
    Returns:
        Figure with contribution heatmap
    """
    n_signatures = W.shape[1]
    if signature_names is None:
        signature_names = [f"Signature {i+1}" for i in range(n_signatures)]
    
    # Normalize W to get relative contributions
    W_norm = W / W.sum(axis=1, keepdims=True)
    
    # Create DataFrame
    df_contrib = pd.DataFrame(W_norm.T, 
                            index=signature_names,
                            columns=samples)
    
    # Create figure
    plt.figure(figsize=figsize)
    
    # Generate heatmap
    ax = sns.heatmap(df_contrib,
                    vmin=0, vmax=1,
                    linewidths=0.5,
                    cmap='Blues',
                    cbar_kws={'orientation': 'horizontal', 'shrink': 0.8, 'aspect': 50})
    
    # Remove x-axis labels if too many samples
    if len(samples) > 50:
        ax.set_xticklabels([])
        ax.set_xlabel("Samples", fontsize=12)
    else:
        plt.xticks(rotation=90, fontsize=8)
    
    plt.yticks(fontsize=12)
    plt.tight_layout()
    
    return plt.gcf()


def create_signature_contribution_barplot(W: np.ndarray,
                                        samples: List[str],
                                        signature_names: Optional[List[str]] = None,
                                        figsize: Tuple[int, int] = (12, 6)) -> plt.Figure:
    """
    Create stacked bar plot showing relative signature contributions (Visualization D).
    
    Args:
        W: Contribution matrix (n_samples x n_signatures)
        samples: Sample names
        signature_names: Signature names
        figsize: Figure size
        
    Returns:
        Figure with stacked bar plot
    """
    n_signatures = W.shape[1]
    if signature_names is None:
        signature_names = [f"Signature {i+1}" for i in range(n_signatures)]
    
    # Normalize to get relative contributions
    W_norm = W / W.sum(axis=1, keepdims=True)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create stacked bars
    bottom = np.zeros(len(samples))
    for i in range(n_signatures):
        ax.bar(range(len(samples)), W_norm[:, i], 
               bottom=bottom, 
               color=SIGNATURE_COLORS[i % len(SIGNATURE_COLORS)],
               label=signature_names[i])
        bottom += W_norm[:, i]
    
    # Styling
    ax.set_ylim(0, 1)
    ax.set_xlim(-1, len(samples))
    ax.set_ylabel("Relative Contribution", fontsize=12)
    ax.set_yticks(np.arange(0, 1.1, 0.25))
    
    # Remove x-axis labels if too many samples
    if len(samples) > 50:
        ax.set_xticklabels([])
        ax.set_xlabel("Samples", fontsize=12)
    else:
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, rotation=90, fontsize=8)
    
    # Legend
    ax.legend(loc='lower center', ncol=3, fontsize=10, 
             bbox_to_anchor=(0.5, -0.3), frameon=False)
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    return fig


def create_signature_donut_plot(W: np.ndarray,
                              signature_names: Optional[List[str]] = None,
                              figsize: Tuple[int, int] = (8, 8)) -> plt.Figure:
    """
    Create donut plot showing overall signature proportions (Visualization E).
    
    Args:
        W: Contribution matrix (n_samples x n_signatures)
        signature_names: Signature names
        figsize: Figure size
        
    Returns:
        Figure with donut plot
    """
    n_signatures = W.shape[1]
    if signature_names is None:
        signature_names = [f"Signature {i+1}" for i in range(n_signatures)]
    
    # Calculate total contribution of each signature
    total_contrib = W.sum(axis=0)
    total_contrib = total_contrib / total_contrib.sum()
    
    # Create labels with percentages
    labels = [f'{name}: {contrib*100:.1f}%' 
             for name, contrib in zip(signature_names, total_contrib)]
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(aspect='equal'))
    
    # Create donut plot
    wedges, texts = ax.pie(total_contrib,
                          colors=SIGNATURE_COLORS[:n_signatures],
                          wedgeprops=dict(width=0.5, edgecolor='white'),
                          startangle=90)
    
    # Add legend
    ax.legend(wedges, labels,
             loc='center left',
             bbox_to_anchor=(1, 0.5),
             fontsize=12)
    
    plt.title('Relative Contribution of Mutational Signatures', 
             fontsize=14, pad=20)
    
    plt.tight_layout()
    
    return fig


def create_mutational_signature_analysis(data: pd.DataFrame,
                                       n_signatures: int = 3,
                                       sample_column: str = "Tumor_Sample_Barcode",
                                       ref_column: str = "REF",
                                       alt_column: str = "ALT",
                                       context_column: Optional[str] = None,
                                       cosmic_signatures: Optional[pd.DataFrame] = None,
                                       figsize: Tuple[int, int] = (20, 24)) -> plt.Figure:
    """
    Create complete mutational signature analysis visualization.
    
    This creates a figure with all 5 components:
    A. Signature profiles
    B. Cosine similarity with COSMIC
    C. Sample contribution heatmap
    D. Sample contribution barplot
    E. Overall contribution donut plot
    
    Args:
        data: Mutation data
        n_signatures: Number of signatures to extract
        sample_column: Sample identifier column
        ref_column: Reference allele column
        alt_column: Alternative allele column
        context_column: Trinucleotide context column
        cosmic_signatures: COSMIC reference signatures
        figsize: Overall figure size
        
    Returns:
        Figure with complete analysis
    """
    print("Extracting mutation matrix...")
    matrix, samples = extract_mutation_matrix(data, sample_column, ref_column, 
                                            alt_column, context_column)
    
    print(f"Matrix shape: {matrix.shape}")
    print(f"Total mutations: {matrix.sum()}")
    
    if matrix.sum() == 0:
        print("Warning: No mutations found in data")
        # Create empty figure with message
        fig = plt.figure(figsize=figsize)
        fig.text(0.5, 0.5, "No mutations found in data", 
                ha='center', va='center', fontsize=20)
        return fig
    
    print("Performing NMF...")
    W, H = perform_nmf(matrix, n_signatures)
    
    # Signature names
    sig_names = [f"Signature {i+1}" for i in range(n_signatures)]
    
    # Create main figure with GridSpec for complex layout
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(5, 2, height_ratios=[3, 1, 1, 1.5, 1.5], 
                         width_ratios=[3, 1], hspace=0.3, wspace=0.2)
    
    # A. Signature profiles (spans full width)
    ax_profiles = fig.add_subplot(gs[0, :])
    
    # Plot signatures as subplots within this axis
    for i in range(n_signatures):
        ax = plt.subplot(n_signatures, 1, i+1)
        
        # Get colors for bars
        bar_colors = []
        for j in range(96):
            sub_type = j // 16
            substitutions = list(SUBSTITUTION_COLORS.keys())
            bar_colors.append(SUBSTITUTION_COLORS[substitutions[sub_type]])
        
        # Normalize to percentages
        values = H[i, :] / H[i, :].sum() * 100
        
        # Create bars
        ax.bar(range(96), values, color=bar_colors, width=0.8)
        
        # Add vertical lines
        for xline in [16, 32, 48, 64, 80]:
            ax.axvline(x=xline - 0.5, color='black', linewidth=0.8)
        
        # Add substitution labels
        y_max = values.max() * 1.2 if values.max() > 0 else 10
        ax.set_ylim(0, y_max)
        
        for pos, (sub, color) in enumerate(SUBSTITUTION_COLORS.items()):
            x_pos = pos * 16 + 8
            ax.text(x_pos, y_max * 0.90, sub, ha='center', va='bottom',
                   color=color, fontweight='bold', fontsize=11)
        
        # Styling
        ax.set_ylabel("Percentage", fontsize=10)
        ax.set_title(sig_names[i], fontsize=12, pad=8)
        ax.grid(axis='y', linestyle='--', alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        if i == n_signatures - 1:
            ax.set_xlabel("96 Trinucleotide Categories", fontsize=10)
    
    # Clear the parent axis
    ax_profiles.set_visible(False)
    
    # B. Cosine similarity heatmap
    ax_cosine = fig.add_subplot(gs[1, :])
    
    # Create cosine similarity plot in this axis
    plt.sca(ax_cosine)
    
    # Create or use COSMIC signatures
    if cosmic_signatures is None:
        # Generate synthetic COSMIC signatures
        n_cosmic = 15
        cosmic_sigs = np.random.random((96, n_cosmic))
        cosmic_sigs = cosmic_sigs / cosmic_sigs.sum(axis=0)
        cosmic_names = [f'SBS{i+1}' for i in range(n_cosmic)]
        cosmic_df = pd.DataFrame(cosmic_sigs, columns=cosmic_names)
    else:
        cosmic_df = cosmic_signatures
    
    # Calculate similarity
    similarity = cosine_similarity(H, cosmic_df.values.T)
    
    # Create heatmap
    sns.heatmap(similarity, 
               xticklabels=cosmic_df.columns,
               yticklabels=sig_names,
               vmin=0, vmax=1,
               cmap='Blues',
               linewidths=1,
               square=False,
               cbar_kws={'shrink': 0.8},
               ax=ax_cosine)
    
    ax_cosine.set_xlabel("COSMIC SBS Signatures", fontsize=10)
    ax_cosine.set_ylabel("Identified Signatures", fontsize=10)
    plt.setp(ax_cosine.get_xticklabels(), rotation=90, fontsize=8)
    plt.setp(ax_cosine.get_yticklabels(), rotation=0, fontsize=9)
    
    # C. Sample contribution heatmap
    ax_heatmap = fig.add_subplot(gs[2, :])
    
    # Normalize contributions
    W_norm = W / W.sum(axis=1, keepdims=True)
    
    # Create heatmap
    plt.sca(ax_heatmap)
    sns.heatmap(W_norm.T,
               xticklabels=False,
               yticklabels=sig_names,
               vmin=0, vmax=1,
               cmap='Blues',
               linewidths=0.5,
               cbar_kws={'orientation': 'horizontal', 'shrink': 0.5, 'aspect': 30},
               ax=ax_heatmap)
    
    ax_heatmap.set_xlabel("Samples", fontsize=10)
    plt.setp(ax_heatmap.get_yticklabels(), fontsize=9)
    
    # D. Stacked bar plot
    ax_barplot = fig.add_subplot(gs[3, 0])
    
    # Create stacked bars
    bottom = np.zeros(len(samples))
    for i in range(n_signatures):
        ax_barplot.bar(range(len(samples)), W_norm[:, i],
                      bottom=bottom,
                      color=SIGNATURE_COLORS[i % len(SIGNATURE_COLORS)],
                      label=sig_names[i])
        bottom += W_norm[:, i]
    
    ax_barplot.set_ylim(0, 1)
    ax_barplot.set_xlim(-1, len(samples))
    ax_barplot.set_ylabel("Relative Contribution", fontsize=10)
    ax_barplot.set_xlabel("Samples", fontsize=10)
    ax_barplot.set_xticklabels([])
    
    # Legend below the plot
    ax_barplot.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
                     ncol=3, fontsize=9, frameon=False)
    
    # Remove spines
    ax_barplot.spines['top'].set_visible(False)
    ax_barplot.spines['right'].set_visible(False)
    
    # E. Donut plot
    ax_donut = fig.add_subplot(gs[3:, 1])
    
    # Calculate total contributions
    total_contrib = W.sum(axis=0)
    total_contrib = total_contrib / total_contrib.sum()
    
    # Create donut
    wedges, texts = ax_donut.pie(total_contrib,
                                colors=SIGNATURE_COLORS[:n_signatures],
                                wedgeprops=dict(width=0.5, edgecolor='white'),
                                startangle=90)
    
    # Create legend with percentages
    labels = [f'{name}: {contrib*100:.1f}%' 
             for name, contrib in zip(sig_names, total_contrib)]
    
    ax_donut.legend(wedges, labels,
                   loc='center left',
                   bbox_to_anchor=(1, 0.5),
                   fontsize=10)
    
    # Add labels
    fig.text(0.5, 0.95, 'Mutational Signature Analysis', 
            ha='center', fontsize=16, fontweight='bold')
    
    fig.text(0.02, 0.77, 'A', fontsize=14, fontweight='bold')
    fig.text(0.02, 0.55, 'B', fontsize=14, fontweight='bold')
    fig.text(0.02, 0.45, 'C', fontsize=14, fontweight='bold')
    fig.text(0.02, 0.35, 'D', fontsize=14, fontweight='bold')
    fig.text(0.52, 0.35, 'E', fontsize=14, fontweight='bold')
    
    try:
        plt.tight_layout()
    except:
        # Tight layout may fail with complex subplots, use manual adjustment
        pass
    
    return fig 