"""
Mutational signature analysis module for PyMutation.

This module provides functionality for analyzing trinucleotide contexts
and generating mutational signature matrices.
"""

import pandas as pd
import numpy as np
from typing import Tuple, Dict, Optional
import logging

# Set up logging
logger = logging.getLogger(__name__)

# Define the 96 trinucleotide contexts in standard order
TRINUCLEOTIDE_CONTEXTS = [
    # C>A mutations
    "A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T",
    "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T",
    "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T",
    "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T",

    # C>G mutations
    "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T",
    "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T",
    "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T",
    "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T",

    # C>T mutations
    "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
    "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
    "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T",
    "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",

    # T>A mutations
    "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T",
    "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T",
    "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T",
    "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T",

    # T>C mutations
    "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T",
    "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T",
    "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T",
    "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T",

    # T>G mutations
    "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
    "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
    "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T",
    "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T"
]

# Create mapping from context to index
CONTEXT_TO_INDEX = {context: idx for idx, context in enumerate(TRINUCLEOTIDE_CONTEXTS)}

# Complement mapping
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def _get_reverse_complement(sequence: str) -> str:
    """Get reverse complement of a DNA sequence."""
    return ''.join(COMPLEMENT[base] for base in reversed(sequence))


def _normalize_to_pyrimidine(ref: str, alt: str, trinuc: str) -> Tuple[str, str, str]:
    """
    Normalize mutation to pyrimidine context (C or T as reference).

    Parameters
    ----------
    ref : str
        Reference allele
    alt : str
        Alternative allele
    trinuc : str
        Trinucleotide context

    Returns
    -------
    Tuple[str, str, str]
        Normalized (ref, alt, trinuc)
    """
    if ref in ['C', 'T']:
        return ref, alt, trinuc
    else:
        # Convert to reverse complement
        new_ref = COMPLEMENT[ref]
        new_alt = COMPLEMENT[alt]
        new_trinuc = _get_reverse_complement(trinuc)
        return new_ref, new_alt, new_trinuc


def _get_trinucleotide_context(fasta, chrom: str, pos: int) -> Optional[str]:
    """
    Extract trinucleotide context from FASTA file.

    Parameters
    ----------
    fasta : pyfaidx.Fasta
        FASTA file object
    chrom : str
        Chromosome name
    pos : int
        Position (1-based)

    Returns
    -------
    Optional[str]
        Trinucleotide context or None if not available
    """
    try:
        # Handle chromosome name variations
        chrom_key = chrom
        if chrom_key not in fasta.keys():
            # Try with 'chr' prefix
            if not chrom_key.startswith('chr'):
                chrom_key = f'chr{chrom}'
            # Try without 'chr' prefix
            elif chrom_key.startswith('chr'):
                chrom_key = chrom_key[3:]

        if chrom_key not in fasta.keys():
            logger.warning(f"Chromosome {chrom} not found in FASTA file")
            return None

        # Extract trinucleotide (1-based coordinates)
        trinuc = fasta[chrom_key][pos-2:pos+1].seq.upper()

        if len(trinuc) != 3 or 'N' in trinuc:
            return None

        return trinuc

    except Exception as e:
        logger.warning(f"Error extracting trinucleotide for {chrom}:{pos}: {e}")
        return None


def _create_context_label(ref: str, alt: str, trinuc: str) -> str:
    """
    Create the 96-context label.

    Parameters
    ----------
    ref : str
        Reference allele
    alt : str
        Alternative allele
    trinuc : str
        Trinucleotide context

    Returns
    -------
    str
        Context label in format "X[REF>ALT]Z"
    """
    return f"{trinuc[0]}[{ref}>{alt}]{trinuc[2]}"


def trinucleotideMatrix(self, fasta_file: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Generate trinucleotide context matrix for mutational signature analysis.

    This method calculates the 96 trinucleotide contexts for all SNVs in the dataset
    and returns both a contexts matrix (96 x samples) and the original data enriched
    with trinucleotide information.

    Parameters
    ----------
    fasta_file : str
        Path to the reference genome FASTA file

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        - contexts_df: 96 x samples matrix with trinucleotide context counts
        - enriched_data: Original data with added trinuc, class96, and idx96 columns

    Raises
    ------
    ImportError
        If pyfaidx is not installed
    ValueError
        If required columns are missing or no valid SNVs found
    """
    try:
        import pyfaidx
    except ImportError:
        raise ImportError("pyfaidx is required for trinucleotide analysis. Install with: pip install pyfaidx")

    # Import field utilities
    from ..utils.fields import col

    # Get required columns using the field mapping system
    chrom_col = col(self.data, "Chromosome", required=True)
    pos_col = col(self.data, "Start_Position", required=True) 
    ref_col = col(self.data, "Reference_Allele", required=True)
    alt_col = col(self.data, "Tumor_Seq_Allele2", required=True)
    sample_col = col(self.data, "Tumor_Sample_Barcode", required=False)

    if any(col is None for col in [chrom_col, pos_col, ref_col, alt_col]):
        raise ValueError("Required columns not found. Need: Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2")

    # Detect data format (wide vs long)
    if sample_col is None:
        # Wide format: each sample is a separate column
        # Identify sample columns (exclude standard VCF/MAF columns)
        standard_cols = {'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                        'Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Start_Position', 
                        'End_position', 'Strand', 'Variant_Classification', 'Variant_Type', 
                        'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 
                        'Tumor_Sample_Barcode', 'Protein_Change', 'i_TumorVAF_WU', 'i_transcript_name'}

        sample_columns = [col for col in self.data.columns if col not in standard_cols]
        logger.info(f"Detected wide format with {len(sample_columns)} sample columns")
        data_format = "wide"
    else:
        # Long format: single sample column
        sample_columns = None
        logger.info("Detected long format with Tumor_Sample_Barcode column")
        data_format = "long"

    # Filter for SNVs only
    snv_mask = (
        (ref_col.isin(['A', 'C', 'G', 'T'])) & 
        (alt_col.isin(['A', 'C', 'G', 'T'])) &
        (ref_col != alt_col)
    )

    if not snv_mask.any():
        raise ValueError("No valid SNVs found in the dataset")

    # Create working dataframe with SNVs only
    snv_data = self.data[snv_mask].copy()
    logger.info(f"Processing {len(snv_data)} SNVs from {len(self.data)} total mutations")

    # Load FASTA file
    try:
        fasta = pyfaidx.Fasta(fasta_file)
        logger.info(f"Loaded FASTA file: {fasta_file}")
    except Exception as e:
        raise ValueError(f"Error loading FASTA file {fasta_file}: {e}")

    # Extract trinucleotide contexts
    trinuc_contexts = []
    class96_labels = []
    idx96_values = []

    for idx, row in snv_data.iterrows():
        chrom = str(row[chrom_col.name])
        pos = int(row[pos_col.name])
        ref = str(row[ref_col.name]).upper()
        alt = str(row[alt_col.name]).upper()

        # Get trinucleotide context
        trinuc = _get_trinucleotide_context(fasta, chrom, pos)

        if trinuc is None:
            trinuc_contexts.append(None)
            class96_labels.append(None)
            idx96_values.append(None)
            continue

        # Normalize to pyrimidine context
        norm_ref, norm_alt, norm_trinuc = _normalize_to_pyrimidine(ref, alt, trinuc)

        # Create context label
        context_label = _create_context_label(norm_ref, norm_alt, norm_trinuc)

        # Get index
        context_idx = CONTEXT_TO_INDEX.get(context_label)

        trinuc_contexts.append(norm_trinuc)
        class96_labels.append(context_label)
        idx96_values.append(context_idx)

    # Add new columns to the SNV data
    snv_data = snv_data.copy()
    snv_data['trinuc'] = trinuc_contexts
    snv_data['class96'] = class96_labels
    snv_data['idx96'] = idx96_values

    # Filter out rows with missing context information
    valid_mask = snv_data['idx96'].notna()
    valid_data = snv_data[valid_mask].copy()

    logger.info(f"Successfully processed {len(valid_data)} SNVs with valid trinucleotide contexts")

    if len(valid_data) == 0:
        raise ValueError("No SNVs with valid trinucleotide contexts found")

    # Create the 96 x samples matrix
    if data_format == "long":
        # Long format: use sample column
        sample_names = valid_data[sample_col.name].unique()
        contexts_matrix = np.zeros((96, len(sample_names)), dtype=int)

        # Fill the matrix
        for sample_idx, sample in enumerate(sample_names):
            sample_data = valid_data[valid_data[sample_col.name] == sample]
            context_counts = sample_data['idx96'].value_counts()

            for context_idx, count in context_counts.items():
                if pd.notna(context_idx):
                    contexts_matrix[int(context_idx), sample_idx] = count

    else:
        # Wide format: each sample is a column with genotype calls
        sample_names = sample_columns
        contexts_matrix = np.zeros((96, len(sample_names)), dtype=int)

        # For each variant (row), check which samples have the mutation
        for idx, row in valid_data.iterrows():
            context_idx = row['idx96']
            if pd.notna(context_idx):
                context_idx = int(context_idx)

                # Check each sample column for this variant
                for sample_idx, sample_name in enumerate(sample_names):
                    genotype = str(row[sample_name])

                    # Count mutations based on genotype
                    # Genotypes like "G|A", "A|G" indicate heterozygous mutation
                    # Genotypes like "A|A" indicate homozygous mutation
                    # Genotypes like "G|G" indicate no mutation
                    if '|' in genotype:
                        alleles = genotype.split('|')
                        ref_allele = str(row[ref_col.name]).upper()
                        alt_allele = str(row[alt_col.name]).upper()

                        # Count how many copies of the alternative allele
                        alt_count = sum(1 for allele in alleles if allele == alt_allele)
                        contexts_matrix[context_idx, sample_idx] += alt_count

    # Create contexts DataFrame
    contexts_df = pd.DataFrame(
        contexts_matrix,
        index=TRINUCLEOTIDE_CONTEXTS,
        columns=sample_names
    )

    logger.info(f"Generated {contexts_df.shape[0]} x {contexts_df.shape[1]} trinucleotide context matrix")

    return contexts_df, valid_data


# Add the method to PyMutation class
def add_trinucleotide_method_to_pymutation():
    """Add trinucleotideMatrix method to PyMutation class."""
    from ..core import PyMutation
    PyMutation.trinucleotideMatrix = trinucleotideMatrix
