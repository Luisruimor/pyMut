import logging
import re
from typing import Optional, Dict, Tuple

import duckdb
import pandas as pd

from ..utils.database import (
    PfamAnnotationError,
    connect_db
)
from ..utils.fields import col, find_alias

# Configure logger
logger = logging.getLogger(__name__)


# VCF ANNOTATION LOGIC
# =============================================================================
# TODO: VCF annotation logic has been removed
# All VCF processing and VariantAnnotator class functionality is now disabled.
# Only MAF processing with existing VEP annotations is supported.

def annotate_vcf_variants(vcf_df: pd.DataFrame) -> pd.DataFrame:
    """VCF annotation functionality removed."""
    pass
    return vcf_df


class VariantAnnotator:
    """
    VariantAnnotator class removed

    This class handled VCF annotation using VEP CLI.
    All VCF annotation logic has been disabled.
    """

    def __init__(self, *args, **kwargs):
        """VCF annotation disabled."""
        pass

    def annotate(self, df: pd.DataFrame) -> pd.DataFrame:
        """VCF annotation disabled - returns empty DataFrame."""
        # TODO: VCF annotation logic removed
        pass
        return pd.DataFrame()


# UNIPROT IDENTIFIER RESOLUTION
def resolve_uniprot_identifiers(df: pd.DataFrame, uniprot_column: str, db_conn: duckdb.DuckDBPyConnection) -> Tuple[
    pd.DataFrame, Dict[str, int]]:
    """
    Resolve UniProt identifiers to canonical accessions.
    
    Handles three types of values:
    1. Accession (format "P31946") - use as is
    2. Short name (format "1433B_HUMAN") - resolve via short_name column
    3. External identifiers (NP_*, ENSP_*, etc.) - resolve via prot_id column
    
    Args:
        df: DataFrame with UniProt identifiers
        uniprot_column: Name of column containing UniProt identifiers
        db_conn: DuckDB connection
        
    Returns:
        Tuple of (DataFrame with uniprot_resolved column, resolution statistics)
    """
    logger.debug("Resolving UniProt identifiers to canonical accessions...")

    # Initialize statistics
    stats = {
        'total': 0,
        'direct_accession': 0,
        'via_short_name': 0,
        'via_external_id': 0,
        'unresolved': 0
    }

    df_work = df.copy()
    df_work['uniprot_resolved'] = None
    df_work['resolution_method'] = None

    # Process each unique identifier
    unique_ids = df_work[uniprot_column].dropna().unique()
    resolution_cache = {}

    for uniprot_id in unique_ids:
        if pd.isna(uniprot_id) or str(uniprot_id).strip() == '':
            continue

        uniprot_str = str(uniprot_id).strip()
        stats['total'] += 1

        # Check if already a valid accession (starts with letter, contains digits)
        if re.match(r'^[A-Z][0-9A-Z]{5}$', uniprot_str) or re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$', uniprot_str):
            resolution_cache[uniprot_id] = (uniprot_str, 'direct_accession')
            stats['direct_accession'] += 1
            continue

        # Try to resolve via short_name column
        try:
            result = db_conn.execute(
                "SELECT uniprot FROM xref WHERE short_name = ? LIMIT 1",
                [uniprot_str]
            ).fetchone()

            if result:
                resolved_accession = result[0]
                resolution_cache[uniprot_id] = (resolved_accession, 'via_short_name')
                stats['via_short_name'] += 1
                continue
        except Exception as e:
            logger.warning(f"Error querying short_name for {uniprot_str}: {e}")

        # Try to resolve via prot_id column (external identifiers)
        try:
            result = db_conn.execute(
                "SELECT uniprot FROM xref WHERE prot_id = ? LIMIT 1",
                [uniprot_str]
            ).fetchone()

            if result:
                resolved_accession = result[0]
                resolution_cache[uniprot_id] = (resolved_accession, 'via_external_id')
                stats['via_external_id'] += 1
                continue
        except Exception as e:
            logger.warning(f"Error querying prot_id for {uniprot_str}: {e}")

        # Mark as unresolved
        resolution_cache[uniprot_id] = (None, 'unresolved')
        stats['unresolved'] += 1

    # Apply resolutions to DataFrame
    for idx, row in df_work.iterrows():
        uniprot_id = row[uniprot_column]
        if uniprot_id in resolution_cache:
            resolved_accession, method = resolution_cache[uniprot_id]
            df_work.loc[idx, 'uniprot_resolved'] = resolved_accession
            df_work.loc[idx, 'resolution_method'] = method

    # Log resolution summary
    logger.info("UniProt resolution summary:")
    logger.info(f"   Total identifiers processed: {stats['total']:,}")
    logger.info(f"   Direct accessions: {stats['direct_accession']:,}")
    logger.info(f"   Resolved via short_name: {stats['via_short_name']:,}")
    logger.info(f"   Resolved via external ID: {stats['via_external_id']:,}")
    logger.info(f"   Unresolved: {stats['unresolved']:,}")

    return df_work, stats


# MAF PFAM ANNOTATION LOGIC
def _annotate_pfam_sql(df: pd.DataFrame, db_conn: duckdb.DuckDBPyConnection, aa_column: str,
                       uniprot_alias: str) -> pd.DataFrame:
    """Annotate PFAM domains using SQL for larger datasets."""
    logger.debug("Using SQL for PFAM annotation...")

    db_conn.register('variants_temp', df)

    # SQL query for range join - includes seq_start and seq_end coordinates
    query = f"""
    SELECT v.*, p.pfam_id, p.pfam_name, p.seq_start, p.seq_end
    FROM variants_temp v
    LEFT JOIN pfam p ON v.{uniprot_alias} = p.uniprot 
                    AND v.{aa_column} BETWEEN p.seq_start AND p.seq_end
    """

    result_df = db_conn.execute(query).df()
    db_conn.unregister('variants_temp')

    return result_df


def _extract_uniprot_id(self, row):
    """Extract UniProt ID from VEP columns"""
    uniprot_series = col(pd.DataFrame([row]), 'UNIPROT')
    if uniprot_series is not None and pd.notna(uniprot_series.iloc[0]):
        return str(uniprot_series.iloc[0]).split('.')[0]  # Remove version
    return None


def _extract_aa_position(self, row):
    """Extract amino acid position from Protein_Change"""
    protein_change_series = col(pd.DataFrame([row]), 'Protein_Change')
    if protein_change_series is not None and pd.notna(protein_change_series.iloc[0]):
        match = re.search(r'p\.[A-Za-z]*?(\d+)', str(protein_change_series.iloc[0]))
        if match:
            return int(match.group(1))
    return None


def _annotate_with_database(self, df, db_conn, aa_column, uniprot_alias):
    """Use database for precise PFAM annotation with enhanced UniProt resolution"""
    # Create a working copy with extracted columns
    df_work = df.copy()

    # Extract uniprot and aa_pos if they don't exist
    if 'uniprot' not in df_work.columns:
        df_work['uniprot'] = df_work[uniprot_alias].apply(
            lambda x: str(x).split('.')[0] if pd.notna(x) and x != '' else None
        )

    # Resolve UniProt identifiers to canonical accessions
    df_work, resolution_stats = resolve_uniprot_identifiers(df_work, 'uniprot', db_conn)

    # Filter for valid data (must have resolved accession and amino acid position)
    df_valid = df_work.dropna(subset=['uniprot_resolved', aa_column])
    df_valid = df_valid[df_valid['uniprot_resolved'] != '']

    if len(df_valid) == 0:
        logger.warning("No variants with resolved UniProt accessions and amino acid positions")
        # Add empty PFAM columns
        df_work['pfam_id'] = None
        df_work['pfam_name'] = None
        df_work['seq_start'] = None
        df_work['seq_end'] = None
        return df_work

    logger.debug(f"Annotating {len(df_valid)} variants with PFAM domains...")

    # Get all variants with PFAM annotations from SQL query using resolved accessions
    result_df = _annotate_pfam_sql(df_valid, db_conn, aa_column, 'uniprot_resolved')

    # Count successful annotations
    pfam_annotated_count = result_df['pfam_id'].notna().sum()
    logger.info(f"Variantes anotadas con PFAM: {pfam_annotated_count}/{len(result_df)}")

    # Add PFAM columns to the working dataframe
    df_work['pfam_id'] = None
    df_work['pfam_name'] = None
    df_work['seq_start'] = None
    df_work['seq_end'] = None

    # Update working dataframe with PFAM annotations
    for idx in result_df.index:
        if pd.notna(result_df.loc[idx, 'pfam_id']):
            df_work.loc[idx, 'pfam_id'] = result_df.loc[idx, 'pfam_id']
            df_work.loc[idx, 'pfam_name'] = result_df.loc[idx, 'pfam_name']
            df_work.loc[idx, 'seq_start'] = result_df.loc[idx, 'seq_start']
            df_work.loc[idx, 'seq_end'] = result_df.loc[idx, 'seq_end']

    # Store resolution statistics for reporting
    df_work.attrs['resolution_stats'] = resolution_stats

    return df_work


def _annotate_with_vep_domains(self, df):
    """Parse PFAM from VEP_DOMAINS as fallback"""
    logger.debug("Extracting PFAM domains from VEP_DOMAINS column...")

    # Check if VEP_DOMAINS column exists
    domains_series = col(df, 'Domains')
    if domains_series is None:
        logger.warning("VEP_DOMAINS column not found")
        result_data = df.copy()
        for pfam_col in ['pfam_id', 'pfam_name', 'seq_start', 'seq_end']:
            result_data[pfam_col] = None
        return result_data

    # Extract UniProt IDs and PFAM domains from VEP columns
    result_rows = []

    for idx, row in df.iterrows():
        new_row = row.to_dict()

        # Extract UniProt ID if missing
        if 'uniprot' not in new_row or pd.isna(new_row.get('uniprot')):
            uniprot_id = None
            uniprot_series = col(pd.DataFrame([row]), 'UNIPROT')
            if uniprot_series is not None and pd.notna(uniprot_series.iloc[0]) and uniprot_series.iloc[0] != '':
                uniprot_id = str(uniprot_series.iloc[0]).split('.')[0]  # Remove version
            new_row['uniprot'] = uniprot_id

        # Extract amino acid position if missing
        if 'aa_pos' not in new_row or pd.isna(new_row.get('aa_pos')):
            aa_pos = None
            protein_change_series = col(pd.DataFrame([row]), 'Protein_Change')
            if protein_change_series is not None and pd.notna(protein_change_series.iloc[0]):
                match = re.search(r'p\.[A-Za-z]*?(\d+)', str(protein_change_series.iloc[0]))
                if match:
                    aa_pos = int(match.group(1))
            new_row['aa_pos'] = aa_pos

        # Extract PFAM domains from VEP_DOMAINS
        pfam_domains = []
        domains_series = col(pd.DataFrame([row]), 'Domains')
        if domains_series is not None and pd.notna(domains_series.iloc[0]) and domains_series.iloc[0] != '':
            domains_str = str(domains_series.iloc[0])
            if 'Pfam:' in domains_str:
                pfam_matches = re.findall(r'Pfam:([^,\s]+)', domains_str)
                pfam_domains.extend(pfam_matches)

        # Set PFAM information
        if pfam_domains:
            new_row['pfam_id'] = pfam_domains[0]  # Take first domain
            new_row['pfam_name'] = pfam_domains[0]  # Use ID as name for now
            new_row['seq_start'] = None
            new_row['seq_end'] = None
        else:
            new_row['pfam_id'] = None
            new_row['pfam_name'] = None
            new_row['seq_start'] = None
            new_row['seq_end'] = None

        result_rows.append(new_row)

    result_df = pd.DataFrame(result_rows)

    # Show summary
    total_variants = len(result_df)
    with_uniprot = result_df['uniprot'].notna().sum()
    with_aa_pos = result_df['aa_pos'].notna().sum()
    with_pfam = result_df['pfam_id'].notna().sum()

    logger.debug("Processing summary:")
    logger.debug(f"   Total variants: {total_variants:,}")
    logger.debug(f"   With UniProt ID: {with_uniprot:,}")
    logger.debug(f"   With amino acid position: {with_aa_pos:,}")
    logger.debug(f"   With PFAM domains: {with_pfam:,}")

    return result_df


def annotate_pfam(self,
                  db_conn: Optional[duckdb.DuckDBPyConnection] = None,
                  *, aa_column: str = 'aa_pos',
                  auto_extract: bool = True,
                  prefer_database: bool = True):
    """
    Annotate PyMutation data with PFAM domains.

    Automatically detects available data and chooses the best annotation strategy:
    1. If uniprot + aa_pos exist → Use database annotation (most precise)
    2. If VEP data available → Extract uniprot + aa_pos, then use database
    3. Fallback → Extract basic PFAM from VEP_DOMAINS

    Args:
        db_conn: DuckDB connection (if None, will create one)
        aa_column: Name of the column containing amino acid positions
        auto_extract: If True, automatically extract uniprot/aa_pos from VEP data
        prefer_database: If True, prefer database annotation over VEP parsing

    Returns:
        PyMutation: New PyMutation object with PFAM domain annotations
    """
    if db_conn is None:
        db_conn = connect_db()
        close_conn = True
    else:
        close_conn = False

    try:
        df = self.data.copy()

        # Check what data we have
        uniprot_alias = find_alias(df.columns, 'UNIPROT')
        has_aa_pos = aa_column in df.columns
        has_vep_domains = col(df, 'Domains') is not None
        has_protein_change = col(df, 'Protein_Change') is not None

        logger.debug("Data availability check:")
        logger.debug(f"   UniProt column: {'Yes' if uniprot_alias else 'No'}")
        logger.debug(f"   AA position column: {'Yes' if has_aa_pos else 'No'}")
        logger.debug(f"   VEP_DOMAINS: {'Yes' if has_vep_domains else 'No'}")
        logger.debug(f"   Protein_Change: {'Yes' if has_protein_change else 'No'}")

        # Extract missing columns if auto_extract=True
        if auto_extract:
            if uniprot_alias is None and has_vep_domains:
                logger.debug("Extracting UniProt IDs from VEP columns...")
                df['uniprot'] = df.apply(lambda row: _extract_uniprot_id(self, row), axis=1)
                uniprot_alias = 'uniprot'

            if not has_aa_pos and has_protein_change:
                logger.debug("Extracting amino acid positions from Protein_Change...")
                df[aa_column] = df.apply(lambda row: _extract_aa_position(self, row), axis=1)
                has_aa_pos = True

        # Choose annotation strategy
        if uniprot_alias and has_aa_pos and prefer_database:
            logger.debug("Using database annotation (most precise)")
            result_df = _annotate_with_database(self, df, db_conn, aa_column, uniprot_alias)

            # Display resolution summary if available
            if hasattr(result_df, 'attrs') and 'resolution_stats' in result_df.attrs:
                stats = result_df.attrs['resolution_stats']
                logger.info("\nFinal annotation summary:")
                logger.info(f"   Total variants processed: {len(df):,}")
                logger.info(
                    f"   UniProt identifiers resolved: {stats['total'] - stats['unresolved']:,}/{stats['total']:,}")
                logger.info(f"   Variants with PFAM annotations: {result_df['pfam_id'].notna().sum():,}")

        elif has_vep_domains:
            logger.debug("Using VEP_DOMAINS parsing (fallback)")
            result_df = _annotate_with_vep_domains(self, df)

        else:
            logger.warning("No suitable data found for PFAM annotation")
            result_df = df.copy()
            for pfam_col in ['pfam_id', 'pfam_name', 'seq_start', 'seq_end']:
                result_df[pfam_col] = None

        # Return new PyMutation object
        from ..core import PyMutation
        new_pymut = PyMutation(result_df, metadata=self.metadata, samples=self.samples)

        # Store resolution statistics in metadata if available
        if hasattr(result_df, 'attrs') and 'resolution_stats' in result_df.attrs:
            if new_pymut.metadata is not None:
                # Add pfam_resolution_stats as an attribute to existing MutationMetadata object
                new_pymut.metadata.pfam_resolution_stats = result_df.attrs['resolution_stats']

        return new_pymut

    finally:
        if close_conn:
            db_conn.close()


# PFAM DOMAIN SUMMARY FUNCTION
def pfam_domains(self, *, aa_column: str = 'aa_pos', summarize_by: str = 'PfamDomain',
                 top_n: int = 10, include_synonymous: bool = False, plot: bool = False) -> pd.DataFrame:
    """
    Summarize PFAM domain annotations similar to maftools pfamDomains function.

    Args:
        aa_column: Column name containing amino acid positions
        summarize_by: 'PfamDomain' or 'AAPos' - how to group results
        top_n: Number of top results to return
        include_synonymous: Whether to include synonymous variants
        plot: Whether to generate plots (placeholder for now)

    Returns:
        DataFrame with summarized PFAM domain information
    """
    logger.debug(f"Summarizing PFAM domains (summarize_by={summarize_by}, top_n={top_n})")

    # Use self.data instead of df parameter
    df = self.data

    # Detectar la columna pfam_id correcta
    pfam_id_col = None
    pfam_name_col = None

    if 'pfam_id' in df.columns:
        pfam_id_col = 'pfam_id'
        pfam_name_col = 'pfam_name'
    else:
        raise PfamAnnotationError("No se encontraron columnas PFAM. Ejecute primero annotate_pfam()")

    logger.debug(f"Using PFAM columns: {pfam_id_col}, {pfam_name_col}")

    # Filter data if needed
    df_work = df.copy()

    if not include_synonymous:
        # Filter out synonymous variants if Variant_Classification column exists
        variant_class_candidates = ['Variant_Classification', 'variant_classification', 'Mutation_Type']
        variant_class_col = None
        for candidate in variant_class_candidates:
            if candidate in df_work.columns:
                variant_class_col = candidate
                break

        if variant_class_col is not None:
            df_work = df_work[df_work[variant_class_col] != 'Silent']

    # Filter for variants with PFAM annotations using detected column
    df_pfam = df_work.dropna(subset=[pfam_id_col])

    if len(df_pfam) == 0:
        logger.warning("No variants with PFAM domain annotations found")
        return pd.DataFrame()

    logger.debug(f"Found {len(df_pfam)} variants with PFAM domain annotations")

    if summarize_by == 'PfamDomain':
        # Group by PFAM domain
        hugo_candidates = ['Hugo_Symbol', 'hugo_symbol', 'Gene_Symbol', 'gene_symbol', 'gene']
        hugo_alias = None
        for candidate in hugo_candidates:
            if candidate in df_pfam.columns:
                hugo_alias = candidate
                break

        if hugo_alias is None:
            raise PfamAnnotationError("Hugo_Symbol column not found (no alias found)")

        summary = df_pfam.groupby([pfam_id_col, pfam_name_col]).agg({
            hugo_alias: 'nunique',  # Number of unique genes
            aa_column: 'count'  # Number of variants
        }).reset_index()

        summary.columns = ['pfam_id', 'pfam_name', 'n_genes', 'n_variants']
        summary = summary.sort_values('n_variants', ascending=False)

    elif summarize_by == 'AAPos':
        # Group by amino acid position
        uniprot_candidates = ['uniprot', 'UniProt', 'UNIPROT', 'uniprot_id']
        uniprot_alias = None
        for candidate in uniprot_candidates:
            if candidate in df_pfam.columns:
                uniprot_alias = candidate
                break

        hugo_candidates = ['Hugo_Symbol', 'hugo_symbol', 'Gene_Symbol', 'gene_symbol', 'gene']
        hugo_alias = None
        for candidate in hugo_candidates:
            if candidate in df_pfam.columns:
                hugo_alias = candidate
                break

        if uniprot_alias is not None and hugo_alias is not None:
            # Fix the aggregation
            summary = df_pfam.groupby([uniprot_alias, aa_column, pfam_id_col, pfam_name_col]).size().reset_index(
                name='n_variants')
            gene_counts = df_pfam.groupby([uniprot_alias, aa_column, pfam_id_col, pfam_name_col])[
                hugo_alias].nunique().reset_index(name='n_genes')
            summary = summary.merge(gene_counts, on=[uniprot_alias, aa_column, pfam_id_col, pfam_name_col])
            summary = summary.sort_values('n_variants', ascending=False)
        else:
            logger.warning("'UNIPROT' or 'Hugo_Symbol' column not found, cannot group by amino acid position")
            return pd.DataFrame()

    else:
        raise ValueError(f"Invalid summarize_by value: {summarize_by}. Must be 'PfamDomain' or 'AAPos'")

    # Return top N results
    result = summary.head(top_n)

    if plot:
        logger.debug("Plot generation requested but not implemented yet")
        # TODO: Implement plotting functionality
        pass

    return result
