import pandas as pd
import duckdb
from typing import Optional
import re

from .utils.database import (
    PfamAnnotationError,
    connect_db
)
from .utils.fields import col, find_alias



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


# MAF PFAM ANNOTATION LOGIC
# =============================================================================
def _annotate_pfam_sql(df: pd.DataFrame, db_conn: duckdb.DuckDBPyConnection, aa_column: str, uniprot_alias: str) -> pd.DataFrame:
    """Annotate PFAM domains using SQL for larger datasets."""
    print("🔍 Using SQL for PFAM annotation...")

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
    """Use database for precise PFAM annotation"""
    # Create a working copy with extracted columns
    df_work = df.copy()

    # Extract uniprot and aa_pos if they don't exist
    if 'uniprot' not in df_work.columns:
        df_work['uniprot'] = df_work[uniprot_alias].apply(
            lambda x: str(x).split('.')[0] if pd.notna(x) and x != '' else None
        )

    # Filter for valid data
    df_valid = df_work.dropna(subset=['uniprot', aa_column])
    df_valid = df_valid[df_valid['uniprot'] != '']

    if len(df_valid) == 0:
        print("⚠️  No variants with valid UniProt IDs and amino acid positions")
        # Add empty PFAM columns
        df_work['pfam_id'] = None
        df_work['pfam_name'] = None
        df_work['seq_start'] = None
        df_work['seq_end'] = None
        return df_work

    print(f"📊 Annotating {len(df_valid)} variants with PFAM domains...")

    # Get all variants with PFAM annotations from SQL query
    result_df = _annotate_pfam_sql(df_valid, db_conn, aa_column, 'uniprot')

    # Count successful annotations
    pfam_annotated_count = result_df['pfam_id'].notna().sum()
    print(f"✓ Variantes anotadas con PFAM: {pfam_annotated_count}/{len(result_df)}")

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

    return df_work

def _annotate_with_vep_domains(self, df):
    """Parse PFAM from VEP_DOMAINS as fallback"""
    print("🔗 Extracting PFAM domains from VEP_DOMAINS column...")

    # Check if VEP_DOMAINS column exists
    domains_series = col(df, 'Domains')
    if domains_series is None:
        print("⚠️  VEP_DOMAINS column not found")
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

    print(f"📊 Processing summary:")
    print(f"   Total variants: {total_variants:,}")
    print(f"   With UniProt ID: {with_uniprot:,}")
    print(f"   With amino acid position: {with_aa_pos:,}")
    print(f"   With PFAM domains: {with_pfam:,}")

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

        print(f"📊 Data availability check:")
        print(f"   UniProt column: {'✓' if uniprot_alias else '✗'}")
        print(f"   AA position column: {'✓' if has_aa_pos else '✗'}")
        print(f"   VEP_DOMAINS: {'✓' if has_vep_domains else '✗'}")
        print(f"   Protein_Change: {'✓' if has_protein_change else '✗'}")

        # Extract missing columns if auto_extract=True
        if auto_extract:
            if uniprot_alias is None and has_vep_domains:
                print("🔗 Extracting UniProt IDs from VEP columns...")
                df['uniprot'] = df.apply(lambda row: _extract_uniprot_id(self, row), axis=1)
                uniprot_alias = 'uniprot'

            if not has_aa_pos and has_protein_change:
                print("🔗 Extracting amino acid positions from Protein_Change...")
                df[aa_column] = df.apply(lambda row: _extract_aa_position(self, row), axis=1)
                has_aa_pos = True

        # Choose annotation strategy
        if uniprot_alias and has_aa_pos and prefer_database:
            print("🎯 Using database annotation (most precise)")
            result_df = _annotate_with_database(self, df, db_conn, aa_column, uniprot_alias)

        elif has_vep_domains:
            print("🔗 Using VEP_DOMAINS parsing (fallback)")
            result_df = _annotate_with_vep_domains(self, df)

        else:
            print("⚠️  No suitable data found for PFAM annotation")
            result_df = df.copy()
            for pfam_col in ['pfam_id', 'pfam_name', 'seq_start', 'seq_end']:
                result_df[pfam_col] = None

        # Return new PyMutation object
        from .core import PyMutation
        return PyMutation(result_df, metadata=self.metadata, samples=self.samples)

    finally:
        if close_conn:
            db_conn.close()


# PFAM DOMAIN SUMMARY FUNCTION
# =============================================================================
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
    print(f"📊 Summarizing PFAM domains (summarize_by={summarize_by}, top_n={top_n})")

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

    print(f"📊 Using PFAM columns: {pfam_id_col}, {pfam_name_col}")

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
        print("⚠️  No variants with PFAM domain annotations found")
        return pd.DataFrame()

    print(f"📊 Found {len(df_pfam)} variants with PFAM domain annotations")

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
            aa_column: 'count'         # Number of variants
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
            summary = df_pfam.groupby([uniprot_alias, aa_column, pfam_id_col, pfam_name_col]).size().reset_index(name='n_variants')
            gene_counts = df_pfam.groupby([uniprot_alias, aa_column, pfam_id_col, pfam_name_col])[hugo_alias].nunique().reset_index(name='n_genes')
            summary = summary.merge(gene_counts, on=[uniprot_alias, aa_column, pfam_id_col, pfam_name_col])
            summary = summary.sort_values('n_variants', ascending=False)
        else:
            print("⚠️  'UNIPROT' or 'Hugo_Symbol' column not found, cannot group by amino acid position")
            return pd.DataFrame()

    else:
        raise ValueError(f"Invalid summarize_by value: {summarize_by}. Must be 'PfamDomain' or 'AAPos'")

    # Return top N results
    result = summary.head(top_n)

    if plot:
        print("📊 Plot generation requested but not implemented yet")
        # TODO: Implement plotting functionality
        pass

    return result