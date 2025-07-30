import gzip
import logging
from pathlib import Path
from typing import Optional, Dict

import duckdb
import pandas as pd

from .format import format_chr

logger = logging.getLogger(__name__)


def _parse_vep_extra_column(extra_str: str) -> Dict[str, str]:
    """
    Parse the VEP Extra column which contains semicolon-separated key-value pairs.

    Parameters
    ----------
    extra_str : str
        The Extra column content from VEP output

    Returns
    -------
    Dict[str, str]
        Dictionary with parsed key-value pairs
    """
    if pd.isna(extra_str) or extra_str == "-" or extra_str == "":
        return {}

    result = {}
    pairs = extra_str.split(';')

    for pair in pairs:
        if '=' in pair:
            key, value = pair.split('=', 1)
            result[key] = value

    return result


def _create_region_key_from_maf(row: pd.Series) -> str:
    """
    Create a region key from MAF row that matches the VEP Uploaded_variation format.
    
    Generates the same key that VEP uses by:
    - Using format_chr for proper chromosome formatting (23→X, 24→Y)
    - Using the actual End_Position for indels
    - Using proper reference/alternative allele logic

    Parameters
    ----------
    row : pd.Series
        MAF row with Chromosome, Start_Position, End_position, Reference_Allele, 
        Tumor_Seq_Allele1, and Tumor_Seq_Allele2

    Returns
    -------
    str
        Region key in format chr:start-end:len(ref)/alt
    """
    chrom = format_chr(str(row['Chromosome']))  # "chr1…chrX, chrY"
    start = int(row['Start_Position'])
    end = int(row['End_position'])  # Use the real range
    ref = str(row['Reference_Allele'])
    alt = str(row['Tumor_Seq_Allele2'] or row['Tumor_Seq_Allele1'] or "-")

    return f"{chrom}:{start}-{end}:{len(ref)}/{alt}"


def merge_maf_with_vep_annotations(
        maf_file: str | Path,
        vep_file: str | Path,
        output_file: Optional[str | Path] = None,
        compress: bool = False
) -> tuple[pd.DataFrame, Path]:
    """
    Merge MAF file with VEP annotations using pandas and DuckDB for optimization.

    Parameters
    ----------
    maf_file : str | Path
        Path to the original MAF file (.maf or .maf.gz)
    vep_file : str | Path
        Path to the VEP annotation file (.txt)
    output_file : str | Path, optional
        Output file path. If None, creates filename with "_annotated" suffix
    compress : bool, optional
        Whether to compress the output file with gzip (default: False)

    Returns
    -------
    tuple[pd.DataFrame, Path]
        A tuple containing:
        - Merged DataFrame with MAF data and VEP annotations
        - Path to the output file that was created
    """
    maf_file = Path(maf_file)
    vep_file = Path(vep_file)

    if output_file is None:
        # Create output filename with _annotated suffix
        if maf_file.suffix == '.gz':
            stem = maf_file.stem.replace('.maf', '')
            base_name = f"{stem}_VEP_annotated.maf"
        else:
            stem = maf_file.stem
            base_name = f"{stem}_VEP_annotated{maf_file.suffix}"

        # Add .gz extension if compression is requested
        if compress:
            output_file = maf_file.parent / f"{base_name}.gz"
        else:
            output_file = maf_file.parent / base_name
    else:
        output_file = Path(output_file)
        # If compression is requested but output file doesn't end with .gz, add it
        if compress and not str(output_file).endswith('.gz'):
            output_file = output_file.with_suffix(output_file.suffix + '.gz')

    logger.info(f"Reading MAF file: {maf_file}")

    if maf_file.suffix == '.gz':
        with gzip.open(maf_file, 'rt') as f:
            maf_df = pd.read_csv(f, sep='\t', comment='#', low_memory=False)
    else:
        maf_df = pd.read_csv(maf_file, sep='\t', comment='#', low_memory=False)

    logger.info(f"MAF file loaded: {maf_df.shape[0]} rows, {maf_df.shape[1]} columns")

    logger.info(f"Reading VEP file: {vep_file}")

    # Find the header line (starts with #Uploaded_variation)
    header_line = None
    with open(vep_file, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('#Uploaded_variation'):
                header_line = i
                break

    if header_line is None:
        raise ValueError("Could not find VEP header line starting with #Uploaded_variation")

    # Read the file starting from the header line
    vep_df = pd.read_csv(vep_file, sep='\t', skiprows=header_line, low_memory=False)
    if vep_df.columns[0].startswith('#'):
        vep_df.columns = [vep_df.columns[0][1:]] + list(vep_df.columns[1:])

    logger.info(f"VEP file loaded: {vep_df.shape[0]} rows, {vep_df.shape[1]} columns")

    logger.info("Creating region keys for MAF data...")
    maf_df['region_key'] = maf_df.apply(_create_region_key_from_maf, axis=1)

    # Use VEP Uploaded_variation as the key (after removing # prefix)
    vep_df['region_key'] = vep_df['Uploaded_variation']

    logger.info("Parsing VEP Extra column...")
    vep_extra_parsed = vep_df['Extra'].apply(_parse_vep_extra_column)
    extra_df = pd.json_normalize(vep_extra_parsed)

    # Combine VEP data with parsed extra columns
    vep_with_extra = pd.concat([vep_df, extra_df], axis=1)

    # Remove rows without meaningful annotations (only IMPACT=MODIFIER)
    meaningful_annotations = vep_with_extra[
        ~((vep_with_extra['Extra'] == 'IMPACT=MODIFIER') |
          (vep_with_extra['Extra'].isna()) |
          (vep_with_extra['Gene'] == '-'))
    ].copy()

    logger.info(f"Filtered to {meaningful_annotations.shape[0]} meaningful annotations")

    logger.info("Removing VEP duplicates...")
    original_vep_count = len(meaningful_annotations)
    meaningful_annotations = meaningful_annotations.drop_duplicates("region_key", keep="first")
    logger.info(f"Removed {original_vep_count - len(meaningful_annotations)} duplicate VEP entries")

    logger.info("Performing optimized merge with DuckDB...")

    conn = duckdb.connect()
    # Register DataFrames with DuckDB
    conn.register('maf_data', maf_df)
    conn.register('vep_data', meaningful_annotations)

    # Dynamic SQL query
    vep_columns_mapping = {
        'Gene': 'VEP_Gene',
        'Feature': 'VEP_Feature',
        'Feature_type': 'VEP_Feature_type',
        'Consequence': 'VEP_Consequence',
        'cDNA_position': 'VEP_cDNA_position',
        'CDS_position': 'VEP_CDS_position',
        'Protein_position': 'VEP_Protein_position',
        'Amino_acids': 'VEP_Amino_acids',
        'Codons': 'VEP_Codons',
        'Existing_variation': 'VEP_Existing_variation',
        'SYMBOL': 'VEP_SYMBOL',
        'SYMBOL_SOURCE': 'VEP_SYMBOL_SOURCE',
        'HGNC_ID': 'VEP_HGNC_ID',
        'ENSP': 'VEP_ENSP',
        'SWISSPROT': 'VEP_SWISSPROT',
        'TREMBL': 'VEP_TREMBL',
        'UNIPARC': 'VEP_UNIPARC',
        'UNIPROT_ISOFORM': 'VEP_UNIPROT_ISOFORM',
        'DOMAINS': 'VEP_DOMAINS',
        'IMPACT': 'VEP_IMPACT',
        'STRAND': 'VEP_STRAND',
        'DISTANCE': 'VEP_DISTANCE'
    }

    # Select columns that exist in the VEP data
    available_vep_columns = meaningful_annotations.columns.tolist()
    vep_select_clauses = []

    for vep_col, alias in vep_columns_mapping.items():
        if vep_col in available_vep_columns:
            vep_select_clauses.append(f"v.{vep_col} as {alias}")

    vep_select_str = ",\n        ".join(vep_select_clauses)

    # Perform the merge using SQL
    merge_query = f"""
    SELECT 
        m.*,
        {vep_select_str}
    FROM maf_data m
    LEFT JOIN vep_data v ON m.region_key = v.region_key
    """

    result_df = conn.execute(merge_query).df()

    # Remove the temporary region_key column
    result_df = result_df.drop('region_key', axis=1)

    # Remove duplicate information - if VEP_SYMBOL is the same as Hugo_Symbol, don't duplicate
    if 'Hugo_Symbol' in result_df.columns and 'VEP_SYMBOL' in result_df.columns:
        mask = result_df['Hugo_Symbol'] == result_df['VEP_SYMBOL']
        result_df.loc[mask, 'VEP_SYMBOL'] = None

    # Replace missing VEP annotations with empty strings instead of "-" or NaN
    vep_columns = [col for col in result_df.columns if col.startswith('VEP_')]
    for col in vep_columns:
        result_df[col] = result_df[col].fillna("")
        result_df[col] = result_df[col].replace("-", "")

    logger.info(f"Merge completed: {result_df.shape[0]} rows, {result_df.shape[1]} columns")

    logger.info(f"Saving annotated file to: {output_file}")

    # Save file with or without compression
    if compress or str(output_file).endswith('.gz'):
        with gzip.open(output_file, 'wt') as f:
            result_df.to_csv(f, sep='\t', index=False)
    else:
        result_df.to_csv(output_file, sep='\t', index=False)

    conn.close()

    return result_df, output_file
