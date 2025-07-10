import pandas as pd
import duckdb
from pathlib import Path
from typing import Optional, Dict, Any
import gzip
import logging

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

    Parameters
    ----------
    row : pd.Series
        MAF row with Chromosome, Start_Position, Reference_Allele, and Tumor_Seq_Allele2

    Returns
    -------
    str
        Region key in format chr:start-end:ref/alt
    """
    chrom = str(row['Chromosome']).replace('chr', '')
    if not chrom.startswith('chr'):
        chrom = f'chr{chrom}'

    start = int(row['Start_Position'])
    ref = str(row['Reference_Allele'])
    alt = str(row['Tumor_Seq_Allele2']) if pd.notna(row['Tumor_Seq_Allele2']) else str(row['Tumor_Seq_Allele1'])

    # For single nucleotide variants, start and end are the same
    return f"{chrom}:{start}-{start}:1/{alt}"

def merge_maf_with_vep_annotations(
    maf_file: str | Path,
    vep_file: str | Path,
    output_file: Optional[str | Path] = None
) -> pd.DataFrame:
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

    Returns
    -------
    pd.DataFrame
        Merged DataFrame with MAF data and VEP annotations
    """
    maf_file = Path(maf_file)
    vep_file = Path(vep_file)

    if output_file is None:
        # Create output filename with _annotated suffix
        if maf_file.suffix == '.gz':
            stem = maf_file.stem.replace('.maf', '')
            output_file = maf_file.parent / f"{stem}_annotated.maf"
        else:
            stem = maf_file.stem
            output_file = maf_file.parent / f"{stem}_annotated{maf_file.suffix}"
    else:
        output_file = Path(output_file)

    logger.info(f"Reading MAF file: {maf_file}")

    # Read MAF file
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

    # Create region keys for MAF data
    logger.info("Creating region keys for MAF data...")
    maf_df['region_key'] = maf_df.apply(_create_region_key_from_maf, axis=1)

    # Use VEP Uploaded_variation as the key (after removing # prefix)
    vep_df['region_key'] = vep_df['Uploaded_variation']

    # Parse VEP Extra column to extract annotations
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

    logger.info("Performing optimized merge with DuckDB...")

    conn = duckdb.connect()
    # Register DataFrames with DuckDB
    conn.register('maf_data', maf_df)
    conn.register('vep_data', meaningful_annotations)

    # Perform the merge using SQL
    merge_query = """
    SELECT 
        m.*,
        v.Gene as VEP_Gene,
        v.Feature as VEP_Feature,
        v.Feature_type as VEP_Feature_type,
        v.Consequence as VEP_Consequence,
        v.cDNA_position as VEP_cDNA_position,
        v.CDS_position as VEP_CDS_position,
        v.Protein_position as VEP_Protein_position,
        v.Amino_acids as VEP_Amino_acids,
        v.Codons as VEP_Codons,
        v.Existing_variation as VEP_Existing_variation,
        v.SYMBOL as VEP_SYMBOL,
        v.SYMBOL_SOURCE as VEP_SYMBOL_SOURCE,
        v.HGNC_ID as VEP_HGNC_ID,
        v.ENSP as VEP_ENSP,
        v.SWISSPROT as VEP_SWISSPROT,
        v.TREMBL as VEP_TREMBL,
        v.UNIPARC as VEP_UNIPARC,
        v.UNIPROT_ISOFORM as VEP_UNIPROT_ISOFORM,
        v.DOMAINS as VEP_DOMAINS,
        v.IMPACT as VEP_IMPACT,
        v.STRAND as VEP_STRAND,
        v.DISTANCE as VEP_DISTANCE
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

    logger.info(f"Merge completed: {result_df.shape[0]} rows, {result_df.shape[1]} columns")

    logger.info(f"Saving annotated file to: {output_file}")
    result_df.to_csv(output_file, sep='\t', index=False)

    conn.close()

    return result_df