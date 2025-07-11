import pandas as pd
import duckdb
import gzip
import logging
from pathlib import Path
from typing import Optional, Union

from ..utils.fields import find_alias, col

logger = logging.getLogger(__name__)


def maf_COSMIC_annotation(
    maf_file: Union[str, Path],
    annotation_table: Union[str, Path],
    output_path: Optional[Union[str, Path]] = None,
    compress_output: bool = True,
    join_column: str = "Hugo_Symbol",
    synonyms_column: str = "SYNONIMS"
) -> tuple[pd.DataFrame, Path]:
    """
    Annotate MAF file with COSMIC Cancer Gene Census data using optimized approach.

    Uses a hybrid approach: DuckDB for large files (>2GB), pandas+pyarrow for smaller files.

    Parameters
    ----------
    maf_file : str | Path
        Path to the MAF file (.maf or .maf.gz)
    annotation_table : str | Path
        Path to the COSMIC annotation table (.tsv or .tsv.gz)
    output_path : str | Path, optional
        Output file path. If None, creates filename with "_COSMIC_annotated.maf" suffix
    compress_output : bool, default True
        Whether to compress the output file with gzip
    join_column : str, default "Hugo_Symbol"
        Column name to use for joining (canonical name from fields.py)
    synonyms_column : str, default "SYNONIMS"
        Column name containing gene synonyms separated by commas. If this column
        exists in the annotation table, a synonym dictionary will be created to
        improve gene matching.

    Returns
    -------
    tuple[pd.DataFrame, Path]
        A tuple containing:
        - Annotated DataFrame with MAF data and COSMIC annotations
        - Path to the output file that was created

    Raises
    ------
    FileNotFoundError
        If input files don't exist
    ValueError
        If join column is not found in either file
    """

    maf_file = Path(maf_file)
    annotation_table = Path(annotation_table)

    logger.info(f"Starting COSMIC annotation process")
    logger.info(f"MAF file: {maf_file}")
    logger.info(f"Annotation table: {annotation_table}")

    if not maf_file.exists():
        logger.error(f"MAF file not found: {maf_file}")
        raise FileNotFoundError(f"MAF file not found: {maf_file}")
    if not annotation_table.exists():
        logger.error(f"Annotation table not found: {annotation_table}")
        raise FileNotFoundError(f"Annotation table not found: {annotation_table}")
    if output_path is None:
        if maf_file.suffix == '.gz':
            stem = maf_file.stem.replace('.maf', '')
            base_name = f"{stem}_COSMIC_annotated.maf"
        else:
            stem = maf_file.stem
            base_name = f"{stem}_COSMIC_annotated{maf_file.suffix}"

        if compress_output:
            output_path = maf_file.parent / f"{base_name}.gz"
        else:
            output_path = maf_file.parent / base_name
    else:
        output_path = Path(output_path)
        if compress_output and not str(output_path).endswith('.gz'):
            output_path = output_path.with_suffix(output_path.suffix + '.gz')

    maf_size_gb = maf_file.stat().st_size / (1024**3)
    use_duckdb = maf_size_gb > 2.0

    logger.info(f"MAF file size: {maf_size_gb:.2f} GB")
    logger.info(f"Using {'DuckDB' if use_duckdb else 'pandas+pyarrow'} approach")

    if use_duckdb:
        return _annotate_with_duckdb(
            maf_file, annotation_table, output_path, compress_output, join_column, synonyms_column
        )
    else:
        return _annotate_with_pandas(
            maf_file, annotation_table, output_path, compress_output, join_column, synonyms_column
        )


def _read_file_auto(file_path: Path, **kwargs) -> pd.DataFrame:
    """
    Automatically read file detecting .gz compression.

    Parameters
    ----------
    file_path : Path
        Path to the file to read
    **kwargs
        Additional arguments passed to pd.read_csv

    Returns
    -------
    pd.DataFrame
        Loaded DataFrame
    """
    if file_path.suffix == '.gz':
        with gzip.open(file_path, 'rt') as f:
            return pd.read_csv(f, **kwargs)
    else:
        return pd.read_csv(file_path, **kwargs)


def _write_file_auto(df: pd.DataFrame, file_path: Path, compress: bool = False) -> None:
    """
    Automatically write file with optional compression.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to write
    file_path : Path
        Output file path
    compress : bool, default False
        Whether to compress the output
    """
    if compress or str(file_path).endswith('.gz'):
        with gzip.open(file_path, 'wt') as f:
            df.to_csv(f, sep='\t', index=False)
    else:
        df.to_csv(file_path, sep='\t', index=False)


def _create_synonyms_dict(annotation_df: pd.DataFrame, synonyms_column: str) -> dict:
    """
    Create a dictionary mapping gene synonyms to their main gene symbol.

    Parameters
    ----------
    annotation_df : pd.DataFrame
        Annotation DataFrame containing gene symbols and synonyms
    synonyms_column : str
        Name of the column containing synonyms separated by commas

    Returns
    -------
    dict
        Dictionary mapping synonyms to main gene symbols
    """
    synonyms_dict = {}

    if synonyms_column not in annotation_df.columns:
        logger.info(f"Synonyms column '{synonyms_column}' not found in annotation table. Skipping synonym mapping.")
        return synonyms_dict

    logger.info(f"Creating synonyms dictionary from column '{synonyms_column}'...")

    for _, row in annotation_df.iterrows():
        gene_symbol = row['GENE_SYMBOL']
        synonyms_str = row[synonyms_column]

        if pd.notna(synonyms_str) and synonyms_str.strip():
            synonyms = [syn.strip() for syn in str(synonyms_str).split(',') if syn.strip()]
            for synonym in synonyms:
                synonyms_dict[synonym] = gene_symbol

    logger.info(f"Created synonyms dictionary with {len(synonyms_dict)} mappings")
    return synonyms_dict


def _apply_synonyms_mapping(maf_df: pd.DataFrame, maf_join_col: str, synonyms_dict: dict) -> pd.DataFrame:
    """
    Apply synonyms mapping to MAF DataFrame to improve gene matching.

    Parameters
    ----------
    maf_df : pd.DataFrame
        MAF DataFrame
    maf_join_col : str
        Column name used for joining in MAF
    synonyms_dict : dict
        Dictionary mapping synonyms to main gene symbols

    Returns
    -------
    pd.DataFrame
        MAF DataFrame with an additional column for mapped gene symbols
    """
    if not synonyms_dict:
        maf_df = maf_df.copy()
        maf_df['_mapped_gene_symbol'] = maf_df[maf_join_col]
        return maf_df

    logger.info("Applying synonyms mapping to MAF data...")

    maf_df = maf_df.copy()
    maf_df['_mapped_gene_symbol'] = maf_df[maf_join_col].map(
        lambda x: synonyms_dict.get(x, x) if pd.notna(x) else x
    )

    direct_matches = (maf_df[maf_join_col] == maf_df['_mapped_gene_symbol']).sum()
    synonym_matches = len(maf_df) - direct_matches

    logger.info(f"Gene mapping results: {direct_matches} direct matches, {synonym_matches} synonym matches")

    return maf_df


def _annotate_with_pandas(
    maf_file: Path,
    annotation_table: Path,
    output_path: Path,
    compress_output: bool,
    join_column: str,
    synonyms_column: str
) -> tuple[pd.DataFrame, Path]:
    """
    Annotate using pandas with pyarrow optimization for smaller files.
    """
    logger.info(f"Reading MAF file: {maf_file}")

    maf_df = _read_file_auto(maf_file, sep='\t', comment='#', low_memory=False)
    logger.info(f"MAF file loaded: {maf_df.shape[0]} rows, {maf_df.shape[1]} columns")

    maf_join_col = find_alias(maf_df.columns, join_column)
    if maf_join_col is None:
        logger.error(f"Join column '{join_column}' not found in MAF file. Available columns: {list(maf_df.columns)}")
        raise ValueError(f"Join column '{join_column}' not found in MAF file. Available columns: {list(maf_df.columns)}")

    logger.info(f"Using MAF join column: {maf_join_col}")

    logger.info(f"Reading annotation table: {annotation_table}")
    annotation_df = _read_file_auto(annotation_table, sep='\t', low_memory=False)
    logger.info(f"Annotation table loaded: {annotation_df.shape[0]} rows, {annotation_df.shape[1]} columns")

    if 'GENE_SYMBOL' not in annotation_df.columns:
        logger.error(f"'GENE_SYMBOL' column not found in annotation table. Available columns: {list(annotation_df.columns)}")
        raise ValueError(f"'GENE_SYMBOL' column not found in annotation table. Available columns: {list(annotation_df.columns)}")

    synonyms_dict = _create_synonyms_dict(annotation_df, synonyms_column)
    maf_df_mapped = _apply_synonyms_mapping(maf_df, maf_join_col, synonyms_dict)

    logger.info("Performing annotation merge...")

    # Add prefix to annotation columns to avoid conflicts (except GENE_SYMBOL)
    annotation_cols_to_rename = {col: f"COSMIC_{col}" for col in annotation_df.columns if col != 'GENE_SYMBOL'}
    annotation_df_renamed = annotation_df.rename(columns=annotation_cols_to_rename)

    result_df = maf_df_mapped.merge(
        annotation_df_renamed,
        left_on='_mapped_gene_symbol',
        right_on='GENE_SYMBOL',
        how='left'
    )

    columns_to_drop = []
    if 'GENE_SYMBOL' in result_df.columns:
        columns_to_drop.append('GENE_SYMBOL')
    if '_mapped_gene_symbol' in result_df.columns:
        columns_to_drop.append('_mapped_gene_symbol')

    if columns_to_drop:
        result_df = result_df.drop(columns_to_drop, axis=1)

    cosmic_columns = [col for col in result_df.columns if col.startswith('COSMIC_')]
    for col in cosmic_columns:
        result_df[col] = result_df[col].fillna("")

    logger.info(f"Annotation completed: {result_df.shape[0]} rows, {result_df.shape[1]} columns")
    logger.info(f"Added {len(cosmic_columns)} COSMIC annotation columns")

    logger.info(f"Saving annotated file to: {output_path}")
    _write_file_auto(result_df, output_path, compress_output)

    logger.info(f"COSMIC annotation completed successfully")
    logger.info(f"Output file saved: {output_path}")

    return result_df, output_path


def _annotate_with_duckdb(
    maf_file: Path,
    annotation_table: Path,
    output_path: Path,
    compress_output: bool,
    join_column: str,
    synonyms_column: str
) -> tuple[pd.DataFrame, Path]:
    """
    Annotate using DuckDB for large files optimization.
    """
    logger.info(f"Reading MAF file: {maf_file}")

    maf_df = _read_file_auto(maf_file, sep='\t', comment='#', low_memory=False)
    logger.info(f"MAF file loaded: {maf_df.shape[0]} rows, {maf_df.shape[1]} columns")

    maf_join_col = find_alias(maf_df.columns, join_column)
    if maf_join_col is None:
        logger.error(f"Join column '{join_column}' not found in MAF file. Available columns: {list(maf_df.columns)}")
        raise ValueError(f"Join column '{join_column}' not found in MAF file. Available columns: {list(maf_df.columns)}")

    logger.info(f"Using MAF join column: {maf_join_col}")

    logger.info(f"Reading annotation table: {annotation_table}")
    annotation_df = _read_file_auto(annotation_table, sep='\t', low_memory=False)
    logger.info(f"Annotation table loaded: {annotation_df.shape[0]} rows, {annotation_df.shape[1]} columns")

    if 'GENE_SYMBOL' not in annotation_df.columns:
        logger.error(f"'GENE_SYMBOL' column not found in annotation table. Available columns: {list(annotation_df.columns)}")
        raise ValueError(f"'GENE_SYMBOL' column not found in annotation table. Available columns: {list(annotation_df.columns)}")

    synonyms_dict = _create_synonyms_dict(annotation_df, synonyms_column)
    maf_df_mapped = _apply_synonyms_mapping(maf_df, maf_join_col, synonyms_dict)

    logger.info("Performing optimized merge with DuckDB...")

    conn = duckdb.connect()

    try:
        conn.register('maf_data', maf_df_mapped)
        conn.register('annotation_data', annotation_df)

        # Build select clauses for annotation columns (with COSMIC_ prefix)
        annotation_columns = [col for col in annotation_df.columns if col != 'GENE_SYMBOL']
        annotation_select_clauses = [f"a.{col} as COSMIC_{col}" for col in annotation_columns]
        annotation_select_str = ",\n        ".join(annotation_select_clauses)

        merge_query = f"""
        SELECT 
            m.*,
            {annotation_select_str}
        FROM maf_data m
        LEFT JOIN annotation_data a ON m._mapped_gene_symbol = a.GENE_SYMBOL
        """

        result_df = conn.execute(merge_query).df()

        if '_mapped_gene_symbol' in result_df.columns:
            result_df = result_df.drop('_mapped_gene_symbol', axis=1)

        cosmic_columns = [col for col in result_df.columns if col.startswith('COSMIC_')]
        for col in cosmic_columns:
            result_df[col] = result_df[col].fillna("")

        logger.info(f"Annotation completed: {result_df.shape[0]} rows, {result_df.shape[1]} columns")
        logger.info(f"Added {len(cosmic_columns)} COSMIC annotation columns")

        logger.info(f"Saving annotated file to: {output_path}")
        _write_file_auto(result_df, output_path, compress_output)

        logger.info(f"COSMIC annotation completed successfully")
        logger.info(f"Output file saved: {output_path}")

        return result_df, output_path

    finally:
        conn.close()


__all__ = [
    "maf_COSMIC_annotation"
]
