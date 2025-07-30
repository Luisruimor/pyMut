import gzip
import logging
from pathlib import Path
from typing import Optional, Union

import duckdb
import pandas as pd

from ..utils.fields import find_alias

logger = logging.getLogger(__name__)


def _maf_COSMIC_OncoKB_annotation_aux(
        maf_file: Union[str, Path],
        annotation_table: Union[str, Path],
        output_path: Optional[Union[str, Path]] = None,
        compress_output: bool = True,
        join_column: str = "Hugo_Symbol",
        oncokb_table: Optional[Union[str, Path]] = None
) -> tuple[pd.DataFrame, Path]:
    """
    Annotate MAF file with COSMIC Cancer Gene Census data and optionally OncoKB data.

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
        Column name to use for joining
    oncokb_table : str | Path, optional
        Path to the OncoKB cancer gene list table (.tsv)

    Returns
    -------
    tuple[pd.DataFrame, Path]
        Annotated DataFrame and output file path
    """

    maf_file = Path(maf_file)
    annotation_table = Path(annotation_table)

    oncokb_table_path = None
    if oncokb_table is not None:
        oncokb_table_path = Path(oncokb_table)

    logger.info("Starting COSMIC annotation process")
    logger.info(f"MAF file: {maf_file}")
    logger.info(f"Annotation table: {annotation_table}")
    if oncokb_table_path:
        logger.info(f"OncoKB table: {oncokb_table_path}")

    if not maf_file.exists():
        logger.error(f"MAF file not found: {maf_file}")
        raise FileNotFoundError(f"MAF file not found: {maf_file}")
    if not annotation_table.exists():
        logger.error(f"Annotation table not found: {annotation_table}")
        raise FileNotFoundError(f"Annotation table not found: {annotation_table}")
    if oncokb_table_path and not oncokb_table_path.exists():
        logger.error(f"OncoKB table not found: {oncokb_table_path}")
        raise FileNotFoundError(f"OncoKB table not found: {oncokb_table_path}")
    if output_path is None:
        if maf_file.suffix == '.gz':
            stem = maf_file.stem.replace('.maf', '')
            if oncokb_table_path:
                base_name = f"{stem}_COSMIC_OncoKB_annotated.maf"
            else:
                base_name = f"{stem}_COSMIC_annotated.maf"
        else:
            stem = maf_file.stem
            if oncokb_table_path:
                base_name = f"{stem}_COSMIC_OncoKB_annotated{maf_file.suffix}"
            else:
                base_name = f"{stem}_COSMIC_annotated{maf_file.suffix}"

        if compress_output:
            output_path = maf_file.parent / f"{base_name}.gz"
        else:
            output_path = maf_file.parent / base_name
    else:
        output_path = Path(output_path)
        if compress_output and not str(output_path).endswith('.gz'):
            output_path = output_path.with_suffix(output_path.suffix + '.gz')

    maf_size_gb = maf_file.stat().st_size / (1024 ** 3)
    use_duckdb = maf_size_gb > 2.0

    logger.info(f"MAF file size: {maf_size_gb:.2f} GB")
    logger.info(f"Using {'DuckDB' if use_duckdb else 'pandas+pyarrow'} approach")

    if use_duckdb:
        return _annotate_with_duckdb(
            maf_file, annotation_table, output_path, compress_output, join_column, "SYNONIMS",
            oncokb_table_path, "Gene Aliases"
        )
    else:
        return _annotate_with_pandas(
            maf_file, annotation_table, output_path, compress_output, join_column, "SYNONIMS",
            oncokb_table_path, "Gene Aliases"
        )


def _read_file_auto(file_path: Path, **kwargs) -> pd.DataFrame:
    """Automatically read file detecting .gz compression."""
    if file_path.suffix == '.gz':
        with gzip.open(file_path, 'rt') as f:
            return pd.read_csv(f, **kwargs)
    else:
        return pd.read_csv(file_path, **kwargs)


def _write_file_auto(df: pd.DataFrame, file_path: Path, compress: bool = False) -> None:
    """Automatically write file with optional compression."""
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


def _create_oncokb_synonyms_dict(oncokb_df: pd.DataFrame, synonyms_column: str) -> dict:
    """
    Create a dictionary mapping OncoKB gene synonyms to their main gene symbol.

    Parameters
    ----------
    oncokb_df : pd.DataFrame
        OncoKB DataFrame containing gene symbols and synonyms
    synonyms_column : str
        Name of the column containing synonyms separated by commas

    Returns
    -------
    dict
        Dictionary mapping synonyms to main gene symbols
    """
    synonyms_dict = {}

    if synonyms_column not in oncokb_df.columns:
        logger.info(f"OncoKB synonyms column '{synonyms_column}' not found. Skipping OncoKB synonym mapping.")
        return synonyms_dict

    logger.info(f"Creating OncoKB synonyms dictionary from column '{synonyms_column}'...")

    for _, row in oncokb_df.iterrows():
        gene_symbol = row['Hugo Symbol']
        synonyms_str = row[synonyms_column]

        if pd.notna(synonyms_str) and synonyms_str.strip():
            synonyms = [syn.strip() for syn in str(synonyms_str).split(',') if syn.strip()]
            for synonym in synonyms:
                synonyms_dict[synonym] = gene_symbol

    logger.info(f"Created OncoKB synonyms dictionary with {len(synonyms_dict)} mappings")
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

    logger.info("Applying synonyms mapping to PyMutation data...")

    maf_df = maf_df.copy()
    maf_df['_mapped_gene_symbol'] = maf_df[maf_join_col].map(
        lambda x: synonyms_dict.get(x, x) if pd.notna(x) else x
    )

    direct_matches = (maf_df[maf_join_col] == maf_df['_mapped_gene_symbol']).sum()
    synonym_matches = len(maf_df) - direct_matches

    logger.info(f"Gene mapping results: {direct_matches} direct matches, {synonym_matches} synonym matches")

    return maf_df


def _annotate_with_pandas(
        data: pd.DataFrame,
        annotation_table: Path,
        join_column: str,
        synonyms_column: str,
        oncokb_table: Optional[Path] = None,
        oncokb_synonyms_column: str = "Gene Aliases"
) -> pd.DataFrame:
    """
    Annotate using pandas with pyarrow optimization for smaller files.
    """
    logger.info(f"Starting pandas annotation for DataFrame: {data.shape[0]} rows, {data.shape[1]} columns")

    maf_df = data.copy()

    maf_join_col = find_alias(maf_df.columns, join_column)
    if maf_join_col is None:
        logger.error(f"Join column '{join_column}' not found in data file. Available columns: {list(maf_df.columns)}")
        raise ValueError(
            f"Join column '{join_column}' not found in data file. Available columns: {list(maf_df.columns)}")

    logger.info(f"Using join column: {maf_join_col}")

    logger.info(f"Reading annotation table: {annotation_table}")
    annotation_df = _read_file_auto(annotation_table, sep='\t', low_memory=False)
    logger.info(f"Annotation table loaded: {annotation_df.shape[0]} rows, {annotation_df.shape[1]} columns")

    if 'GENE_SYMBOL' not in annotation_df.columns:
        logger.error(
            f"'GENE_SYMBOL' column not found in annotation table. Available columns: {list(annotation_df.columns)}")
        raise ValueError(
            f"'GENE_SYMBOL' column not found in annotation table. Available columns: {list(annotation_df.columns)}")

    synonyms_dict = _create_synonyms_dict(annotation_df, synonyms_column)
    maf_df_mapped = _apply_synonyms_mapping(maf_df, maf_join_col, synonyms_dict)

    logger.info("Performing annotation merge...")

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
        # Always use empty string for COSMIC_TIER to ensure Is_Oncogene_any calculates correctly
        if col == 'COSMIC_TIER':
            result_df[col] = result_df[col].fillna("")
        elif pd.api.types.is_numeric_dtype(result_df[col]):
            result_df[col] = result_df[col].fillna(pd.NA)
        else:
            result_df[col] = result_df[col].fillna("")

    logger.info(f"COSMIC annotation completed: {result_df.shape[0]} rows, {result_df.shape[1]} columns")
    logger.info(f"Added {len(cosmic_columns)} COSMIC annotation columns")

    # Process OncoKB annotation if provided
    if oncokb_table is not None:
        logger.info(f"Reading OncoKB table: {oncokb_table}")
        oncokb_df = _read_file_auto(oncokb_table, sep='\t', low_memory=False)
        logger.info(f"OncoKB table loaded: {oncokb_df.shape[0]} rows, {oncokb_df.shape[1]} columns")

        if 'Hugo Symbol' not in oncokb_df.columns:
            logger.error(
                f"'Hugo Symbol' column not found in OncoKB table. Available columns: {list(oncokb_df.columns)}")
            raise ValueError(
                f"'Hugo Symbol' column not found in OncoKB table. Available columns: {list(oncokb_df.columns)}")

        oncokb_synonyms_dict = _create_oncokb_synonyms_dict(oncokb_df, oncokb_synonyms_column)
        result_df_oncokb_mapped = _apply_synonyms_mapping(result_df, maf_join_col, oncokb_synonyms_dict)

        logger.info("Performing OncoKB annotation merge...")

        oncokb_cols_to_rename = {col: f"OncoKB_{col}" for col in oncokb_df.columns if col != 'Hugo Symbol'}
        oncokb_df_renamed = oncokb_df.rename(columns=oncokb_cols_to_rename)

        result_df = result_df_oncokb_mapped.merge(
            oncokb_df_renamed,
            left_on='_mapped_gene_symbol',
            right_on='Hugo Symbol',
            how='left'
        )

        columns_to_drop = []
        if 'Hugo Symbol' in result_df.columns:
            columns_to_drop.append('Hugo Symbol')
        if '_mapped_gene_symbol' in result_df.columns:
            columns_to_drop.append('_mapped_gene_symbol')

        if columns_to_drop:
            result_df = result_df.drop(columns_to_drop, axis=1)

        oncokb_columns = [col for col in result_df.columns if col.startswith('OncoKB_')]
        for col in oncokb_columns:
            if pd.api.types.is_numeric_dtype(result_df[col]):
                result_df[col] = result_df[col].fillna(pd.NA)
            else:
                result_df[col] = result_df[col].fillna("")

        logger.info(f"OncoKB annotation completed: {result_df.shape[0]} rows, {result_df.shape[1]} columns")
        logger.info(f"Added {len(oncokb_columns)} OncoKB annotation columns")

    logger.info(f"Total annotation completed: {result_df.shape[0]} rows, {result_df.shape[1]} columns")
    logger.info("Pandas annotation completed successfully")

    return result_df


def _annotate_with_duckdb(
        data: pd.DataFrame,
        annotation_table: Path,
        join_column: str,
        synonyms_column: str,
        oncokb_table: Optional[Path] = None,
        oncokb_synonyms_column: str = "Gene Aliases"
) -> pd.DataFrame:
    """
    Annotate using DuckDB for large files optimization.
    """
    logger.info(f"Starting DuckDB annotation for DataFrame: {data.shape[0]} rows, {data.shape[1]} columns")

    maf_df = data.copy()

    maf_join_col = find_alias(maf_df.columns, join_column)
    if maf_join_col is None:
        logger.error(f"Join column '{join_column}' not found in data file. Available columns: {list(maf_df.columns)}")
        raise ValueError(
            f"Join column '{join_column}' not found in data file. Available columns: {list(maf_df.columns)}")

    logger.info(f"Using join column: {maf_join_col}")

    logger.info(f"Reading annotation table: {annotation_table}")
    annotation_df = _read_file_auto(annotation_table, sep='\t', low_memory=False)
    logger.info(f"Annotation table loaded: {annotation_df.shape[0]} rows, {annotation_df.shape[1]} columns")

    if 'GENE_SYMBOL' not in annotation_df.columns:
        logger.error(
            f"'GENE_SYMBOL' column not found in annotation table. Available columns: {list(annotation_df.columns)}")
        raise ValueError(
            f"'GENE_SYMBOL' column not found in annotation table. Available columns: {list(annotation_df.columns)}")

    synonyms_dict = _create_synonyms_dict(annotation_df, synonyms_column)
    maf_df_mapped = _apply_synonyms_mapping(maf_df, maf_join_col, synonyms_dict)

    logger.info("Performing optimized merge with DuckDB...")

    conn = duckdb.connect()

    try:
        conn.register('maf_data', maf_df_mapped)
        conn.register('annotation_data', annotation_df)

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
            # Always use empty string for COSMIC_TIER to ensure Is_Oncogene_any calculates correctly
            if col == 'COSMIC_TIER':
                result_df[col] = result_df[col].fillna("")
            elif pd.api.types.is_numeric_dtype(result_df[col]):
                result_df[col] = result_df[col].fillna(pd.NA)
            else:
                result_df[col] = result_df[col].fillna("")

        logger.info(f"COSMIC annotation completed: {result_df.shape[0]} rows, {result_df.shape[1]} columns")
        logger.info(f"Added {len(cosmic_columns)} COSMIC annotation columns")

        # Process OncoKB annotation if provided
        if oncokb_table is not None:
            logger.info(f"Reading OncoKB table: {oncokb_table}")
            oncokb_df = _read_file_auto(oncokb_table, sep='\t', low_memory=False)
            logger.info(f"OncoKB table loaded: {oncokb_df.shape[0]} rows, {oncokb_df.shape[1]} columns")

            if 'Hugo Symbol' not in oncokb_df.columns:
                logger.error(
                    f"'Hugo Symbol' column not found in OncoKB table. Available columns: {list(oncokb_df.columns)}")
                raise ValueError(
                    f"'Hugo Symbol' column not found in OncoKB table. Available columns: {list(oncokb_df.columns)}")

            oncokb_synonyms_dict = _create_oncokb_synonyms_dict(oncokb_df, oncokb_synonyms_column)
            result_df_oncokb_mapped = _apply_synonyms_mapping(result_df, maf_join_col, oncokb_synonyms_dict)

            logger.info("Performing OncoKB annotation merge with DuckDB...")

            conn.register('oncokb_data', oncokb_df)
            conn.register('result_data', result_df_oncokb_mapped)

            oncokb_columns = [col for col in oncokb_df.columns if col != 'Hugo Symbol']
            oncokb_select_clauses = [f"o.{col} as OncoKB_{col}" for col in oncokb_columns]
            oncokb_select_str = ",\n        ".join(oncokb_select_clauses)

            oncokb_merge_query = f"""
            SELECT 
                r.*,
                {oncokb_select_str}
            FROM result_data r
            LEFT JOIN oncokb_data o ON r._mapped_gene_symbol = o."Hugo Symbol"
            """

            result_df = conn.execute(oncokb_merge_query).df()

            if '_mapped_gene_symbol' in result_df.columns:
                result_df = result_df.drop('_mapped_gene_symbol', axis=1)

            oncokb_columns = [col for col in result_df.columns if col.startswith('OncoKB_')]
            for col in oncokb_columns:
                if pd.api.types.is_numeric_dtype(result_df[col]):
                    result_df[col] = result_df[col].fillna(pd.NA)
                else:
                    result_df[col] = result_df[col].fillna("")

            logger.info(f"OncoKB annotation completed: {result_df.shape[0]} rows, {result_df.shape[1]} columns")
            logger.info(f"Added {len(oncokb_columns)} OncoKB annotation columns")

        logger.info(f"Total annotation completed: {result_df.shape[0]} rows, {result_df.shape[1]} columns")
        logger.info("DuckDB annotation completed successfully")

        return result_df

    finally:
        conn.close()


class CancerAnnotateMixin:
    def knownCancer(self, annotation_table, output_path=None, compress_output=True, join_column="Hugo_Symbol",
                    oncokb_table=None, in_place=False):
        """
        Annotate mutations with COSMIC and OncoKB cancer-related annotations.
        
        Parameters
        ----------
        annotation_table : str | Path
            Path to the COSMIC annotation table (.tsv or .tsv.gz)
        output_path : str | Path, optional
            Output file path. If provided, saves the annotated DataFrame to file
        compress_output : bool, default True
            Whether to compress the output file with gzip
        join_column : str, default "Hugo_Symbol"
            Column name to use for joining (canonical name from fields.py)
        oncokb_table : str | Path, optional
            Path to the OncoKB cancer gene list table (.tsv). If provided, OncoKB
            annotations will be added to the output.
        in_place : bool, default False
            If True, replaces self.data with annotated data. If False, returns
            annotated DataFrame for external use.
            
        Returns
        -------
        pd.DataFrame or None
            If in_place=False, returns annotated DataFrame.
            If in_place=True, returns None and updates self.data.

        Raises
        ------
        FileNotFoundError
            If annotation files don't exist
        ValueError
            If join column is not found in DataFrame
        """
        from datetime import datetime

        data_memory_mb = self.data.memory_usage(deep=True).sum() / (1024 * 1024)
        data_memory_gb = data_memory_mb / 1024
        # use_duckdb = data_memory_gb > 2.0  # DuckDB option disabled

        logger.info(f"DataFrame memory usage: {data_memory_gb:.2f} GB")
        logger.info("Using pandas backend for annotation")

        # Get full annotations using pandas backend
        # if use_duckdb:  # DuckDB option disabled
        #     full_df = _annotate_with_duckdb(
        #         data=self.data,
        #         annotation_table=Path(annotation_table),
        #         join_column=join_column,
        #         synonyms_column="SYNONYMS",
        #         oncokb_table=Path(oncokb_table) if oncokb_table else None
        #     )
        # else:
        full_df = _annotate_with_pandas(
            data=self.data,
            annotation_table=Path(annotation_table),
            join_column=join_column,
            synonyms_column="SYNONYMS",
            oncokb_table=Path(oncokb_table) if oncokb_table else None
        )

        # Define the specific columns we want to keep
        target_columns = [
            "COSMIC_ROLE_IN_CANCER",
            "COSMIC_TIER",
            "OncoKB_Is Oncogene",
            "OncoKB_Is Tumor Suppressor Gene",
            "OncoKB_OncoKB Annotated",
            "OncoKB_MSK-IMPACT",
            "OncoKB_MSK-HEME",
            "OncoKB_FOUNDATION ONE",
            "OncoKB_FOUNDATION ONE HEME",
            "OncoKB_Vogelstein"
        ]

        # Get original data columns
        original_columns = [col for col in full_df.columns if
                            not (col.startswith('COSMIC_') or col.startswith('OncoKB_'))]

        # Filter to keep only original columns plus target annotation columns
        available_target_columns = [col for col in target_columns if col in full_df.columns]
        columns_to_keep = original_columns + available_target_columns

        filtered_df = full_df[columns_to_keep].copy()

        # Is_Oncogene_any -> True if OncoKB_Is Oncogene is True OR COSMIC has any annotation (not empty)
        oncokb_oncogene = filtered_df.get('OncoKB_Is Oncogene', pd.Series([False] * len(filtered_df)))
        cosmic_role = filtered_df.get('COSMIC_ROLE_IN_CANCER', pd.Series([''] * len(filtered_df)))
        cosmic_tier = filtered_df.get('COSMIC_TIER', pd.Series([''] * len(filtered_df)))

        # Convert OncoKB_Is Oncogene to boolean (handle string representations)
        oncokb_is_oncogene_bool = oncokb_oncogene.astype(str).str.lower().isin(['true', '1', 'yes'])

        # COSMIC annotation exists if ROLE_IN_CANCER or TIER is not empty
        # Handle pd.NA values and "<NA>" strings as empty
        def _is_not_empty_or_na(series):
            """Check if series has meaningful values (not empty, not NA, not '<NA>')"""
            return (
                    series.notna() &
                    (series.astype(str).str.strip() != '') &
                    (series.astype(str).str.strip() != '<NA>')
            )

        cosmic_role_has_annotation = _is_not_empty_or_na(cosmic_role)
        cosmic_tier_has_annotation = _is_not_empty_or_na(cosmic_tier)
        cosmic_has_annotation = cosmic_role_has_annotation | cosmic_tier_has_annotation

        # Is_Oncogene_any is True if either source indicates oncogene status
        filtered_df['Is_Oncogene_any'] = oncokb_is_oncogene_bool | cosmic_has_annotation

        # Optional persistence: save to file if output_path is provided
        if output_path is not None:
            final_output_path = Path(output_path)
            if compress_output and not str(final_output_path).endswith('.gz'):
                final_output_path = final_output_path.with_suffix(final_output_path.suffix + '.gz')

            logger.info(f"Saving KnownCancer annotated file to: {final_output_path}")
            _write_file_auto(filtered_df, final_output_path, compress_output)
            logger.info(f"Output file saved: {final_output_path}")

        logger.info("KnownCancer annotation completed successfully")
        logger.info(f"Filtered to {len(available_target_columns)} annotation columns plus Is_Oncogene_any field")

        # Update metadata
        if self.metadata is not None:
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            annotation_note = f"COSMIC cancer annotation applied at {timestamp}"
            if oncokb_table:
                annotation_note += " (with OncoKB data)"

            if self.metadata.notes:
                self.metadata.notes += f"; {annotation_note}"
            else:
                self.metadata.notes = annotation_note

        # Handle in_place flag
        if in_place:
            self.data = filtered_df
            return None
        else:
            return filtered_df


__all__ = [
    "CancerAnnotateMixin"
]
