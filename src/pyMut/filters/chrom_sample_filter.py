import pandas as pd
import logging
from copy import deepcopy
from typing import Optional, List, Union
from ..utils.format import format_chr
from ..utils.constants import SAMPLE_COLUMN

# ────────────────────────────────────────────────────────────────
# Logging configuration
# ────────────────────────────────────────────────────────────────
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Change to DEBUG for more verbosity
if not logger.handlers:
    logging.basicConfig(
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        level=logging.INFO,
    )


def filter_by_chrom_sample(
    self, 
    chrom: Optional[Union[str, List[str]]] = None, 
    sample: Optional[Union[str, List[str]]] = None,
    sample_column: str = SAMPLE_COLUMN
):
    """
    Filter this PyMutation by chromosome and/or sample.

    This method allows filtering by chromosome, sample, or both. At least one
    parameter must be provided. Both parameters accept single values or lists
    of values for multiple filtering.

    When filtering by sample, this method performs both row and column filtering:
    - Row filtering: Keeps only rows where the sample column matches the specified sample(s)
    - Column filtering: Removes columns corresponding to samples not in the filter list
      (preserves VCF-like columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER and other non-sample columns)

    Parameters
    ----------
    chrom : str, list of str, or None, optional
        Chromosome(s) to filter by (e.g., 'chr1', '1', 'X', 'Y', ['chr1', 'chr2']).
        If None, no chromosome filtering is applied.
    sample : str, list of str, or None, optional
        Sample(s) to filter by (e.g., 'TCGA-AB-2802', ['TCGA-AB-2802', 'TCGA-AB-2803']).
        If None, no sample filtering is applied.
        When provided, both rows and columns are filtered to include only the specified samples.
    sample_column : str, optional
        Name of the column containing sample information. 
        Defaults to 'Tumor_Sample_Barcode'.

    Returns
    -------
    PyMutation
        A new instance of PyMutation containing only the rows and columns that match
        the specified chromosome and/or sample criteria. The metadata is 
        copied and updated to record the applied filter.

    Raises
    ------
    ValueError
        If both chrom and sample are None, or if required columns are missing.
    KeyError
        If the DataFrame does not contain the required 'CHROM' column when
        chromosome filtering is requested, or the sample column when sample
        filtering is requested.

    Examples
    --------
    >>> # Filter by chromosome only (preserves all columns)
    >>> filtered = py_mutation.filter_by_chrom_sample(chrom='chr17')

    >>> # Filter by sample (removes both rows and sample columns)
    >>> filtered = py_mutation.filter_by_chrom_sample(sample='TCGA-AB-2988')

    >>> # Filter by multiple samples
    >>> filtered = py_mutation.filter_by_chrom_sample(sample=['TCGA-AB-2988', 'TCGA-AB-2869'])

    >>> # Combined filtering
    >>> filtered = py_mutation.filter_by_chrom_sample(chrom='chr17', sample='TCGA-AB-2988')
    """
    # Validate that at least one filter parameter is provided
    if chrom is None and sample is None:
        logger.error("At least one of 'chrom' or 'sample' must be provided")
        raise ValueError("At least one of 'chrom' or 'sample' must be provided")

    df = self.data.copy()
    filter_descriptions = []

    # Apply chromosome filter if provided
    if chrom is not None:
        # Verify CHROM column exists
        if "CHROM" not in df.columns:
            logger.error("Column 'CHROM' does not exist in the DataFrame")
            raise KeyError("Column 'CHROM' does not exist in the DataFrame")

        # Convert single chromosome to list for uniform processing
        if isinstance(chrom, str):
            chrom_list = [chrom]
        else:
            chrom_list = list(chrom)

        # Format chromosomes for consistency
        chrom_formatted = [format_chr(str(c)) for c in chrom_list]
        logger.info(f"Chromosomes to filter: {chrom_formatted}")

        # Format CHROM column in DataFrame
        df["CHROM"] = df["CHROM"].astype(str).map(format_chr)

        # Apply chromosome filter
        chrom_mask = df["CHROM"].isin(chrom_formatted)
        df = df[chrom_mask]

        # Add to filter description
        chrom_desc = ",".join(chrom_formatted)
        filter_descriptions.append(f"chromosome:{chrom_desc}")
        logger.info(f"Applied chromosome filter: {chrom_desc}")

    # Apply sample filter if provided
    if sample is not None:
        # Verify sample column exists
        if sample_column not in df.columns:
            logger.error(f"Column '{sample_column}' does not exist in the DataFrame")
            raise KeyError(f"Column '{sample_column}' does not exist in the DataFrame")

        # Convert single sample to list for uniform processing
        if isinstance(sample, str):
            sample_list = [sample]
        else:
            sample_list = list(sample)

        logger.info(f"Samples to filter: {sample_list}")

        # Apply sample filter (rows)
        sample_mask = df[sample_column].isin(sample_list)
        df = df[sample_mask]

        # Apply sample filter (columns) - Remove columns for samples not in the filter
        # Define VCF-like columns that should always be kept
        vcf_like_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]

        # Get all available samples from the PyMutation object
        all_samples = self.samples if hasattr(self, 'samples') and self.samples else []

        # Identify which columns are sample columns vs other columns
        sample_columns_to_keep = [s for s in sample_list if s in df.columns]
        sample_columns_to_remove = [s for s in all_samples if s in df.columns and s not in sample_list]

        # Keep VCF-like columns, filtered sample columns, and any other non-sample columns
        columns_to_keep = (
            vcf_like_cols + 
            sample_columns_to_keep + 
            [col for col in df.columns if col not in vcf_like_cols + all_samples]
        )

        # Filter columns to only those that actually exist in the DataFrame
        columns_to_keep = [col for col in columns_to_keep if col in df.columns]

        # Apply column filtering
        df = df[columns_to_keep]

        logger.info(f"Sample columns kept: {sample_columns_to_keep}")
        logger.info(f"Sample columns removed: {sample_columns_to_remove}")

        # Add to filter description
        sample_desc = ",".join(sample_list)
        filter_descriptions.append(f"sample:{sample_desc}")
        logger.info(f"Applied sample filter: {sample_desc}")

    # Create updated metadata with information about the applied filters
    updated_metadata = deepcopy(self.metadata)

    # Create combined filter description
    combined_filter_description = "|".join(filter_descriptions)

    # Add the filter to the existing filters list
    if hasattr(updated_metadata, 'filters') and updated_metadata.filters:
        updated_metadata.filters = updated_metadata.filters + [combined_filter_description]
    else:
        updated_metadata.filters = [combined_filter_description]

    # Log filtering information
    original_count = len(self.data)
    filtered_count = len(df)

    logger.info(f"Combined filter applied: {combined_filter_description}")
    logger.info(f"Variants before filter: {original_count}")
    logger.info(f"Variants after filter: {filtered_count}")
    logger.info(f"Variants filtered out: {original_count - filtered_count}")

    if filtered_count == 0:
        logger.warning(f"No variants found matching the filter criteria: {combined_filter_description}")
    elif filtered_count == original_count:
        logger.warning(f"Filter did not remove any variants - check filter criteria: {combined_filter_description}")
    else:
        logger.info(f"Successfully applied filter: {combined_filter_description}")

    # Update samples list if sample filtering was applied
    if sample is not None:
        # Only keep samples that were in the filter and exist in the data
        filtered_samples = [s for s in sample_list if s in df.columns]
    else:
        # If no sample filtering, keep original samples
        filtered_samples = self.samples

    # Return new PyMutation instance with filtered data and updated metadata
    return type(self)(data=df, metadata=updated_metadata, samples=filtered_samples)
