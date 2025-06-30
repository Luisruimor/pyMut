import pandas as pd
import logging
from copy import deepcopy
from typing import Optional, List, Union
from ..utils.format import formatear_chr
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
    
    Parameters
    ----------
    chrom : str, list of str, or None, optional
        Chromosome(s) to filter by (e.g., 'chr1', '1', 'X', 'Y', ['chr1', 'chr2']).
        If None, no chromosome filtering is applied.
    sample : str, list of str, or None, optional
        Sample(s) to filter by (e.g., 'TCGA-AB-2802', ['TCGA-AB-2802', 'TCGA-AB-2803']).
        If None, no sample filtering is applied.
    sample_column : str, optional
        Name of the column containing sample information. 
        Defaults to 'Tumor_Sample_Barcode'.
    
    Returns
    -------
    PyMutation
        A new instance of PyMutation containing only the rows that match
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
        chrom_formatted = [formatear_chr(str(c)) for c in chrom_list]
        logger.info(f"Chromosomes to filter: {chrom_formatted}")
        
        # Format CHROM column in DataFrame
        df["CHROM"] = df["CHROM"].astype(str).map(formatear_chr)
        
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
        
        # Apply sample filter
        sample_mask = df[sample_column].isin(sample_list)
        df = df[sample_mask]
        
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
    
    # Return new PyMutation instance with filtered data and updated metadata
    return type(self)(data=df, metadata=updated_metadata, samples=self.samples)