import logging
import math
import time
import requests
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Union, Optional, List, Dict, Any
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from ..utils.fields import col
from ..core import PyMutation

# Configure logger
logger = logging.getLogger(__name__)

# OncoKB API constants
ONCOKB_ENDPOINT = "https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange"
VALID_REFERENCE_GENOMES = ["GRCh37", "GRCh38"]

def export_oncokb_input(self, token: str, batch_size: int = 5000, timeout: int = 30, 
                      max_retries: int = 3, retry_backoff: float = 1.0) -> pd.DataFrame:
    """
    Export mutation data to OncoKB API and add annotations to self.data.
    
    This method extracts the required columns from self.data, sends the data to the OncoKB API
    in batches, and adds the annotations as columns to self.data.
    
    Parameters
    ----------
    token : str
        OncoKB API authentication token
    batch_size : int, optional
        Maximum number of variants to send in a single API request (default: 5000)
    timeout : int, optional
        Timeout for API requests in seconds (default: 30)
    max_retries : int, optional
        Maximum number of retries for failed API requests (default: 3)
    retry_backoff : float, optional
        Backoff factor for retries (default: 1.0)
        
    Returns
    -------
    pd.DataFrame
        The original self.data DataFrame with OncoKB annotations added as columns.
        
    Raises
    ------
    ValueError
        If the DataFrame doesn't contain the necessary data for export or if the reference genome is invalid.
    requests.exceptions.RequestException
        If there's an error with the API request that can't be resolved with retries.
    """
    # Generate referenceGenome variable
    referenceGenome = f"GRCh{self.metadata.assembly}"
    logger.info("Using reference genome: %s", referenceGenome)

    # Validate reference genome
    if referenceGenome not in VALID_REFERENCE_GENOMES:
        raise ValueError(f"Invalid reference genome: {referenceGenome}. Must be one of {VALID_REFERENCE_GENOMES}")

    # Extract required columns using the col function to handle aliases
    oncokb_input_df = pd.DataFrame()
    
    # Extract Chromosome -> CHROM
    chromosome = col(self.data, "Chromosome", required=False)
    if chromosome is None:
        raise ValueError("Missing required column 'Chromosome' for OncoKB input")
    # Remove 'chr' prefix using str.lstrip
    oncokb_input_df["CHROM"] = chromosome.str.lstrip("chr")
    
    # Extract Start_Position -> POS
    start_position = col(self.data, "Start_Position", required=False)
    if start_position is None:
        raise ValueError("Missing required column 'Start_Position' for OncoKB input")
    # Convert to int32 for memory efficiency
    oncokb_input_df["POS"] = start_position.astype('int32')
    
    # Extract Reference_Allele -> REF
    reference_allele = col(self.data, "Reference_Allele", required=False)
    if reference_allele is None:
        raise ValueError("Missing required column 'Reference_Allele' for OncoKB input")
    oncokb_input_df["REF"] = reference_allele
    
    # Calculate End_Position as POS + len(REF) and convert to int32
    oncokb_input_df["END"] = (oncokb_input_df["POS"] + oncokb_input_df["REF"].str.len()).astype('int32')
    
    # Extract Tumor_Seq_Allele2 -> ALT
    tumor_seq_allele2 = col(self.data, "Tumor_Seq_Allele2", required=False)
    if tumor_seq_allele2 is None:
        # Try to use ALT column if Tumor_Seq_Allele2 is not present
        alt = col(self.data, "ALT", required=False)
        if alt is not None:
            logger.info("Using ALT column for Tumor_Seq_Allele2")
            oncokb_input_df["ALT"] = alt
        else:
            raise ValueError("Missing required column 'Tumor_Seq_Allele2' for OncoKB input")
    else:
        oncokb_input_df["ALT"] = tumor_seq_allele2
    
    # Validate the DataFrame
    # Check for nulls in required columns
    required_columns = ["CHROM", "POS", "END", "REF", "ALT"]
    initial_row_count = len(oncokb_input_df)

    # Create a mask for rows with null values in any required column
    null_mask = pd.Series(False, index=oncokb_input_df.index)
    for col_name in required_columns:
        col_nulls = oncokb_input_df[col_name].isnull()
        if col_nulls.any():
            # Log which column has nulls and how many
            null_count = col_nulls.sum()
            logger.warning("Column '%s' contains %d null values, these rows will be removed", col_name, null_count)
            null_mask = null_mask | col_nulls

    # Remove rows with null values
    if null_mask.any():
        oncokb_input_df = oncokb_input_df[~null_mask]
        removed_count = null_mask.sum()
        logger.warning("Removed %d rows with null values from input data", removed_count)
        logger.info("Remaining rows after null removal: %d", len(oncokb_input_df))
    
    # Add index to preserve original order
    oncokb_input_df = oncokb_input_df.reset_index().rename(columns={"index": "_original_index"})
    
    # Check if TUMOR_TYPE column exists in the original data
    has_tumor_type = "TUMOR_TYPE" in self.data.columns
    if has_tumor_type:
        oncokb_input_df["TUMOR_TYPE"] = self.data["TUMOR_TYPE"].reset_index(drop=True)
        logger.info("TUMOR_TYPE column found and will be included in the OncoKB request")
    
    # Split data into batches
    num_variants = len(oncokb_input_df)
    num_batches = math.ceil(num_variants / batch_size)
    logger.info("Splitting %d variants into %d batches of max %d variants each", 
                num_variants, num_batches, batch_size)
    
    batches = np.array_split(oncokb_input_df, num_batches)
    
    # Set up API request headers
    headers = {
        "Authorization": f"Bearer {token}",
        "Content-Type": "application/json"
    }
    
    # Set up retry strategy
    retry_strategy = Retry(
        total=max_retries,
        backoff_factor=retry_backoff,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["POST"]
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session = requests.Session()
    session.mount("https://", adapter)
    
    # Define the fields to extract from the OncoKB API response
    oncokb_fields = [
        "highestSensitiveLevel",
        "highestResistanceLevel",
        "highestDiagnosticImplicationLevel",
        "highestPrognosticImplicationLevel",
        "otherSignificantSensitiveLevels",
        "otherSignificantResistanceLevels",
        "hotspot",
        "geneSummary",
        "variantSummary",
        "tumorTypeSummary",
        "prognosticSummary",
        "diagnosticSummary",
        "diagnosticImplications",
        "prognosticImplications",
        "treatments",
        "dataVersion",
        "lastUpdate",
        "vus"
    ]
    
    # Create a dictionary to store annotations for each variant
    annotations_dict = {}
    
    for i, batch_df in enumerate(batches):
        logger.info("Processing batch %d/%d with %d variants", i+1, num_batches, len(batch_df))
        
        # Construct payload for this batch
        payload = []
        for _, row in batch_df.iterrows():
            # Create genomic location string
            loc_string = f"{row.CHROM},{row.POS},{row.END},{row.REF},{row.ALT}"
            
            # Create variant entry
            variant_entry = {
                "referenceGenome": referenceGenome,
                "genomicLocation": loc_string
            }
            
            # Add tumor type if available
            if has_tumor_type and not pd.isna(row.get("TUMOR_TYPE", None)):
                variant_entry["tumorType"] = row.TUMOR_TYPE
                
            # Add Hugo_Symbol if available
            # This can help the OncoKB API identify genes more accurately
            hugo_symbol = col(self.data, "Hugo_Symbol", required=False)
            if hugo_symbol is not None:
                # Get the Hugo_Symbol for this row
                orig_idx = row["_original_index"]
                gene_symbol = hugo_symbol.iloc[orig_idx]
                if not pd.isna(gene_symbol):
                    variant_entry["hugoSymbol"] = gene_symbol
                    logger.info("Added Hugo_Symbol '%s' to variant entry", gene_symbol)
                
            payload.append(variant_entry)
        
        # Send API request with retries
        retry_count = 0
        max_backoff = 30  # Maximum backoff time in seconds
        
        while True:
            try:
                logger.info("Sending batch %d/%d to OncoKB API", i+1, num_batches)
                response = session.post(
                    ONCOKB_ENDPOINT,
                    headers=headers,
                    json=payload,
                    timeout=timeout
                )
                
                # Handle response based on status code
                if response.status_code == 200:
                    # Success - parse the response
                    resp_json = response.json()
                    logger.info("Batch %d/%d processed successfully", i+1, num_batches)
                    
                    # Check for warnings in the response
                    for j, anno in enumerate(resp_json):
                        if anno.get("geneExist") is False:
                            logger.warning("Gene does not exist for variant at index %d in batch %d", j, i+1)
                        if anno.get("oncogenic") == "Unknown":
                            logger.warning("Unknown oncogenicity for variant at index %d in batch %d", j, i+1)
                    
                    # Process the response and extract the specified fields
                    if resp_json:
                        for j, anno in enumerate(resp_json):
                            # Get the original index for this variant
                            orig_idx = batch_df.iloc[j]["_original_index"]
                            
                            # Extract only the specified fields
                            variant_annotations = {}
                            for field in oncokb_fields:
                                variant_annotations[f"oncokb_{field}"] = anno.get(field)
                            
                            # Store the annotations for this variant
                            annotations_dict[orig_idx] = variant_annotations
                    else:
                        logger.warning("Empty response for batch %d/%d", i+1, num_batches)
                    
                    # Break out of the retry loop
                    break
                    
                elif response.status_code == 401:
                    # Authentication error - abort
                    raise ValueError(f"Authentication error: Invalid OncoKB token. Status code: {response.status_code}")
                    
                elif response.status_code in [429, 500, 502, 503, 504]:
                    # Retry with exponential backoff for these status codes
                    retry_count += 1
                    if retry_count > max_retries:
                        raise ValueError(f"Maximum retries exceeded for batch {i+1}. Last status code: {response.status_code}")
                    
                    # Calculate backoff time
                    backoff_time = min(max_backoff, retry_backoff * (2 ** (retry_count - 1)))
                    logger.warning("Received status code %d for batch %d/%d. Retrying in %.1f seconds (retry %d/%d)",
                                  response.status_code, i+1, num_batches, backoff_time, retry_count, max_retries)
                    time.sleep(backoff_time)
                    continue
                    
                else:
                    # Other error - abort
                    raise ValueError(f"Error from OncoKB API for batch {i+1}. Status code: {response.status_code}, Response: {response.text}")
                    
            except requests.exceptions.RequestException as e:
                # Handle network errors
                retry_count += 1
                if retry_count > max_retries:
                    raise ValueError(f"Maximum retries exceeded for batch {i+1} due to network error: {str(e)}")
                
                # Calculate backoff time
                backoff_time = min(max_backoff, retry_backoff * (2 ** (retry_count - 1)))
                logger.warning("Network error for batch %d/%d: %s. Retrying in %.1f seconds (retry %d/%d)",
                              i+1, num_batches, str(e), backoff_time, retry_count, max_retries)
                time.sleep(backoff_time)
                continue
        
        # Sleep between batches to avoid rate limiting
        if i < num_batches - 1:  # Don't sleep after the last batch
            logger.info("Sleeping for 1 second before processing next batch")
            time.sleep(1)

    # Add annotations to self.data
    if annotations_dict:
        logger.info("Adding OncoKB annotations to self.data")
    
        # Initialize new columns in self.data with None values
        for field in oncokb_fields:
            column_name = f"oncokb_{field}"
            if column_name not in self.data.columns:
                self.data[column_name] = None
    
        # Add annotations to self.data based on the original indices
        for orig_idx, annotations in annotations_dict.items():
            for field, value in annotations.items():
                self.data.at[orig_idx, field] = value
    
        logger.info("OncoKB annotation completed successfully: %d variants annotated", len(annotations_dict))
    else:
        logger.warning("No annotations returned from OncoKB API")

    # Return self.data with the added annotations
    return self.data

# Add the export_oncokb_input method to the PyMutation class
PyMutation.export_oncokb_input = export_oncokb_input