import logging
import math
import time

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Configure logger
logger = logging.getLogger(__name__)

# OncoKB API constants
ONCOKB_ENDPOINT = "https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange"
VALID_REFERENCE_GENOMES = ["GRCh37", "GRCh38"]


class ActionableMutationMixin:
    def actionable_mutations_oncokb(self, token: str, batch_size: int = 5000, timeout: int = 30,
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
        referenceGenome = f"GRCh{self.metadata.assembly}"
        logger.info("Using reference genome: %s", referenceGenome)

        # Validate reference genome
        if referenceGenome not in VALID_REFERENCE_GENOMES:
            raise ValueError(f"Invalid reference genome: {referenceGenome}. Must be one of {VALID_REFERENCE_GENOMES}")

        # Extract required columns directly from PyMutation.data
        oncokb_input_df = pd.DataFrame()

        if "CHROM" not in self.data.columns:
            raise ValueError("Missing required column 'CHROM' for OncoKB input")
        oncokb_input_df["CHROM"] = self.data["CHROM"].str.lstrip("chr")

        if "POS" not in self.data.columns:
            raise ValueError("Missing required column 'POS' for OncoKB input")
        # Convert to int32 for memory efficiency
        oncokb_input_df["POS"] = self.data["POS"].astype('int32')

        if "REF" not in self.data.columns:
            raise ValueError("Missing required column 'REF' for OncoKB input")
        oncokb_input_df["REF"] = self.data["REF"]

        # Calculate End_Position as POS + len(REF) and convert to int32
        oncokb_input_df["END"] = (oncokb_input_df["POS"] + oncokb_input_df["REF"].str.len()).astype('int32')

        if "ALT" in self.data.columns:
            oncokb_input_df["ALT"] = self.data["ALT"]
        elif "ALT" in self.data.columns:
            logger.info("Using ALT column for ALT")
            oncokb_input_df["ALT"] = self.data["ALT"]
        else:
            raise ValueError("Missing required column 'ALT' for OncoKB input")

        # Check for nulls in required columns
        required_columns = ["CHROM", "POS", "END", "REF", "ALT"]
        initial_row_count = len(oncokb_input_df)

        # Create a mask for rows with null values in any required column
        null_mask = pd.Series(False, index=oncokb_input_df.index)
        for col_name in required_columns:
            col_nulls = oncokb_input_df[col_name].isnull()
            if col_nulls.any():
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
            logger.info("Processing batch %d/%d with %d variants", i + 1, num_batches, len(batch_df))

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

                payload.append(variant_entry)

            # Send API request with retries
            retry_count = 0
            max_backoff = 30

            while True:
                try:
                    logger.info("Sending batch %d/%d to OncoKB API", i + 1, num_batches)
                    response = session.post(
                        ONCOKB_ENDPOINT,
                        headers=headers,
                        json=payload,
                        timeout=timeout
                    )

                    if response.status_code == 200:
                        # Success - parse the response
                        resp_json = response.json()
                        logger.info("Batch %d/%d processed successfully", i + 1, num_batches)

                        for j, anno in enumerate(resp_json):
                            if anno.get("geneExist") is False:
                                logger.warning("Gene does not exist for variant at index %d in batch %d", j, i + 1)
                            if anno.get("oncogenic") == "Unknown":
                                logger.warning("Unknown oncogenicity for variant at index %d in batch %d", j, i + 1)

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
                            logger.warning("Empty response for batch %d/%d", i + 1, num_batches)

                        break

                    elif response.status_code == 401:
                        raise ValueError(
                            f"Authentication error: Invalid OncoKB token. Status code: {response.status_code}")

                    elif response.status_code in [429, 500, 502, 503, 504]:
                        retry_count += 1
                        if retry_count > max_retries:
                            raise ValueError(
                                f"Maximum retries exceeded for batch {i + 1}. Last status code: {response.status_code}")

                        backoff_time = min(max_backoff, retry_backoff * (2 ** (retry_count - 1)))
                        logger.warning(
                            "Received status code %d for batch %d/%d. Retrying in %.1f seconds (retry %d/%d)",
                            response.status_code, i + 1, num_batches, backoff_time, retry_count, max_retries)
                        time.sleep(backoff_time)
                        continue

                    else:
                        raise ValueError(
                            f"Error from OncoKB API for batch {i + 1}. Status code: {response.status_code}, Response: {response.text}")

                except requests.exceptions.RequestException as e:
                    retry_count += 1
                    if retry_count > max_retries:
                        raise ValueError(f"Maximum retries exceeded for batch {i + 1} due to network error: {str(e)}")

                    backoff_time = min(max_backoff, retry_backoff * (2 ** (retry_count - 1)))
                    logger.warning("Network error for batch %d/%d: %s. Retrying in %.1f seconds (retry %d/%d)",
                                   i + 1, num_batches, str(e), backoff_time, retry_count, max_retries)
                    time.sleep(backoff_time)
                    continue

            # Sleep between batches to avoid rate limiting
            if i < num_batches - 1:  # Don't sleep after the last batch
                logger.info("Sleeping for 1 second before processing next batch")
                time.sleep(1)

        # Add annotations to self.data
        if annotations_dict:
            logger.info("Adding OncoKB annotations to self.data")

            for field in oncokb_fields:
                column_name = f"oncokb_{field}"
                if column_name not in self.data.columns:
                    self.data[column_name] = None

            for orig_idx, annotations in annotations_dict.items():
                for field, value in annotations.items():
                    self.data.at[orig_idx, field] = value

            logger.info("OncoKB annotation completed successfully: %d variants annotated", len(annotations_dict))
        else:
            logger.warning("No annotations returned from OncoKB API")

        return self.data


__all__ = [
    "ActionableMutationMixin"
]
