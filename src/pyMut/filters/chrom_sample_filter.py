import logging
from copy import deepcopy
from typing import Optional, List, Union

from ..utils.constants import SAMPLE_COLUMN
from ..utils.format import format_chr

# Logging configuration
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Change to DEBUG for more verbosity
if not logger.handlers:
    logging.basicConfig(
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        level=logging.INFO,
    )


class ChromSampleFilterMixin:
    """
    Mixin class providing chromosome and sample filtering functionality for PyMutation objects.
    
    This mixin adds the ability to filter by chromosome and/or sample,
    following the same architectural pattern as other mixins in the project.
    """

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
        if chrom is None and sample is None:
            logger.error("At least one of 'chrom' or 'sample' must be provided")
            raise ValueError("At least one of 'chrom' or 'sample' must be provided")

        df = self.data.copy()
        filter_descriptions = []

        if chrom is not None:
            if "CHROM" not in df.columns:
                logger.error("Column 'CHROM' does not exist in the DataFrame")
                raise KeyError("Column 'CHROM' does not exist in the DataFrame")

            if isinstance(chrom, str):
                chrom_list = [chrom]
            else:
                chrom_list = list(chrom)

            chrom_formatted = [format_chr(str(c)) for c in chrom_list]
            logger.info(f"Chromosomes to filter: {chrom_formatted}")

            df["CHROM"] = df["CHROM"].astype(str).map(format_chr)
            chrom_mask = df["CHROM"].isin(chrom_formatted)
            df = df[chrom_mask]

            chrom_desc = ",".join(chrom_formatted)
            filter_descriptions.append(f"chromosome:{chrom_desc}")
            logger.info(f"Applied chromosome filter: {chrom_desc}")

        if sample is not None:
            if isinstance(sample, str):
                sample_list = [sample]
            else:
                sample_list = list(sample)

            logger.info(f"Samples to filter: {sample_list}")

            # Check if MAF-style (sample column) or VCF-style (individual sample columns)
            has_sample_column = sample_column in df.columns

            if has_sample_column:
                logger.info(f"Using MAF-style filtering with column '{sample_column}'")
                sample_mask = df[sample_column].isin(sample_list)
                df = df[sample_mask]
            else:
                logger.info("Using VCF-style filtering (no sample column found)")
                all_samples = self.samples if hasattr(self, 'samples') and self.samples else []
                missing_samples = [s for s in sample_list if s not in all_samples]
                if missing_samples:
                    logger.error(f"Requested samples not found in PyMutation.samples: {missing_samples}")
                    raise ValueError(f"Requested samples not found in PyMutation.samples: {missing_samples}")
                logger.info("VCF-style: keeping all variant rows, filtering only sample columns")

            # Column filtering: keep VCF-like columns and requested sample columns
            vcf_like_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
            all_samples = self.samples if hasattr(self, 'samples') and self.samples else []
            sample_columns_to_keep = [s for s in sample_list if s in df.columns]
            sample_columns_to_remove = [s for s in all_samples if s in df.columns and s not in sample_list]

            non_sample_cols = [col for col in df.columns if col not in vcf_like_cols + all_samples]
            if has_sample_column and sample_column not in vcf_like_cols:
                non_sample_cols.append(sample_column)

            columns_to_keep = vcf_like_cols + sample_columns_to_keep + non_sample_cols
            columns_to_keep = [col for col in columns_to_keep if col in df.columns]
            df = df[columns_to_keep]

            logger.info(f"Sample columns kept: {sample_columns_to_keep}")
            logger.info(f"Sample columns removed: {sample_columns_to_remove}")

            sample_desc = ",".join(sample_list)
            filter_descriptions.append(f"sample:{sample_desc}")
            logger.info(f"Applied sample filter: {sample_desc}")

        updated_metadata = deepcopy(self.metadata)
        combined_filter_description = "|".join(filter_descriptions)

        if hasattr(updated_metadata, 'filters') and updated_metadata.filters:
            updated_metadata.filters = updated_metadata.filters + [combined_filter_description]
        else:
            updated_metadata.filters = [combined_filter_description]

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

        if sample is not None:
            filtered_samples = [s for s in sample_list if s in df.columns]
        else:
            filtered_samples = self.samples

        return type(self)(data=df, metadata=updated_metadata, samples=filtered_samples)
