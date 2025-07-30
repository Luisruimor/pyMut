import logging

# Logging configuration
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Change to DEBUG for more verbosity
if not logger.handlers:
    logging.basicConfig(
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        level=logging.INFO,
    )


class PassFilterMixin:
    """
    Mixin class providing PASS filter functionality for PyMutation objects.
    
    This mixin adds the ability to check if specific records have FILTER == "PASS",
    following the same architectural pattern as other mixins in the project.
    """

    def pass_filter(self, chrom: str, pos: int, ref: str, alt: str):
        """
        Filter this PyMutation by checking if a specific record (CHROM, POS, REF, ALT)
        has FILTER == "PASS" using pyarrow for optimized performance.

        Parameters
        ----------
        chrom : str
            Chromosome of the record to check (e.g., 'chr1', '1', 'X', 'Y').
        pos : int
            Position of the record to check.
        ref : str
            Reference allele of the record to check.
        alt : str
            Alternative allele of the record to check.

        Returns
        -------
        bool
            True if the record exists and has FILTER == "PASS", False otherwise.

        Raises
        ------
        KeyError
            If the DataFrame does not contain the required columns.
        """
        df = self.data

        # Verify required columns exist
        required_columns = ["CHROM", "POS", "REF", "ALT", "FILTER"]
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            logger.error(f"Missing required columns: {missing_columns}")
            raise KeyError(f"Missing required columns: {missing_columns}")

        # Format input chromosome to ensure consistency
        from ..utils.format import format_chr
        chrom_formatted = format_chr(str(chrom))
        logger.info(f"Checking PASS filter for: {chrom_formatted}:{pos} {ref}>{alt}")

        # Try to use pyarrow if available
        try:
            logger.info("Attempting to use PyArrow optimization")

            # Convert relevant columns to pyarrow types for optimization
            df_work = df.copy(deep=True)

            # Format CHROM column before optimization
            df_work["CHROM"] = df_work["CHROM"].astype(str).map(format_chr)

            # Optimize columns with pyarrow if not already
            if not str(df_work["CHROM"].dtype).startswith("string"):
                df_work["CHROM"] = df_work["CHROM"].astype("string[pyarrow]")

            if not str(df_work["POS"].dtype).startswith("int") or "pyarrow" not in str(df_work["POS"].dtype):
                df_work["POS"] = df_work["POS"].astype("int64[pyarrow]")

            if not str(df_work["REF"].dtype).startswith("string"):
                df_work["REF"] = df_work["REF"].astype("string[pyarrow]")

            if not str(df_work["ALT"].dtype).startswith("string"):
                df_work["ALT"] = df_work["ALT"].astype("string[pyarrow]")

            if not str(df_work["FILTER"].dtype).startswith("string"):
                df_work["FILTER"] = df_work["FILTER"].astype("string[pyarrow]")

            # Create masks using pyarrow engine internally
            chrom_mask = df_work["CHROM"] == chrom_formatted
            pos_mask = df_work["POS"] == pos
            ref_mask = df_work["REF"] == ref
            alt_mask = df_work["ALT"] == alt

            # Combine masks to find the specific record
            record_mask = chrom_mask & pos_mask & ref_mask & alt_mask

            # Check if record exists
            matching_records = df_work[record_mask]

            if len(matching_records) == 0:
                logger.info(f"Record not found: {chrom_formatted}:{pos} {ref}>{alt}")
                return False
            elif len(matching_records) > 1:
                logger.warning(f"Multiple records found for: {chrom_formatted}:{pos} {ref}>{alt}")

            # Check if FILTER is PASS for the matching record(s)
            pass_records = matching_records[matching_records["FILTER"] == "PASS"]
            result = len(pass_records) > 0

            logger.info(f"PASS filter result: {result}")
            return result

        except ImportError:
            logger.warning("PyArrow not available, using standard pandas operations")

            df_formatted = df.copy()
            df_formatted["CHROM"] = df_formatted["CHROM"].astype(str).map(format_chr)

            mask = (
                    (df_formatted["CHROM"] == chrom_formatted) &
                    (df_formatted["POS"] == pos) &
                    (df_formatted["REF"] == ref) &
                    (df_formatted["ALT"] == alt)
            )

            matching_records = df_formatted[mask]

            if len(matching_records) == 0:
                logger.info(f"Record not found: {chrom_formatted}:{pos} {ref}>{alt}")
                return False
            elif len(matching_records) > 1:
                logger.warning(f"Multiple records found for: {chrom_formatted}:{pos} {ref}>{alt}")

            # Check if FILTER is PASS
            result = (matching_records["FILTER"] == "PASS").any()
            logger.info(f"PASS filter result: {result}")
            return result

        except Exception as e:
            logger.warning(f"PyArrow optimization failed: {e}, falling back to standard operations")

            df_formatted = df.copy()
            df_formatted["CHROM"] = df_formatted["CHROM"].astype(str).map(format_chr)

            mask = (
                    (df_formatted["CHROM"] == chrom_formatted) &
                    (df_formatted["POS"] == pos) &
                    (df_formatted["REF"] == ref) &
                    (df_formatted["ALT"] == alt)
            )

            matching_records = df_formatted[mask]

            if len(matching_records) == 0:
                logger.info(f"Record not found: {chrom_formatted}:{pos} {ref}>{alt}")
                return False
            elif len(matching_records) > 1:
                logger.warning(f"Multiple records found for: {chrom_formatted}:{pos} {ref}>{alt}")

            # Check if FILTER is PASS
            result = (matching_records["FILTER"] == "PASS").any()
            logger.info(f"PASS filter result: {result}")
            return result
