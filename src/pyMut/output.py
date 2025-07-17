import logging
from pathlib import Path
from typing import List
import pandas as pd

from .input import required_columns_MAF
from .utils.format import reverse_format_chr

logger = logging.getLogger(__name__)


def _load_maf_column_order() -> List[str]:
    """
    Load the MAF column order from MAF_COL_ORDER.csv file.

    Returns
    -------
    List[str]
        List of column names in the order specified in MAF_COL_ORDER.csv
    """
    maf_columns_file = Path(__file__).parent / "data" / "MAF_COL_ORDER.csv"

    try:
        if maf_columns_file.exists():
            columns_df = pd.read_csv(maf_columns_file)
            # Extract column names from the 'nombre' column, ordered by 'id'
            ordered_columns = columns_df.sort_values('id')['nombre'].tolist()
            logger.debug(f"Loaded {len(ordered_columns)} column names from MAF_COL_ORDER.csv")
            return ordered_columns
        else:
            logger.warning(f"MAF_COL_ORDER.csv not found at {maf_columns_file}, using default order")
            return []
    except Exception as e:
        logger.warning(f"Error loading MAF_COL_ORDER.csv: {e}, using default order")
        return []




def to_maf(self, output_path: str | Path) -> None:
    """
    Export a PyMutation object back to MAF format.

    This function reverses the transformations done by read_maf() to recreate
    a MAF file from a PyMutation object. The output follows the column order
    specified in MAF_COL_ORDER.csv when available.

    Parameters
    ----------
    output_path : str | Path
        Path where the MAF file will be written.

    Raises
    ------
    ValueError
        If the PyMutation object doesn't contain the necessary data for MAF export.
    """
    output_path = Path(output_path)
    logger.info("Starting MAF export to: %s", output_path)

    # Get the data and samples from PyMutation object
    data = self.data.copy()
    samples = self.samples
    metadata = self.metadata

    # ─── 1) VALIDATE REQUIRED COLUMNS FOR EXPORT ──────────────────────
    vcf_like_cols = ["CHROM", "POS", "REF", "ALT"]
    missing_vcf_cols = [col for col in vcf_like_cols if col not in data.columns]
    if missing_vcf_cols:
        raise ValueError(f"Missing required VCF-style columns for MAF export: {missing_vcf_cols}")

    missing_samples = [sample for sample in samples if sample not in data.columns]
    if missing_samples:
        raise ValueError(f"Missing sample columns for MAF export: {missing_samples}")

    # ─── 2) CONVERT BACK TO MAF FORMAT ────────────────────────────────
    maf_rows = []

    for idx, row in data.iterrows():
        # Extract basic variant information
        chrom = row["CHROM"]
        pos = row["POS"]
        ref = row["REF"]
        alt = row["ALT"]

        # Process each sample
        for sample in samples:
            genotype = row[sample]

            # Skip if this sample doesn't have this variant (REF|REF)
            if genotype == f"{ref}|{ref}":
                continue

            # Parse genotype to get alleles
            if "|" in genotype:
                allele1, allele2 = genotype.split("|")
            elif "/" in genotype:
                allele1, allele2 = genotype.split("/")
            else:
                # Handle cases where genotype might not be properly formatted
                allele1, allele2 = ref, alt

            # Create MAF row
            maf_row = {
                "Chromosome": reverse_format_chr(str(chrom)),
                "Start_Position": int(pos),
                "Reference_Allele": str(ref),
                "Tumor_Seq_Allele1": str(allele1),
                "Tumor_Seq_Allele2": str(allele2),
                "Tumor_Sample_Barcode": str(sample),
            }

            # Add End_Position (for most variants, it's the same as Start_Position)
            if len(ref) == 1 and len(alt) == 1:
                # SNP
                maf_row["End_Position"] = pos
            elif len(ref) > len(alt):
                # Deletion
                maf_row["End_Position"] = pos + len(ref) - 1
            else:
                # Insertion or other
                maf_row["End_Position"] = pos

            # Add other columns that exist in the original data
            for col in data.columns:
                if col not in ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"] + samples:
                    maf_row[col] = row[col]

            maf_rows.append(maf_row)

    # ─── 3) CREATE DATAFRAME AND ENSURE COLUMN ORDER ─────────────────
    if not maf_rows:
        logger.warning("No variant data found to export")
        maf_df = pd.DataFrame(columns=required_columns_MAF)
    else:
        maf_df = pd.DataFrame(maf_rows)

    # Load the preferred column order from MAF_COL_ORDER.csv
    preferred_column_order = _load_maf_column_order()

    # Ensure required columns are present
    for col in required_columns_MAF:
        if col not in maf_df.columns:
            # Add missing required columns with default values
            maf_df[col] = "."

    # Determine final column order
    final_columns = []

    if preferred_column_order:
        # Use the order from MAF_COL_ORDER.csv, but only include columns that exist
        for col in preferred_column_order:
            if col in maf_df.columns:
                final_columns.append(col)

        # Add any remaining columns that weren't in the preferred order
        for col in maf_df.columns:
            if col not in final_columns:
                final_columns.append(col)

        logger.info(f"Using MAF_COL_ORDER.csv column order: {len(final_columns)} columns arranged")
    else:
        # Add required columns first
        for col in required_columns_MAF:
            if col in maf_df.columns:
                final_columns.append(col)

        # Add other common MAF columns if they exist
        common_maf_cols = ["Hugo_Symbol", "End_Position", "Variant_Classification", "Variant_Type"]
        for col in common_maf_cols:
            if col in maf_df.columns and col not in final_columns:
                final_columns.append(col)

        # Add any remaining columns
        for col in maf_df.columns:
            if col not in final_columns:
                final_columns.append(col)

        logger.info("Using default column order (MAF_COL_ORDER.csv not available)")

    maf_df = maf_df[final_columns]

    # ─── 4) WRITE TO FILE WITH COMMENTS ──────────────────────────────
    with open(output_path, 'w', encoding='utf-8') as f:
        # Write comments from metadata if they exist
        if metadata.notes:
            for line in metadata.notes.split('\n'):
                if line.strip():
                    if not line.startswith('#'):
                        f.write(f"#{line}\n")
                    else:
                        f.write(f"{line}\n")

        # Write the DataFrame as TSV without trailing newline
        maf_df.to_csv(f, sep='\t', index=False, lineterminator='\n')

        f.seek(f.tell() - 1)
        f.truncate()

    logger.info("MAF export completed successfully: %d rows written to %s", len(maf_df), output_path)


# Add the method to PyMutation class
def add_to_maf_method_to_pymutation():
    from .core import PyMutation
    PyMutation.to_maf = to_maf

add_to_maf_method_to_pymutation()
