import logging
from pathlib import Path
from typing import List, Optional

import pandas as pd

from .input import required_columns_MAF

# Import PyArrow if available
try:
    import pyarrow as pa
    import pyarrow.compute as pc

    HAS_PYARROW = True
except ImportError:
    HAS_PYARROW = False

logger = logging.getLogger(__name__)

# Global variable to cache the column order
_MAF_COLUMN_ORDER_CACHE: Optional[List[str]] = None


def _load_maf_column_order() -> List[str]:
    """
    Load the MAF column order from MAF_COL_ORDER.csv file.
    Uses a cached version if available to avoid reading the file multiple times.

    Returns
    -------
    List[str]
        List of column names in the order specified in MAF_COL_ORDER.csv
    """
    global _MAF_COLUMN_ORDER_CACHE

    # Use cache if it already exists
    if _MAF_COLUMN_ORDER_CACHE is not None:
        return _MAF_COLUMN_ORDER_CACHE

    maf_columns_file = Path(__file__).parent / "data" / "MAF_COL_ORDER.csv"

    try:
        if maf_columns_file.exists():
            columns_df = pd.read_csv(maf_columns_file)
            # Extract column names from the 'nombre' column, ordered by 'id'
            ordered_columns = columns_df.sort_values('id')['nombre'].tolist()
            logger.debug(f"Loaded {len(ordered_columns)} column names from MAF_COL_ORDER.csv")

            # Save to cache
            _MAF_COLUMN_ORDER_CACHE = ordered_columns
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

    # Log total number of variants to process
    total_variants = len(data)
    logger.info(f"Starting to process {total_variants} variants from {len(samples)} samples")

    # ─── VALIDATE REQUIRED COLUMNS FOR EXPORT ──────────────────────
    vcf_like_cols = ["CHROM", "POS", "REF", "ALT", "ID"]
    missing_vcf_cols = [col for col in vcf_like_cols if col not in data.columns]
    if missing_vcf_cols:
        raise ValueError(f"Missing required VCF-style columns for MAF export: {missing_vcf_cols}")

    missing_samples = [sample for sample in samples if sample not in data.columns]
    if missing_samples:
        raise ValueError(f"Missing sample columns for MAF export: {missing_samples}")

    # ─── CONVERT TO MAF FORMAT ────────────────────────────────
    # Determine whether to use PyArrow for large datasets
    use_pyarrow = HAS_PYARROW and len(data) > 10000

    if use_pyarrow:
        logger.info(f"Using PyArrow for processing large dataset ({len(data)} rows)")

        # Convert to PyArrow table for faster processing
        table = pa.Table.from_pandas(data)

        # Create base columns using PyArrow and convert CHROM to Chromosome (extract chromosome number only)
        chrom_array = pc.cast(table['CHROM'], pa.string())
        # Remove 'chr' prefix if it exists
        chrom_array = pc.utf8_replace_substring(chrom_array, 'chr', '')
        table = table.append_column('Chromosome', chrom_array)

        # Other base columns
        table = table.append_column('Start_Position', table['POS'])
        table = table.append_column('Reference_Allele', pc.cast(table['REF'], pa.string()))
        table = table.append_column('NCBI_Build', pa.array([metadata.assembly] * len(table)))
        table = table.append_column('dbSNP_RS', pc.cast(table['ID'], pa.string()))

        # Convert back to pandas for specific operations
        base_data = table.to_pandas()
    else:
        # Create base columns for all variants using pandas
        base_data = data.copy()
        base_data['Chromosome'] = base_data['CHROM'].str.replace('chr', '', regex=False)
        base_data['Start_Position'] = base_data['POS']
        base_data['Reference_Allele'] = base_data['REF'].astype(str)
        base_data['NCBI_Build'] = metadata.assembly
        base_data['dbSNP_RS'] = base_data['ID'].astype(str)

    # Process each sample and create a DataFrame per sample
    sample_dfs = []
    processed_samples = 0
    total_samples = len(samples)
    log_frequency = max(1, total_samples // 50)  # Log every ~2% of samples, minimum 1

    for sample_idx, sample in enumerate(samples, 1):
        # Log progress for sample processing (only every log_frequency samples)
        if sample_idx % log_frequency == 0 or sample_idx == 1 or sample_idx == total_samples:
            logger.info(
                f"Processing sample {sample_idx}/{total_samples}: {sample} ({sample_idx / total_samples * 100:.1f}%)")

        # Create a copy of the data for this sample
        sample_data = base_data.copy()
        sample_data['Tumor_Sample_Barcode'] = sample

        # Filter only variants where the genotype is not REF|REF
        ref_pattern = sample_data['REF'] + '|' + sample_data['REF']
        sample_data = sample_data[sample_data[sample] != ref_pattern]

        # Increment processed samples counter
        processed_samples += 1

        # Log number of variants for this sample (only every log_frequency samples)
        if sample_idx % log_frequency == 0 or sample_idx == 1 or sample_idx == total_samples:
            sample_variants = len(sample_data)
            logger.info(f"Sample {sample}: {sample_variants} variants found")

        if len(sample_data) == 0:
            continue

        # Process genotypes in a vectorized way
        if '|' in sample_data[sample].iloc[0]:
            sample_data[['Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']] = sample_data[sample].str.split('|', expand=True)
        elif '/' in sample_data[sample].iloc[0]:
            sample_data[['Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']] = sample_data[sample].str.split('/', expand=True)
        else:
            sample_data['Tumor_Seq_Allele1'] = sample_data['REF']
            sample_data['Tumor_Seq_Allele2'] = sample_data['ALT']

        # Calculate End_Position in a vectorized way for different variant types
        snp_mask = (sample_data['Reference_Allele'].str.len() == 1) & (sample_data['Tumor_Seq_Allele2'].str.len() == 1)
        sample_data.loc[snp_mask, 'End_Position'] = sample_data.loc[snp_mask, 'Start_Position'].astype(int)

        del_mask = sample_data['Reference_Allele'].str.len() > sample_data['Tumor_Seq_Allele2'].str.len()
        sample_data.loc[del_mask, 'End_Position'] = (sample_data.loc[del_mask, 'Start_Position'] + sample_data.loc[
            del_mask, 'Reference_Allele'].str.len() - 1).astype(int)

        other_mask = ~(snp_mask | del_mask)
        sample_data.loc[other_mask, 'End_Position'] = sample_data.loc[other_mask, 'Start_Position'].astype(int)

        # Select relevant columns
        maf_cols = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele1',
                    'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 'NCBI_Build', 'dbSNP_RS']

        # Add other columns that are neither VCF nor samples
        other_cols = [col for col in data.columns if col not in vcf_like_cols + samples]
        all_cols = maf_cols + other_cols

        # Select only columns that exist in the data
        existing_cols = [col for col in all_cols if col in sample_data.columns]

        sample_dfs.append(sample_data[existing_cols])

    # Combine all sample DataFrames
    if not sample_dfs:
        logger.warning("No variant data found to export")
        maf_df = pd.DataFrame(columns=required_columns_MAF)
    else:
        maf_df = pd.concat(sample_dfs, ignore_index=True)

    # Log summary after processing all samples
    total_variants_found = len(maf_df)
    logger.info(f"Sample processing completed: {processed_samples}/{total_samples} samples processed")
    logger.info(f"Total variants found: {total_variants_found} variants")

    # ─── ENSURE COLUMN ORDER ─────────────────────

    # Load the preferred column order from MAF_COL_ORDER.csv
    preferred_column_order = _load_maf_column_order()

    # Ensure required columns are present
    for col in required_columns_MAF:
        if col not in maf_df.columns:
            # Add missing required columns with default values
            maf_df[col] = "."

    # Optimized column ordering
    if preferred_column_order:
        # Filter only columns that exist in the DataFrame
        existing_preferred_columns = [col for col in preferred_column_order if col in maf_df.columns]

        # Get columns that are not in the preferred order
        remaining_columns = list(set(maf_df.columns) - set(existing_preferred_columns))

        # Combine the lists
        final_columns = existing_preferred_columns + remaining_columns

        logger.info(f"Using MAF_COL_ORDER.csv column order: {len(final_columns)} columns arranged")

        # Reorder the DataFrame
        maf_df = maf_df[final_columns]

    # Ensure End_Position is an integer
    if 'End_Position' in maf_df.columns:
        maf_df['End_Position'] = maf_df['End_Position'].astype(int)

    # ─── WRITE TO FILE WITH COMMENTS ──────────────────────────────
    chunk_size = 10000

    with open(output_path, 'w', encoding='utf-8') as f:
        # Write only INFO comments from metadata if they exist
        if metadata.notes:
            for line in metadata.notes.split('\n'):
                if line.strip():
                    # Only include lines that start with INFO or ##INFO
                    if line.startswith('##INFO') or line.startswith('INFO'):
                        if not line.startswith('#'):
                            f.write(f"#{line}\n")
                        else:
                            f.write(f"{line}\n")

        # Write header
        f.write('\t'.join(maf_df.columns) + '\n')

        # For small datasets, write directly
        if len(maf_df) <= chunk_size:
            logger.info(f"Writing {len(maf_df)} variants to file")
            maf_df.to_csv(f, sep='\t', index=False, header=False, lineterminator='\n')
            logger.info(f"Progress: {len(maf_df)}/{len(maf_df)} variants written (100.0%)")
        else:
            # For large datasets, write in batches
            total_rows = len(maf_df)
            logger.info(f"Writing large dataset ({total_rows} variants) in chunks of {chunk_size}")

            for i in range(0, total_rows, chunk_size):
                chunk = maf_df.iloc[i:i + chunk_size]
                chunk.to_csv(f, sep='\t', index=False, header=False, lineterminator='\n')

                # Calculate the actual number of rows written (may be less than i+chunk_size in the last batch)
                rows_written = min(i + chunk_size, total_rows)

                # Log progress for each batch
                logger.info(
                    f"Progress: {rows_written}/{total_rows} variants written ({rows_written / total_rows * 100:.1f}%)")

        # Remove the last newline character if it exists
        if f.tell() > 0:
            f.seek(f.tell() - 1)
            f.truncate()

    # Log final summary with detailed statistics
    logger.info(f"MAF export completed successfully: {len(maf_df)} variants processed and written to {output_path}")
    logger.info(
        f"Conversion summary: {len(samples)} samples, {total_variants} input variants, {len(maf_df)} output variants")


# Add the method to PyMutation class
def to_vcf(self, output_path: str | Path) -> None:
    """
    Export a PyMutation object to VCF format.

    This function creates a VCF file from a PyMutation object, including
    a proper VCF header with metadata information.

    Parameters
    ----------
    output_path : str | Path
        Path where the VCF file will be written.

    Raises
    ------
    ValueError
        If the PyMutation object doesn't contain the necessary data for VCF export.
    """
    output_path = Path(output_path)
    logger.info("Starting VCF export to: %s", output_path)

    # Get the data and samples from PyMutation object
    data = self.data.copy()
    samples = self.samples
    metadata = self.metadata

    # Log total number of variants to process
    total_variants = len(data)
    logger.info(f"Starting to process {total_variants} variants from {len(samples)} samples")

    # ─── VALIDATE REQUIRED COLUMNS FOR EXPORT ──────────────────────
    vcf_like_cols = ["CHROM", "POS", "REF", "ALT", "ID"]
    missing_vcf_cols = [col for col in vcf_like_cols if col not in data.columns]
    if missing_vcf_cols:
        raise ValueError(f"Missing required VCF-style columns for VCF export: {missing_vcf_cols}")

    missing_samples = [sample for sample in samples if sample not in data.columns]
    if missing_samples:
        raise ValueError(f"Missing sample columns for VCF export: {missing_samples}")

    # ─── PREPARE DATA FOR VCF FORMAT ────────────────────────────────
    # Get unique chromosomes for contig lines
    unique_chroms = data['CHROM'].unique()

    # Get all columns that are not standard VCF columns and not sample columns
    standard_vcf_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    excluded_cols = standard_vcf_cols + samples
    info_cols = [col for col in data.columns if col not in excluded_cols]
    info_cols_str = "|".join(info_cols)

    # ─── WRITE TO FILE WITH HEADER ──────────────────────────────
    from datetime import datetime

    with open(output_path, 'w', encoding='utf-8') as f:
        # Generate VCF header
        f.write("##fileformat=VCFv4.3\n")
        current_date = datetime.now().strftime("%Y%m%d")
        f.write(f"##fileDate={current_date}\n")
        f.write("##source=https://github.com/Luisruimor/pyMut\n")
        f.write(f"##reference={metadata.assembly}\n")
        f.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")

        # Contig lines for each chromosome
        for chrom in unique_chroms:
            formatted_chrom = chrom
            if chrom.startswith('chr'):
                formatted_chrom = chrom[3:]
            f.write(f"##contig=<ID={formatted_chrom}>\n")

        # Format and INFO fields
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased Genotype\">\n")
        f.write(
            f"##INFO=<ID=PMUT,Number=.,Type=String,Description=\"Consequence annotations columns from PyMut. Format: {info_cols_str}\">\n")

        # Write VCF column headers
        vcf_columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples
        f.write('\t'.join(vcf_columns) + '\n')

        # Prepare data for writing
        vcf_data = data.copy()

        # Format chromosome values to remove 'chr' prefix
        vcf_data['CHROM'] = vcf_data['CHROM'].apply(lambda x: x[3:] if x.startswith('chr') else x)

        # Set default values for missing columns
        if "QUAL" not in vcf_data.columns:
            vcf_data["QUAL"] = "."
        if "FILTER" not in vcf_data.columns:
            vcf_data["FILTER"] = "PASS"
        if "FORMAT" not in vcf_data.columns:
            vcf_data["FORMAT"] = "GT"

        # Populate INFO column with PMUT values
        if info_cols:
            vcf_data["INFO"] = vcf_data.apply(
                lambda row: f"PMUT={('|'.join([str(row[col]) if pd.notna(row[col]) else '.' for col in info_cols]))}",
                axis=1)
        elif "INFO" not in vcf_data.columns:
            vcf_data["INFO"] = "."

        # Process genotype data to replace bases with indices
        logger.info("Processing genotype data to replace bases with indices")
        for i, row in vcf_data.iterrows():
            ref = row['REF']
            alt_alleles = row['ALT'].split(',') if row['ALT'] else []

            # Build allele map {REF:0, ALT1:1, ALT2:2,...}
            allele_map = {ref: '0'}
            for idx, alt in enumerate(alt_alleles, 1):
                allele_map[alt] = str(idx)

            # Process each sample's genotype
            for sample in samples:
                gt = row[sample]

                # Skip if genotype is missing or no-call
                if pd.isna(gt) or gt in [".", "./.", ".|."]:
                    continue

                # Keep only the part before ':' if it exists
                gt_core = gt.split(":", 1)[0]

                # Determine the separator used ('|' preferentially)
                sep = "|" if "|" in gt_core else "/"

                # Skip if not using the '|' separator
                if sep != "|":
                    continue

                # Split genotype by separator
                alleles = gt_core.split(sep)
                new_alleles = []

                # Replace each base with its index according to the map
                for allele in alleles:
                    if allele == ".":
                        new_alleles.append(".")
                    elif allele in allele_map:
                        new_alleles.append(allele_map[allele])
                    else:
                        # If a base appears that is not in REF/ALT, add it to ALT and update indices
                        alt_alleles.append(allele)
                        new_idx = len(alt_alleles)
                        allele_map[allele] = str(new_idx)
                        new_alleles.append(str(new_idx))
                        # Update the ALT column
                        vcf_data.at[i, 'ALT'] = ','.join(alt_alleles)

                # Join the new alleles with the separator
                vcf_data.at[i, sample] = sep.join(new_alleles)

        # Reorder columns to match VCF format
        vcf_columns_data = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples
        existing_columns = [col for col in vcf_columns_data if col in vcf_data.columns]
        vcf_data = vcf_data[existing_columns]

        # Write data
        chunk_size = 10000

        # For small datasets, write directly
        if len(vcf_data) <= chunk_size:
            logger.info(f"Writing {len(vcf_data)} variants to file")
            vcf_data.to_csv(f, sep='\t', index=False, header=False, lineterminator='\n')
            logger.info(f"Progress: {len(vcf_data)}/{len(vcf_data)} variants written (100.0%)")
        else:
            # For large datasets, write in batches
            total_rows = len(vcf_data)
            logger.info(f"Writing large dataset ({total_rows} variants) in chunks of {chunk_size}")

            for i in range(0, total_rows, chunk_size):
                chunk = vcf_data.iloc[i:i + chunk_size]
                chunk.to_csv(f, sep='\t', index=False, header=False, lineterminator='\n')

                # Calculate the actual number of rows written (may be less than i+chunk_size in the last batch)
                rows_written = min(i + chunk_size, total_rows)

                # Log progress for each batch
                logger.info(
                    f"Progress: {rows_written}/{total_rows} variants written ({rows_written / total_rows * 100:.1f}%)")

        # Remove the last newline character if it exists
        if f.tell() > 0:
            f.seek(f.tell() - 1)
            f.truncate()

    # Log final summary with detailed statistics
    logger.info(f"VCF export completed successfully: {len(vcf_data)} variants processed and written to {output_path}")
    logger.info(
        f"Conversion summary: {len(samples)} samples, {total_variants} input variants, {len(vcf_data)} output variants")


def add_to_maf_method_to_pymutation():
    from .core import PyMutation
    PyMutation.to_maf = to_maf


def add_to_vcf_method_to_pymutation():
    from .core import PyMutation
    PyMutation.to_vcf = to_vcf


add_to_maf_method_to_pymutation()
add_to_vcf_method_to_pymutation()
