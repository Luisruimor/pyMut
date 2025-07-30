import gzip
import logging
from pathlib import Path
from typing import Optional

# Logger configuration
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Change to DEBUG for more verbosity
if not logger.handlers:
    logging.basicConfig(
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        level=logging.INFO,
    )


def extract_vcf_subset(
        input_vcf_path: str | Path,
        output_vcf_path: str | Path,
        max_variants: Optional[int] = None,
        chromosome: Optional[str] = None,
        start_pos: Optional[int] = None,
        end_pos: Optional[int] = None
) -> bool:
    """
    Extract a subset of variants from a VCF file

    This function creates a smaller VCF file from a large VCF file without requiring
    external tools like tabix. It works with both compressed (.gz) and uncompressed files.

    Parameters
    ----------
    input_vcf_path : str | Path
        Path to the input VCF file (.vcf or .vcf.gz).
    output_vcf_path : str | Path
        Path where the subset VCF file will be saved.
    max_variants : int, optional
        Maximum number of variants to extract. If specified, will extract the first
        N variants from the specified region or entire file.
    chromosome : str, optional
        Chromosome to extract (e.g., "10", "chr10"). If None, extracts from all chromosomes.
    start_pos : int, optional
        Start position for extraction. Used with chromosome parameter.
    end_pos : int, optional
        End position for extraction. Used with chromosome parameter.

    Returns
    -------
    bool
        True if extraction was successful, False otherwise.

    Examples
    --------
    # Extract first 10000 variants from chromosome 10
    extract_vcf_subset(
        "large_file.vcf.gz",
        "subset.vcf",
        chromosome="10",
        max_variants=10000
    )

    # Extract first 50000 variants from entire file
    extract_vcf_subset(
        "large_file.vcf.gz",
        "subset.vcf",
        max_variants=50000
    )
    """
    input_path = Path(input_vcf_path)
    output_path = Path(output_vcf_path)

    if not input_path.exists():
        logger.error("Input VCF file does not exist: %s", input_path)
        return False

    try:
        logger.info("Extracting VCF subset")
        logger.info("Input: %s", input_path)
        logger.info("Output: %s", output_path)
        logger.info("Max variants: %s", max_variants)
        logger.info("Chromosome filter: %s", chromosome)
        logger.info("Position range: %s-%s", start_pos, end_pos)

        if str(input_path).endswith('.gz'):
            input_file = gzip.open(input_path, 'rt', encoding='utf-8')
        else:
            input_file = open(input_path, 'r', encoding='utf-8')

        output_file = open(output_path, 'w', encoding='utf-8')

        variants_written = 0
        header_written = False

        try:
            for line_num, line in enumerate(input_file, 1):
                line = line.strip()

                # Skip empty lines
                if not line:
                    continue

                # Handle header lines
                if line.startswith('#'):
                    output_file.write(line + '\n')
                    if line.startswith('#CHROM'):
                        header_written = True
                    continue

                # Process variant lines
                if not header_written:
                    logger.warning("No header found in VCF file")
                    return False

                # Check if we've reached the maximum number of variants
                if max_variants and variants_written >= max_variants:
                    break

                # Parse the variant line
                fields = line.split('\t')
                if len(fields) < 8:  # VCF requires at least 8 columns
                    continue

                chrom = fields[0]
                pos = int(fields[1])

                # Apply chromosome filter
                if chromosome:
                    # Handle both "10" and "chr10" formats
                    chrom_normalized = chrom.replace('chr', '')
                    chromosome_normalized = chromosome.replace('chr', '')
                    if chrom_normalized != chromosome_normalized:
                        continue

                # Apply position filter
                if start_pos is not None and pos < start_pos:
                    continue
                if end_pos is not None and pos > end_pos:
                    continue

                # Write the variant
                output_file.write(line + '\n')
                variants_written += 1

                if variants_written % 10000 == 0:
                    logger.info("Processed %d variants...", variants_written)

        finally:
            input_file.close()
            output_file.close()

        input_size_mb = input_path.stat().st_size / (1024 * 1024)
        output_size_mb = output_path.stat().st_size / (1024 * 1024)

        logger.info("VCF subset extraction completed successfully")
        logger.info("Variants extracted: %d", variants_written)
        logger.info("Input file size: %.2f MB", input_size_mb)
        logger.info("Output file size: %.2f MB", output_size_mb)

        if output_size_mb > 0:
            logger.info("Size reduction: %.1fx", input_size_mb / output_size_mb)

        return True

    except Exception as e:
        logger.error("Unexpected error during VCF subset extraction: %s", e)
        return False
