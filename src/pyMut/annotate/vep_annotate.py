import csv
import gzip
import logging
import re
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Union, Optional, Tuple

from ..utils.format import format_chr
from ..utils.merge_vep_annotation import merge_maf_with_vep_annotations

logger = logging.getLogger(__name__)


def _extract_assembly_and_version_from_cache(cache_dir: Union[str, Path]) -> tuple[str, str]:
    """Extract assembly and version information from VEP cache directory name."""
    cache_path = Path(cache_dir)
    cache_name = cache_path.name

    pattern = r'homo_sapiens_vep_(\d+)_([A-Za-z0-9\.]+)'
    match = re.search(pattern, cache_name)

    if not match:
        raise ValueError(
            f"Cache directory name '{cache_name}' doesn't match expected format 'homo_sapiens_vep_{{version}}_{{assembly}}'")

    version = match.group(1)
    assembly = match.group(2)

    return assembly, version


def _get_case_insensitive_column(columns: list, target_column: str) -> str:
    """Find column name case-insensitively."""
    column_map = {col.lower(): col for col in columns}
    target_lower = target_column.lower()
    if target_lower in column_map:
        return column_map[target_lower]
    else:
        raise KeyError(f"Column '{target_column}' not found in MAF file. Available columns: {columns}")


def _maf_to_region(maf_path: Union[str, Path],
                   out_path: Optional[Union[str, Path]] = None) -> Tuple[bool, str]:
    """
    Convert a MAF file to region format.

    Args:
        maf_path: Path to the MAF file (.maf or .maf.gz)
        out_path: Path to the output .region file (optional)

    Returns:
        Tuple[bool, str]: (success_status, output_path)
    """
    maf_file = Path(maf_path)

    if not maf_file.exists():
        logger.error(f"MAF file not found: {maf_file}")
        return False, ""

    if out_path is None:
        if maf_file.suffix == '.gz' and maf_file.stem.endswith('.maf'):
            base_name = maf_file.stem[:-4]
            output_file = maf_file.parent / f"{base_name}.region"
        elif maf_file.suffix == '.maf':
            output_file = maf_file.with_suffix('.region')
        else:
            output_file = maf_file.with_suffix('.region')
    else:
        output_file = Path(out_path)

    logger.info(f"Converting MAF to region format: {maf_file} -> {output_file}")

    try:
        if maf_file.suffix == '.gz':
            file_opener = lambda: gzip.open(maf_file, 'rt', encoding='utf-8', newline='')
        else:
            file_opener = lambda: open(maf_file, 'r', encoding='utf-8', newline='')

        with file_opener() as maf, open(output_file, "w") as out:
            reader = csv.DictReader(maf, delimiter="\t")

            columns = reader.fieldnames
            if not columns:
                raise ValueError("No columns found in MAF file")

            try:
                chrom_col = _get_case_insensitive_column(columns, "Chromosome")
                start_col = _get_case_insensitive_column(columns, "Start_Position")
                end_col = _get_case_insensitive_column(columns, "End_position")
                alt_col = _get_case_insensitive_column(columns, "Tumor_Seq_Allele2")

                try:
                    strand_col = _get_case_insensitive_column(columns, "Strand")
                    has_strand = True
                except KeyError:
                    has_strand = False
                    logger.warning("Strand column not found in MAF file, using default value '+'")

            except KeyError as e:
                logger.error(f"Required column not found: {e}")
                raise

            for row in reader:
                chrom_raw = row[chrom_col]
                chrom = format_chr(chrom_raw)
                start = row[start_col]
                endpos = row[end_col]
                alt = row[alt_col]
                strand = row[strand_col] if has_strand else "+"
                strand_num = "1" if strand == "+" else "-1"

                out.write(f"{chrom}:{start}-{endpos}:{strand_num}/{alt}\n")

        logger.info(f"Successfully converted MAF to region format: {output_file}")
        return True, str(output_file)

    except Exception as e:
        logger.error(f"Error converting MAF to region format: {e}")
        return False, str(output_file) if 'output_file' in locals() else ""


def wrap_maf_vep_annotate_protein(maf_file: Union[str, Path],
                                  cache_dir: Union[str, Path],
                                  fasta: Union[str, Path],
                                  output_file: Optional[Union[str, Path]] = None,
                                  synonyms_file: Optional[Union[str, Path]] = None,
                                  assembly: Optional[str] = None,
                                  version: Optional[str] = None,
                                  compress: bool = True,
                                  no_stats: bool = True) -> Tuple[bool, str]:
    """
    Wrapper method for VEP annotation that accepts MAF files and merges annotations back to MAF.

    This method converts a MAF file to region format internally, runs VEP annotation, and then
    merges the VEP annotations back with the original MAF file. VEP annotation uses the following
    fixed parameters:
    - --offline --cache
    - --protein --uniprot --domains --symbol
    - --synonyms (automatically constructed from cache directory or provided explicitly)
    - --no_stats (only when no_stats=False)

    After successful VEP annotation, the method automatically merges the VEP results with the
    original MAF file, creating an annotated MAF file with VEP_ prefixed columns. The original
    VEP annotation files are preserved and not deleted.

    Assembly and cache version can be provided explicitly or automatically extracted from the cache directory name.
    The chr_synonyms file path can be provided explicitly or automatically constructed as: cache_dir/homo_sapiens/{version}_{assembly}/chr_synonyms.txt

    Args:
        maf_file: Path to the MAF file to annotate (.maf or .maf.gz)
        cache_dir: Path to the VEP cache directory
        fasta: Path to the reference FASTA file
        output_file: Path to the output file (optional). If None, creates a directory
                    in the same location as maf_file with format 'vep_annotation_HHMMDDMMYYYY'
        synonyms_file: Path to the chromosome synonyms file (optional). If None, automatically
                      constructed from cache directory structure
        assembly: Genome assembly name (optional). If None, automatically extracted from cache directory name
        version: VEP cache version (optional). If None, automatically extracted from cache directory name
        compress: Whether to compress the merged output file with gzip (default: True)
        no_stats: Whether to disable VEP statistics generation (default: True). When True, --no_stats flag is omitted

    Returns:
        Tuple[bool, str]: (success_status, output_info) where success_status is True
                         if annotation was successful, and output_info contains information
                         about both the VEP file and the merged MAF file paths

    Raises:
        ValueError: If cache directory name format is invalid and assembly/version not provided
        FileNotFoundError: If required files don't exist
    """
    maf_path = Path(maf_file)
    cache_path = Path(cache_dir)
    fasta_path = Path(fasta)
    if not maf_path.exists():
        raise FileNotFoundError(f"MAF file not found: {maf_path}")
    if not cache_path.exists():
        raise FileNotFoundError(f"Cache directory not found: {cache_path}")
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    logger.info(f"Converting MAF file to region format: {maf_path}")
    region_success, region_path = _maf_to_region(maf_path)

    if not region_success:
        logger.error("Failed to convert MAF to region format")
        return False, ""

    region_file = Path(region_path)
    logger.info(f"Successfully converted MAF to region format: {region_file}")

    if output_file is None:
        timestamp = datetime.now().strftime("%H%M%d%m")
        output_dir_name = f"vep_annotation_{timestamp}"
        output_dir = maf_path.parent / output_dir_name
        output_dir.mkdir(exist_ok=True)

        output_filename = f"{maf_path.stem}_vep_protein.txt"
        output_path = output_dir / output_filename
    else:
        output_path = Path(output_file)

    if assembly is None or version is None:
        try:
            extracted_assembly, extracted_version = _extract_assembly_and_version_from_cache(cache_path)
            if assembly is None:
                assembly = extracted_assembly
            if version is None:
                version = extracted_version
            logger.info(f"Extracted from cache: assembly={assembly}, version={version}")
        except ValueError as e:
            logger.error(f"Failed to extract assembly/version from cache: {e}")
            raise
    else:
        logger.info(f"Using provided: assembly={assembly}, version={version}")

    if synonyms_file is None:
        chr_synonyms_path = cache_path / "homo_sapiens" / f"{version}_{assembly}" / "chr_synonyms.txt"
        logger.info(f"Auto-constructed chr synonyms path: {chr_synonyms_path}")
    else:
        chr_synonyms_path = Path(synonyms_file)
        logger.info(f"Using provided chr synonyms path: {chr_synonyms_path}")

    vep_cmd = [
        "vep",
        "--input_file", str(region_file),
        "--format", "region",
        "--offline", "--cache", "--cache_version", version,
        "--dir_cache", str(cache_path),
        "--assembly", assembly,
        "--synonyms", str(chr_synonyms_path),
        "--fasta", str(fasta_path),
        "--protein", "--uniprot", "--domains", "--symbol",
        "--pick",
        "--keep_csq",
        "--force_overwrite",
        "--output_file", str(output_path)
    ]

    # Add --no_stats when no_stats is True
    if no_stats:
        vep_cmd.insert(-2, "--no_stats")

    try:
        logger.info(f"Running VEP annotation: {' '.join(vep_cmd)}")
        result = subprocess.run(vep_cmd, check=True, capture_output=True, text=True)
        logger.info("VEP annotation completed successfully")

        if result.stderr:
            logger.warning(f"VEP warnings/messages: {result.stderr}")

        logger.info("Merging VEP annotations with original MAF file...")
        try:
            merged_df, merged_output_path = merge_maf_with_vep_annotations(
                maf_file=maf_path,
                vep_file=output_path,
                output_file=None,
                compress=compress
            )
            logger.info(f"Successfully merged VEP annotations. Merged file: {merged_output_path}")

            return True, f"VEP folder: {output_path}, Merged file: {merged_output_path}"

        except Exception as merge_error:
            logger.error(f"Failed to merge VEP annotations: {merge_error}")
            # Success for VEP annotation, merge failure
            return True, f"VEP folder: {output_path}, Merge failed: {merge_error}"

    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.error(f"VEP annotation failed: {e}")
        if hasattr(e, 'stderr') and e.stderr:
            logger.error(f"VEP error output: {e.stderr}")
        return False, str(output_path)
    except Exception as e:
        logger.error(f"Unexpected error during VEP annotation: {e}")
        return False, str(output_path)
    finally:
        # Clean up temporary region file
        try:
            if region_file.exists():
                region_file.unlink()
                logger.info(f"Cleaned up temporary region file: {region_file}")
        except Exception as e:
            logger.warning(f"Failed to clean up temporary region file {region_file}: {e}")


def wrap_vcf_vep_annotate_unified(vcf_file: Union[str, Path],
                                  cache_dir: Union[str, Path],
                                  fasta: Union[str, Path],
                                  output_file: Optional[Union[str, Path]] = None,
                                  synonyms_file: Optional[Union[str, Path]] = None,
                                  assembly: Optional[str] = None,
                                  version: Optional[str] = None,
                                  no_stats: bool = True,
                                  # New parameters to control annotations
                                  annotate_protein: bool = False,
                                  annotate_gene: bool = False,
                                  annotate_variant_class: bool = False,
                                  distance: Optional[int] = None) -> Tuple[bool, str]:
    """
    Unified method for VEP annotation that allows combining different types of annotation.
    
    Args:
        vcf_file: Path to the VCF file to annotate
        cache_dir: Path to the VEP cache directory
        fasta: Path to the reference FASTA file
        output_file: Path to the output file (optional)
        synonyms_file: Path to the chromosome synonyms file (optional)
        assembly: Genome assembly name (optional)
        version: VEP cache version (optional)
        no_stats: Whether to disable VEP statistics generation
        annotate_protein: Whether to include protein annotation (--protein --uniprot --domains --symbol)
        annotate_gene: Whether to include gene annotation (--symbol)
        annotate_variant_class: Whether to include variant classification (--variant_class)
        distance: Distance for nearest gene search (only for annotate_gene)
    
    Returns:
        Tuple[bool, str]: (success_status, output_info)
    """
    # Validation: at least one annotation must be enabled
    if not any([annotate_protein, annotate_gene, annotate_variant_class]):
        raise ValueError("At least one annotation option must be enabled")

    vcf_path = Path(vcf_file)
    cache_path = Path(cache_dir)
    fasta_path = Path(fasta)

    # Validate input files
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    if not cache_path.exists():
        raise FileNotFoundError(f"Cache directory not found: {cache_path}")
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    logger.info(f"Starting unified VEP annotation for VCF file: {vcf_path}")

    # Handle output file creation
    if output_file is None:
        timestamp = datetime.now().strftime("%H%M%d%m")
        output_dir_name = f"vep_annotation_{timestamp}"
        output_dir = vcf_path.parent / output_dir_name
        output_dir.mkdir(exist_ok=True)

        # Create descriptive name based on selected annotations
        annotations = []
        if annotate_protein:
            annotations.append("protein")
        if annotate_gene:
            annotations.append("gene")
        if annotate_variant_class:
            annotations.append("variant_class")

        output_filename = f"{vcf_path.stem}_vep_{'_'.join(annotations)}.vcf"
        output_path = output_dir / output_filename
    else:
        output_path = Path(output_file)

    # Extract assembly and version from cache if not provided
    if assembly is None or version is None:
        try:
            extracted_assembly, extracted_version = _extract_assembly_and_version_from_cache(cache_path)
            if assembly is None:
                assembly = extracted_assembly
            if version is None:
                version = extracted_version
            logger.info(f"Extracted from cache: assembly={assembly}, version={version}")
        except ValueError as e:
            logger.error(f"Failed to extract assembly/version from cache: {e}")
            raise
    else:
        logger.info(f"Using provided: assembly={assembly}, version={version}")

    # Handle chromosome synonyms file
    if synonyms_file is None:
        chr_synonyms_path = cache_path / "homo_sapiens" / f"{version}_{assembly}" / "chr_synonyms.txt"
        logger.info(f"Auto-constructed chr synonyms path: {chr_synonyms_path}")
    else:
        chr_synonyms_path = Path(synonyms_file)
        logger.info(f"Using provided chr synonyms path: {chr_synonyms_path}")

    # Build base VEP command
    vep_cmd = [
        "vep",
        "--input_file", str(vcf_path),
        "--vcf",
        "--offline",
        "--cache",
        "--cache_version", version,
        "--dir_cache", str(cache_path),
        "--assembly", assembly,
        "--synonyms", str(chr_synonyms_path),
        "--fasta", str(fasta_path),
        "--pick",
        "--force_overwrite",
        "--output_file", str(output_path)
    ]

    # Add specific parameters based on selected options
    if annotate_protein:
        vep_cmd.extend(["--protein", "--uniprot", "--domains", "--symbol"])

    if annotate_gene:
        if not annotate_protein:  # Avoid duplicating --symbol
            vep_cmd.append("--symbol")

        # Add distance parameters if specified
        if distance is not None:
            vep_cmd.extend(["--nearest", "symbol", "--distance", str(distance)])

    if annotate_variant_class:
        # Add --variant_class WITHOUT --fields to get all variant class fields
        vep_cmd.append("--variant_class")

    # Add --no_stats when no_stats is True
    if no_stats:
        vep_cmd.append("--no_stats")

    try:
        logger.info(f"Running unified VEP annotation: {' '.join(vep_cmd)}")
        result = subprocess.run(vep_cmd, check=True, capture_output=True, text=True)
        logger.info("Unified VEP annotation completed successfully")

        if result.stderr:
            logger.warning(f"VEP warnings/messages: {result.stderr}")

        return True, f"VEP output file: {output_path}"

    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.error(f"Unified VEP annotation failed: {e}")
        if hasattr(e, 'stderr') and e.stderr:
            logger.error(f"VEP error output: {e.stderr}")
        return False, str(output_path)
    except Exception as e:
        logger.error(f"Unexpected error during unified VEP annotation: {e}")
        return False, str(output_path)
