import gzip
import logging
import os
import subprocess
import shutil
from datetime import datetime
from pathlib import Path
from typing import Union, Literal, Tuple, List, Optional

from ..utils.compatibility import (
    check_numpy_compatibility, 
    check_pandas_compatibility,
    safe_import_with_fallback,
    print_compatibility_report
)

pd, pandas_available = safe_import_with_fallback('pandas')
if not pandas_available:
    numpy_compatible, numpy_error = check_numpy_compatibility()
    pandas_compatible, pandas_error = check_pandas_compatibility()

    if not numpy_compatible:
        logging.warning("NumPy compatibility issue detected!")
        logging.warning(numpy_error)
        logging.info("For a full compatibility report, run:")
        logging.info("python -c 'from pyMut.utils.compatibility import print_compatibility_report; print_compatibility_report()'")
    elif not pandas_compatible:
        logging.warning("Pandas compatibility issue detected!")
        logging.warning(pandas_error)

    raise ImportError("pandas could not be imported due to compatibility issues. See error messages above for solutions.")

try:
    from pyensembl import EnsemblRelease
except ImportError:
    EnsemblRelease = None

from ..utils.fields import col
from ..utils.format import reverse_format_chr


def _generate_transcript_regions(genome_version: str) -> pd.DataFrame:
    """
    Generate transcript regions DataFrame using pyensembl.

    Parameters
    ----------
    genome_version : str
        Genome version ("GRCh37" or "GRCh38")

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: CHROMOSOME, START, END, ELEMENT, SYMBOL, STRAND
    """
    if EnsemblRelease is None:
        raise ImportError("pyensembl is required for transcript regions generation. Install with: pip install pyensembl")

    if genome_version == "GRCh37":
        release_number = 75
    elif genome_version == "GRCh38":
        release_number = 114
    else:
        raise ValueError(f"Unsupported genome version: {genome_version}. Use 'GRCh37' or 'GRCh38'")

    logging.info(f"Initializing Ensembl release {release_number} for {genome_version}...")
    ens = EnsemblRelease(release_number)

    try:
        ens.download()
        ens.index()
    except Exception as e:
        logging.warning(f"Error during download/indexing: {e}")
        logging.info("Attempting to use existing cached data...")

    logging.info("Generating transcript regions...")

    transcript_data = []

    try:
        for transcript in ens.transcripts():
            try:
                transcript_data.append({
                    "CHROMOSOME": transcript.contig,
                    "START": transcript.start,
                    "END": transcript.end,
                    "ELEMENT": transcript.transcript_id,
                    "SYMBOL": transcript.gene_name,
                    "STRAND": "+" if transcript.strand == 1 else "-"
                })
            except Exception as e:
                continue

    except Exception as e:
        raise RuntimeError(f"Error accessing transcript data: {e}")

    if not transcript_data:
        raise RuntimeError("No transcript data could be retrieved")

    regions_df = pd.DataFrame(transcript_data)
    logging.info(f"Generated {len(regions_df)} transcript regions")

    return regions_df


def _validate_chromosome_format(mutations_file: Path, regions_file: Path) -> bool:
    """
    Validate that CHROMOSOME uses the same format in both mutations and regions files.

    Parameters
    ----------
    mutations_file : Path
        Path to mutations.tsv file
    regions_file : Path
        Path to regions.gz file

    Returns
    -------
    bool
        True if formats are consistent, False otherwise
    """
    try:
        mutations_df = pd.read_csv(mutations_file, sep='\t', nrows=100)
        mutations_chrs = set(mutations_df['CHROMOSOME'].astype(str).unique())

        regions_df = pd.read_csv(regions_file, sep='\t', compression='gzip', nrows=100)
        regions_chrs = set(regions_df['CHROMOSOME'].astype(str).unique())

        mutations_has_chr = any(chr_val.startswith('chr') for chr_val in mutations_chrs)
        regions_has_chr = any(chr_val.startswith('chr') for chr_val in regions_chrs)

        if mutations_has_chr != regions_has_chr:
            logging.warning("Chromosome format mismatch detected!")
            logging.warning(f"Mutations file format: {'chr prefix' if mutations_has_chr else 'no chr prefix'}")
            logging.warning(f"Regions file format: {'chr prefix' if regions_has_chr else 'no chr prefix'}")
            return False

        logging.info("✓ Chromosome formats are consistent between files")
        return True

    except Exception as e:
        logging.error(f"Error validating chromosome formats: {e}")
        return False


def _calculate_optimal_permutations(snv_count: int) -> int:
    """
    Calculate optimal number of permutations based on SNV count in cohort.

    Parameters
    ----------
    snv_count : int
        Number of SNVs in the cohort

    Returns
    -------
    int
        Recommended number of permutations
    """
    if snv_count < 10000:
        return 1000
    elif snv_count <= 50000:
        return 10000
    else:
        return 20000


def _run_oncodriveclustl(mutations_file: Path, regions_file: Path, genome_build: str, 
                        output_dir: Path, n_permutations: int, threads: int = 4) -> bool:
    """
    Run OncodriveCLUSTL analysis with real-time progress monitoring.

    Parameters
    ----------
    mutations_file : Path
        Path to mutations.tsv file
    regions_file : Path
        Path to regions.gz file
    genome_build : str
        Genome build (hg19 or hg38)
    output_dir : Path
        Output directory for results
    n_permutations : int
        Number of permutations
    threads : int, default 4
        Number of threads to use

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    try:
        if not shutil.which('oncodriveclustl'):
            raise RuntimeError("OncodriveCLUSTL not found in PATH. Please install OncodriveCLUSTL.")

        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            'oncodriveclustl',
            '-i', str(mutations_file),
            '-r', str(regions_file),
            '-g', genome_build,
            '-o', str(output_dir),
            '--concatenate',
            '-n', str(n_permutations),
            '--element-mutations', '3',
            '-c', str(threads)
        ]

        logging.info(f"Running OncodriveCLUSTL with command:")
        logging.info(' '.join(cmd))

        logger = logging.getLogger("oncodrive")
        logger.setLevel(logging.INFO)
        logger.addHandler(logging.StreamHandler())

        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in proc.stdout:
            logger.info(line.rstrip())
        proc.wait()

        if proc.returncode == 0:
            logging.info("✓ OncodriveCLUSTL completed successfully")
            return True
        else:
            logging.error(f"OncodriveCLUSTL failed with return code {proc.returncode}")
            return False

    except Exception as e:
        logging.error(f"Error running OncodriveCLUSTL: {e}")
        return False


def process_oncodriveclustl_results(results_dir: Path, threshold: float = 0.10) -> pd.DataFrame:
    """
    Process OncodriveCLUSTL results and filter by Q_EMPIRICAL threshold.

    Parameters
    ----------
    results_dir : Path
        OncodriveCLUSTL results directory
    threshold : float, default 0.10
        Q_EMPIRICAL threshold for filtering

    Returns
    -------
    pd.DataFrame
        Filtered results with Q_EMPIRICAL ≤ threshold, sorted by q-value
    """
    try:
        results_file = results_dir / "elements_results.txt"

        if not results_file.exists():
            raise FileNotFoundError(f"Results file not found: {results_file}")

        results_df = pd.read_csv(results_file, sep='\t')

        if 'Q_EMPIRICAL' not in results_df.columns:
            raise ValueError("Q_EMPIRICAL column not found in results")

        filtered_df = results_df[results_df['Q_EMPIRICAL'] <= threshold].copy()
        filtered_df = filtered_df.sort_values('Q_EMPIRICAL')

        logging.info(f"Found {len(filtered_df)} genes with Q_EMPIRICAL ≤ {threshold}")

        return filtered_df

    except Exception as e:
        logging.error(f"Error processing OncodriveCLUSTL results: {e}")
        return pd.DataFrame()




def detect_smg_oncodriveclustl(input_file: Union[str, Path], genome_version: Literal["GRCh37", "GRCh38"] = "GRCh37", 
                              threads: int = 4, run_analysis: bool = True, threshold: Optional[float] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Detects significantly mutated genes (SMG) using OncodriveCLUSTL pipeline.

    This function processes a MAF file to generate mutations and regions data, then runs
    OncodriveCLUSTL analysis to identify significantly mutated genes. It includes:
    1. MAF processing and SNP filtering
    2. Transcript regions generation using pyensembl
    3. Chromosome format validation
    4. Optimal permutation calculation based on SNV count
    5. OncodriveCLUSTL execution
    6. Results filtering by Q_EMPIRICAL threshold

    Parameters
    ----------
    input_file : str or Path
        Path to the input MAF file (.maf or .maf.gz)
    genome_version : {"GRCh37", "GRCh38"}, default "GRCh37"
        Genome version to use for transcript annotation.
        - "GRCh37": Uses Ensembl release 75 (hg19)
        - "GRCh38": Uses Ensembl release 114 (hg38)
    threads : int, default 4
        Number of threads to use for OncodriveCLUSTL analysis
    run_analysis : bool, default True
        Whether to run the complete OncodriveCLUSTL analysis pipeline
    threshold : float, optional, default None
        Q_EMPIRICAL threshold for filtering significant genes. If None, 
        process_oncodriveclustl_results will not be executed and an empty 
        DataFrame will be returned for significant genes.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        - mutations_df: DataFrame with columns CHROMOSOME, POSITION, REF, ALT, SAMPLE
        - significant_genes_df: DataFrame with significant genes (Q_EMPIRICAL ≤ threshold)

    Examples
    --------
    >>> mutations_df, results_df = detect_smg_oncodriveclustl("tcga_laml.maf.gz", genome_version="GRCh37")
    >>> print(f"Found {len(results_df)} significant genes")
    """
    input_file = Path(input_file)

    if input_file.suffix == '.gz':
        with gzip.open(input_file, 'rt') as f:
            df = pd.read_csv(f, sep='\t', low_memory=False)
    else:
        df = pd.read_csv(input_file, sep='\t', low_memory=False)

    # Filter for SNP variants only
    variant_type_col = col(df, "Variant_Type", required=True)
    df_snp = df[variant_type_col == "SNP"].copy()

    if df_snp.empty:
        raise ValueError("No SNP variants found in the input file")

    output_df = pd.DataFrame()

    chromosome_col = col(df_snp, "Chromosome", required=True)
    output_df['CHROMOSOME'] = chromosome_col.astype(str).apply(reverse_format_chr)

    position_col = col(df_snp, "Start_Position", required=True)
    output_df['POSITION'] = position_col

    ref_col = col(df_snp, "Reference_Allele", required=True)
    output_df['REF'] = ref_col

    allele2_col = col(df_snp, "Tumor_Seq_Allele2", required=True)
    allele1_col = col(df_snp, "Tumor_Seq_Allele1", required=True)

    output_df['ALT'] = allele2_col.fillna('').replace('', None)
    mask_empty_allele2 = output_df['ALT'].isna() | (output_df['ALT'] == '')
    output_df.loc[mask_empty_allele2, 'ALT'] = allele1_col[mask_empty_allele2]

    sample_col = col(df_snp, "Tumor_Sample_Barcode", required=True)
    output_df['SAMPLE'] = sample_col

    logging.info(f"Generating transcript regions for {genome_version}...")
    regions_df = _generate_transcript_regions(genome_version)

    if input_file.suffix == '.gz':
        base_name = input_file.stem
        base_name = Path(base_name).stem
    else:
         base_name = input_file.stem

    now = datetime.now()
    timestamp = f"_aux_{now.hour:02d}{now.minute:02d}{now.day:02d}{now.month:02d}"
    
    output_folder = input_file.parent / f"{base_name}{timestamp}"
    output_folder.mkdir(exist_ok=True)

    mutations_file = output_folder / "tcga_laml_mutations.tsv"
    output_df.to_csv(mutations_file, sep='\t', index=False)

    # Save regions data with gzip compression in the new folder
    regions_file = output_folder / f"cds.{genome_version.lower()}.regions.gz"
    regions_df.to_csv(
        regions_file,
        sep="\t",
        index=False,
        compression="gzip"
    )

    logging.info(f"Processed {len(df)} total variants")
    logging.info(f"Filtered to {len(df_snp)} SNP variants")
    logging.info(f"Mutations output saved to: {mutations_file}")
    logging.info(f"Regions output saved to: {regions_file}")

    significant_genes_df = pd.DataFrame()

    # Run OncodriveCLUSTL analysis pipeline if requested
    if run_analysis:
        logging.info("Starting OncodriveCLUSTL analysis...")

        if not _validate_chromosome_format(mutations_file, regions_file):
            logging.warning("Chromosome format inconsistency detected, but continuing...")

        snv_count = len(output_df)
        n_permutations = _calculate_optimal_permutations(snv_count)
        logging.info(f"Using {n_permutations} permutations for {snv_count} SNVs")

        genome_build = "hg19" if genome_version == "GRCh37" else "hg38"

        # Run OncodriveCLUSTL
        output_dir = output_folder / f"{base_name}_oncodriveclustl_results"

        success = _run_oncodriveclustl(
            mutations_file=mutations_file,
            regions_file=regions_file,
            genome_build=genome_build,
            output_dir=output_dir,
            n_permutations=n_permutations,
            threads=threads
        )

        if success:
            if threshold is not None:
                significant_genes_df = process_oncodriveclustl_results(results_dir=output_dir, threshold=threshold)

                if not significant_genes_df.empty:
                    logging.info(f"Analysis completed: {len(significant_genes_df)} significant genes found")
                else:
                    logging.warning(f"No significant genes found with Q_EMPIRICAL ≤ {threshold}")
            else:
                logging.info("Analysis completed (results not processed - no threshold provided)")
        else:
            logging.error("OncodriveCLUSTL analysis failed")
    else:
        logging.info("Skipping OncodriveCLUSTL analysis (run_analysis=False)")

    return output_df, significant_genes_df
