import logging
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Tuple, Optional

import pandas as pd

try:
    from pyensembl import EnsemblRelease
except ImportError:
    EnsemblRelease = None

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
        raise ImportError(
            "pyensembl is required for transcript regions generation. Install with: pip install pyensembl")

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
            except Exception:
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

        logging.info("Chromosome formats are consistent between files")
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

        logging.info("Running OncodriveCLUSTL with command:")
        logging.info(' '.join(cmd))

        logger = logging.getLogger("oncodrive")
        logger.setLevel(logging.INFO)
        logger.addHandler(logging.StreamHandler())

        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in proc.stdout:
            logger.info(line.rstrip())
        proc.wait()

        if proc.returncode == 0:
            logging.info("OncodriveCLUSTL completed successfully")
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


class SmgDetectionMixin:
    """
    Mixin class providing significantly mutated genes (SMG) detection functionality for PyMutation objects.
    
    This mixin adds OncodriveCLUSTL-based SMG detection capabilities to PyMutation,
    following the same architectural pattern as other mixins in the project.
    """

    def detect_smg_oncodriveclustl(self, threads: int = 4, run_analysis: bool = True, threshold: Optional[float] = None) -> \
    Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Detects significantly mutated genes (SMG) using OncodriveCLUSTL pipeline.

        This method processes VCF-like mutation data to generate mutations and regions data, then runs
        OncodriveCLUSTL analysis to identify significantly mutated genes. It includes:
        1. VCF-like data processing and SNP filtering
        2. Transcript regions generation using pyensembl
        3. Chromosome format validation
        4. Optimal permutation calculation based on SNV count
        5. OncodriveCLUSTL execution
        6. Results filtering by Q_EMPIRICAL threshold

        Parameters
        ----------
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
        >>> mutations_df, results_df = py_mutation.detect_smg_oncodriveclustl()
        >>> print(f"Found {len(results_df)} significant genes")
        """
        # Use self.data directly (VCF-like format with CHROM, POS, REF, ALT, SAMPLE columns)
        df = self.data.copy()

        # Get genome version from metadata
        if self.metadata is None or self.metadata.assembly is None:
            raise ValueError("Metadata with assembly information is required")

        # Convert assembly format: "37" -> "GRCh37", "38" -> "GRCh38"
        if self.metadata.assembly == "37":
            genome_version = "GRCh37"
        elif self.metadata.assembly == "38":
            genome_version = "GRCh38"
        else:
            raise ValueError(f"Unsupported assembly version: {self.metadata.assembly}. Expected '37' or '38'")

        # Check if we have MAF format first (has Tumor_Sample_Barcode)
        maf_cols = ['CHROM', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode']
        has_maf_format = all(col in df.columns for col in maf_cols)

        # Check if we have basic VCF-like columns (CHROM, POS, REF, ALT)
        vcf_basic_cols = ['CHROM', 'POS', 'REF', 'ALT']
        has_vcf_basic = all(col in df.columns for col in vcf_basic_cols)

        if has_maf_format:
            # Convert MAF format to VCF-like format
            logging.info("Converting MAF format to VCF-like format...")
            df = df.copy()
            # CHROM is already present in MAF data, no need to copy
            df['POS'] = df['Start_Position']
            df['REF'] = df['Reference_Allele']
            # Use Tumor_Seq_Allele2, fallback to Tumor_Seq_Allele1 if empty
            df['ALT'] = df['Tumor_Seq_Allele2'].fillna('')
            if 'Tumor_Seq_Allele1' in df.columns:
                mask_empty = (df['ALT'] == '') | df['ALT'].isna()
                df.loc[mask_empty, 'ALT'] = df.loc[mask_empty, 'Tumor_Seq_Allele1']
            df['SAMPLE'] = df['Tumor_Sample_Barcode']
            is_vcf_like = False
        elif has_vcf_basic:
            # VCF-like format with samples in separate columns
            logging.info("Processing VCF-like format with sample columns...")
            is_vcf_like = True
        else:
            # Neither format found
            missing_maf = [col for col in maf_cols if col not in df.columns]
            missing_vcf = [col for col in vcf_basic_cols if col not in df.columns]
            raise ValueError(
                f"Data format not recognized. Missing MAF columns: {missing_maf} or VCF basic columns: {missing_vcf}")

        # Filter for SNP variants only (REF and ALT should be single nucleotides)
        df_snp = df[
            (df['REF'].str.len() == 1) &
            (df['ALT'].str.len() == 1) &
            (df['REF'] != df['ALT'])
            ].copy()

        if df_snp.empty:
            raise ValueError("No SNP variants found in the data")

        # Process data to create mutations DataFrame
        output_df = pd.DataFrame()
        expanded_rows = []

        if is_vcf_like:
            # Handle VCF-like format with samples in separate columns
            # Get sample columns from self.samples or detect them
            if hasattr(self, 'samples') and self.samples:
                sample_columns = self.samples
            else:
                # Detect sample columns (exclude standard VCF columns)
                standard_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
                sample_columns = [col for col in df_snp.columns if col not in standard_cols]

            logging.info(f"Processing {len(sample_columns)} sample columns: {sample_columns[:5]}...")

            for _, row in df_snp.iterrows():
                for sample_col in sample_columns:
                    genotype = row[sample_col]

                    # Check if this sample has a mutation (non-reference genotype)
                    if pd.notna(genotype) and genotype not in ['0/0', '0|0', './.', '.|.', '']:
                        # Parse genotype to check if it contains the ALT allele
                        if isinstance(genotype, str):
                            # Handle various genotype formats
                            if ('1' in genotype or  # 0/1, 1/0, 1/1, 1|1, etc.
                                    '/' in genotype or '|' in genotype or  # Any phased/unphased genotype
                                    genotype == row['ALT'] or  # Direct ALT match
                                    (len(genotype) > 1 and row['ALT'] in genotype)):  # ALT in compound genotype

                                expanded_rows.append({
                                    'CHROMOSOME': str(row['CHROM']).replace('chr', ''),
                                    'POSITION': row['POS'],
                                    'REF': row['REF'],
                                    'ALT': row['ALT'],
                                    'SAMPLE': sample_col
                                })
        else:
            # Handle MAF format (already has SAMPLE column)
            for _, row in df_snp.iterrows():
                expanded_rows.append({
                    'CHROMOSOME': str(row['CHROM']).replace('chr', ''),
                    'POSITION': row['POS'],
                    'REF': row['REF'],
                    'ALT': row['ALT'],
                    'SAMPLE': row['SAMPLE']
                })

        output_df = pd.DataFrame(expanded_rows)

        # Apply reverse_format_chr to ensure consistent chromosome format
        output_df['CHROMOSOME'] = output_df['CHROMOSOME'].astype(str).apply(reverse_format_chr)

        logging.info(f"Generating transcript regions for {genome_version}...")
        regions_df = _generate_transcript_regions(genome_version)

        # Create temporary directory for analysis files
        now = datetime.now()
        timestamp = f"_aux_{now.hour:02d}{now.minute:02d}{now.day:02d}{now.month:02d}"
        base_name = "pymutation_data"

        output_folder = Path.cwd() / f"{base_name}{timestamp}"
        output_folder.mkdir(exist_ok=True)

        mutations_file = output_folder / "mutations.tsv"
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
