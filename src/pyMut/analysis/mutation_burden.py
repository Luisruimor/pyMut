import logging
import os
import re
from typing import Optional, Dict

import pandas as pd

# Logger configuration
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Change to DEBUG for more verbosity
if not logger.handlers:
    logging.basicConfig(
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        level=logging.INFO,
    )


class MutationBurdenMixin:
    """
    Mixin class providing mutation burden analysis functionality for PyMutation objects.
    
    This mixin adds TMB (Tumor Mutation Burden) analysis capabilities to PyMutation,
    following the same architectural pattern as other mixins in the project.
    """

    def calculate_tmb_analysis(self,
                               variant_classification_column: Optional[str] = None,
                               genome_size_bp: int = 60456963,
                               output_dir: str = ".",
                               save_files: bool = True) -> Dict[str, pd.DataFrame]:
        """
        Calculate Tumor Mutation Burden (TMB) analysis for each sample in a PyMutation object.

        This method analyzes mutation data from PyMutation objects to calculate TMB metrics 
        for each sample, including total mutations, non-synonymous mutations, and normalized 
        TMB values. It's designed to work with the VCF-like structure that PyMutation uses.

        Parameters
        ----------
        variant_classification_column : str, optional
            Name of the column containing variant classification information.
            If None, will automatically detect variant classification columns.
        genome_size_bp : int, default 60456963
            Size of the interrogated region in base pairs for exact TMB normalization.
            Default 60,456,963 bp is the standard for Whole Exome Sequencing (WES).
            Use approximately 3,000,000,000 bp for Whole Genome Sequencing (WGS).
        output_dir : str, default "."
            Directory where output files will be saved.
        save_files : bool, default True
            Whether to save the results to TSV files.

        Returns
        -------
        Dict[str, pd.DataFrame]
            Dictionary containing:
            - 'analysis': Per-sample TMB analysis DataFrame
            - 'statistics': Global TMB statistics DataFrame

        Notes
        -----
        Non-synonymous mutations are defined as those with variant classifications
        that typically have biological impact, including:
        - MISSENSE_MUTATION, NONSENSE_MUTATION, FRAME_SHIFT_DEL, FRAME_SHIFT_INS
        - NONSTOP_MUTATION, TRANSLATION_START_SITE, SPLICE_SITE

        The method generates two output files:
        - TMB_analysis.tsv: Per-sample analysis with mutation counts and normalized TMB
        - TMB_statistics.tsv: Global statistics (mean, median, quartiles, etc.)
        """

        # Define non-synonymous mutation types (more comprehensive list)
        non_synonymous_types = {
            'MISSENSE_MUTATION', 'NONSENSE_MUTATION', 'FRAME_SHIFT_DEL', 'FRAME_SHIFT_INS',
            'NONSTOP_MUTATION', 'TRANSLATION_START_SITE', 'SPLICE_SITE', 'IN_FRAME_DEL',
            'IN_FRAME_INS', 'START_CODON_SNP', 'START_CODON_DEL', 'START_CODON_INS',
            'STOP_CODON_DEL', 'STOP_CODON_INS', 'DE_NOVO_START_IN_FRAME',
            'DE_NOVO_START_OUT_FRAME'
        }

        # Validate PyMutation structure
        if not hasattr(self, 'data') or not hasattr(self, 'samples'):
            raise ValueError("Invalid PyMutation object: missing 'data' or 'samples' attributes")

        if self.data.empty:
            raise ValueError("PyMutation data is empty")

        if not self.samples:
            raise ValueError("No samples found in PyMutation object")

        # Check for required columns
        required_cols = ['REF', 'ALT']
        missing_cols = [col for col in required_cols if col not in self.data.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in PyMutation data: {missing_cols}")

        # Auto-detect variant classification column if not provided
        if variant_classification_column is None:
            # Look for variant classification columns using regex pattern
            pattern = re.compile(r'^(gencode_\d+_)?variant[_]?classification$', flags=re.IGNORECASE)
            variant_cols = [col for col in self.data.columns if pattern.match(col)]

            if variant_cols:
                variant_classification_column = variant_cols[0]
                logger.info(f"Auto-detected variant classification column: {variant_classification_column}")
            else:
                logger.warning("No variant classification column found. Non-synonymous counts will be 0.")
                variant_classification_column = None

        # Validate variant classification column if provided
        if variant_classification_column and variant_classification_column not in self.data.columns:
            raise ValueError(f"Column provide '{variant_classification_column}' not found in data")

        # Initialize results
        results = []

        # Process each sample
        for sample in self.samples:
            if sample not in self.data.columns:
                logger.warning(f"Sample '{sample}' not found in data columns. Skipping.")
                continue

            # Find mutations for this sample, when the genotype is not REF|REF
            sample_genotypes = self.data[sample].astype(str)

            # Handle missing values and invalid genotypes
            valid_genotypes = sample_genotypes.notna() & (sample_genotypes != '') & (sample_genotypes != '.')

            # Create mutation mask by checking if genotype contains ALT allele
            mutation_mask = valid_genotypes.copy()

            for idx in valid_genotypes[valid_genotypes].index:
                genotype = sample_genotypes.loc[idx]
                ref_allele = str(self.data.loc[idx, 'REF'])
                alt_allele = str(self.data.loc[idx, 'ALT'])

                # Check if genotype contains the ALT allele
                if '|' in genotype:
                    alleles = genotype.split('|')
                elif '/' in genotype:
                    alleles = genotype.split('/')
                else:
                    # Single allele or malformed genotype
                    alleles = [genotype]

                # A mutation exists if any allele matches ALT or is not REF
                has_mutation = False
                for allele in alleles:
                    allele = allele.strip()
                    if allele == alt_allele or (allele != ref_allele and allele != '.' and allele != ''):
                        has_mutation = True
                        break

                mutation_mask.loc[idx] = has_mutation

            # Get mutations for this sample
            sample_mutations = self.data[mutation_mask]
            total_mutations = len(sample_mutations)

            # Count non-synonymous mutations
            non_synonymous_mutations = 0
            if variant_classification_column and not sample_mutations.empty:
                if variant_classification_column in sample_mutations.columns:
                    # Handle missing values in variant classification
                    valid_classifications = sample_mutations[variant_classification_column].notna()
                    valid_sample_mutations = sample_mutations[valid_classifications]

                    if not valid_sample_mutations.empty:
                        # Convert to uppercase and check against non-synonymous types
                        classifications = valid_sample_mutations[variant_classification_column].astype(str).str.upper()
                        non_synonymous_mask = classifications.isin(non_synonymous_types)
                        non_synonymous_mutations = non_synonymous_mask.sum()

            # Calculate normalized TMB (mutations per million bases)
            # TMB = (mutations / genome_size_bp) * 1,000,000
            tmb_total_normalized = (total_mutations / genome_size_bp) * 1_000_000 if genome_size_bp > 0 else 0
            tmb_non_synonymous_normalized = (
                                                    non_synonymous_mutations / genome_size_bp) * 1_000_000 if genome_size_bp > 0 else 0

            # Store results (keep full precision until final output)
            results.append({
                'Sample': sample,
                'Total_Mutations': total_mutations,
                'Non_Synonymous_Mutations': non_synonymous_mutations,
                'TMB_Total_Normalized': tmb_total_normalized,
                'TMB_Non_Synonymous_Normalized': tmb_non_synonymous_normalized
            })

        # Create analysis DataFrame
        if not results:
            raise ValueError("No valid samples found for TMB analysis")

        analysis_df = pd.DataFrame(results)

        # Calculate global statistics
        stats_data = []
        metrics = ['Total_Mutations', 'Non_Synonymous_Mutations',
                   'TMB_Total_Normalized', 'TMB_Non_Synonymous_Normalized']

        for metric in metrics:
            if metric in analysis_df.columns:
                values = analysis_df[metric]
                stats_data.append({
                    'Metric': metric,
                    'Count': len(values),
                    'Mean': values.mean(),
                    'Median': values.median(),
                    'Min': values.min(),
                    'Max': values.max(),
                    'Q1': values.quantile(0.25),
                    'Q3': values.quantile(0.75),
                    'Std': values.std()
                })

        statistics_df = pd.DataFrame(stats_data)

        # Create output directory if it doesn't exist
        if save_files:
            os.makedirs(output_dir, exist_ok=True)

            analysis_path = os.path.join(output_dir, "TMB_analysis.tsv")
            statistics_path = os.path.join(output_dir, "TMB_statistics.tsv")

            # Save files with proper formatting (maintain precision)
            analysis_df.to_csv(analysis_path, sep='\t', index=False, float_format='%.6f')
            statistics_df.to_csv(statistics_path, sep='\t', index=False, float_format='%.6f')

            logger.info(f"TMB analysis saved to: {analysis_path}")
            logger.info(f"TMB statistics saved to: {statistics_path}")
            logger.info(f"Analyzed {len(analysis_df)} samples with {len(self.data)} total mutations")

        log_tmb_summary(analysis_df)

        return {
            'analysis': analysis_df,
            'statistics': statistics_df
        }


def log_tmb_summary(analysis_df: pd.DataFrame) -> None:
    """
    Log a comprehensive summary of TMB analysis results.

    This function provides a detailed summary of the TMB analysis using proper logging
    instead of print statements. It displays key insights about the mutation burden
    across all samples in the analysis.

    Parameters
    ----------
    analysis_df : pd.DataFrame
        The analysis DataFrame returned by calculate_tmb_analysis containing
        per-sample TMB metrics.
    """
    if analysis_df.empty:
        logger.warning("No analysis data provided for summary")
        return

    logger.info("TMB ANALYSIS SUMMARY")

    # Basic statistics
    total_samples = len(analysis_df)
    avg_total_mutations = analysis_df['Total_Mutations'].mean()
    avg_non_synonymous = analysis_df['Non_Synonymous_Mutations'].mean()
    avg_tmb_total = analysis_df['TMB_Total_Normalized'].mean()
    avg_tmb_non_synonymous = analysis_df['TMB_Non_Synonymous_Normalized'].mean()

    logger.info(f"• Total samples analyzed: {total_samples}")
    logger.info(f"• Average total mutations per sample: {avg_total_mutations:.1f}")
    logger.info(f"• Average non-synonymous mutations per sample: {avg_non_synonymous:.1f}")
    logger.info(f"• Average normalized TMB (total): {avg_tmb_total:.6f} mutations/Mb")
    logger.info(f"• Average normalized TMB (non-synonymous): {avg_tmb_non_synonymous:.6f} mutations/Mb")

    # Extreme values
    max_tmb_idx = analysis_df['TMB_Total_Normalized'].idxmax()
    min_tmb_idx = analysis_df['TMB_Total_Normalized'].idxmin()

    max_tmb_sample = analysis_df.loc[max_tmb_idx, 'Sample']
    max_tmb_value = analysis_df['TMB_Total_Normalized'].max()

    min_tmb_sample = analysis_df.loc[min_tmb_idx, 'Sample']
    min_tmb_value = analysis_df['TMB_Total_Normalized'].min()

    logger.info(f"• Sample with highest TMB: {max_tmb_sample}")
    logger.info(f"  - TMB value: {max_tmb_value:.6f} mutations/Mb")
    logger.info(f"• Sample with lowest TMB: {min_tmb_sample}")
    logger.info(f"  - TMB value: {min_tmb_value:.6f} mutations/Mb")
