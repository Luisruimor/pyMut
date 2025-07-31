"""
PyMutation core module.

Este módulo contiene la clase principal PyMutation que sirve como API principal
para la librería pyMut. Proporciona métodos para generar todos los tipos de
visualizaciones a partir de datos de mutación.
"""

from datetime import datetime
from typing import List
from typing import Tuple, Optional

import matplotlib.pyplot as plt
import pandas as pd

from .analysis.pfam_annotation import PfamAnnotationMixin
from .analysis.mutation_burden import MutationBurdenMixin
from .analysis.mutational_signature import MutationalSignatureMixin
from .analysis.smg_detection import SmgDetectionMixin
from .annotate.actionable_mutation import ActionableMutationMixin
from .annotate.cosmic_cancer_annotate import CancerAnnotateMixin
from .output import OutputMixin
from .filters.chrom_sample_filter import ChromSampleFilterMixin
from .filters.genomic_range import GenomicRangeMixin
from .filters.pass_filter import PassFilterMixin
from .filters.tissue_expression import TissueExpressionMixin
from .utils.constants import (DEFAULT_PLOT_FIGSIZE, DEFAULT_SUMMARY_FIGSIZE, DEFAULT_PLOT_TITLE,
                              DEFAULT_TOP_GENES_COUNT, GENE_COLUMN, VARIANT_CLASSIFICATION_COLUMN, SAMPLE_COLUMN,
                              REF_COLUMN, ALT_COLUMN, MODE_VARIANTS, DEFAULT_ONCOPLOT_FIGSIZE,
                              DEFAULT_ONCOPLOT_TOP_GENES, DEFAULT_ONCOPLOT_MAX_SAMPLES, FUNCOTATION_COLUMN,
                              VARIANT_TYPE_COLUMN, VALID_PLOT_MODES)


class MutationMetadata:
    """
    Clase para almacenar metadatos de mutaciones.

    Atributos:
        source_format (str): Formato de origen (VCF, MAF, etc.).
        file_path (str): Ruta del archivo de origen.
        loaded_at (datetime): Fecha y hora de carga.
        filters (List[str]): Filtros aplicados al archivo.
        assembly (str): Versión del genoma (37 o 38).
        notes (Optional[str]): Notas adicionales.
    """

    def __init__(self, source_format: str, file_path: str, filters: List[str], assembly: str,
                 notes: Optional[str] = None):
        self.source_format = source_format
        self.file_path = file_path
        self.loaded_at = datetime.now()
        self.filters = filters
        self.notes = notes
        self.assembly = assembly


class PyMutation(CancerAnnotateMixin, ActionableMutationMixin, MutationBurdenMixin, MutationalSignatureMixin, PfamAnnotationMixin, SmgDetectionMixin, OutputMixin, ChromSampleFilterMixin, GenomicRangeMixin, PassFilterMixin, TissueExpressionMixin):
    def __init__(self, data: pd.DataFrame, metadata: Optional[MutationMetadata] = None,
                 samples: Optional[List[str]] = None):
        self.data = data
        self.samples = samples if samples is not None else []
        self.metadata = metadata

    def head(self, n: int = 5):
        """
        Return the first n rows of the mutation data.
        This method delegates to the pandas DataFrame head() method.
        """
        return self.data.head(n)

    def info(self):
        """
        Print a concise summary of the mutation data.
        This method delegates to the pandas DataFrame info() method
        """
        return self.data.info()

    def save_figure(self, figure: plt.Figure, filename: str, dpi: int = 300, bbox_inches: str = 'tight',
                    **kwargs) -> None:
        """
        Save a figure with high-quality configuration by default.

        This method centralizes figure saving to ensure all visualizations
        are saved with the best possible quality.

        Args:
            figure: The matplotlib figure to save.
            filename: Filename where to save the figure.
            dpi: Resolution in dots per inch (300 = high quality).
            bbox_inches: Margin adjustment ('tight' = no unnecessary spaces).
            **kwargs: Additional parameters for matplotlib.savefig().

        Examples:
            >>> py_mut = PyMutation(data)
            >>> fig = py_mut.summary_plot()
            >>> py_mut.save_figure(fig, 'my_summary.png')  # Automatic high quality
            >>> py_mut.save_figure(fig, 'my_summary.pdf', dpi=600)  # Very high quality
        """
        figure.savefig(filename, dpi=dpi, bbox_inches=bbox_inches, **kwargs)
        print(f"📁 Figure saved: {filename} (DPI: {dpi}, margins: {bbox_inches})")

    @staticmethod
    def configure_high_quality_plots():
        """
        Configure matplotlib to generate high-quality plots by default.

        This function modifies matplotlib's global configuration so that
        ALL figures are automatically saved with high quality, without
        needing to specify parameters each time.

        Applied configurations:
        - DPI: 300 (high resolution)
        - bbox_inches: 'tight' (optimized margins)
        - Format: PNG with optimized compression

        Examples:
            >>> PyMutation.configure_high_quality_plots()  # Configure once
            >>> py_mut = PyMutation(data)
            >>> fig = py_mut.summary_plot()
            >>> fig.savefig('plot.png')  # Automatically high quality!

        Note:
            This configuration affects ALL matplotlib figures in the session.
            It's recommended to call this function at the beginning of the script.
        """
        import matplotlib as mpl

        # Configure default DPI for high resolution
        mpl.rcParams['figure.dpi'] = 300
        mpl.rcParams['savefig.dpi'] = 300

        # Configure automatic margins
        mpl.rcParams['savefig.bbox'] = 'tight'

        # Configure format and compression
        mpl.rcParams['savefig.format'] = 'png'
        mpl.rcParams['savefig.transparent'] = False

        # Improve text quality
        mpl.rcParams['savefig.facecolor'] = 'white'
        mpl.rcParams['savefig.edgecolor'] = 'none'

    def summary_plot(self, figsize: Tuple[int, int] = DEFAULT_SUMMARY_FIGSIZE, title: str = DEFAULT_PLOT_TITLE,
                     max_samples: Optional[int] = 200, top_genes_count: int = DEFAULT_TOP_GENES_COUNT) -> plt.Figure:
        """
        Generate a summary plot with general mutation statistics.

        This visualization includes multiple plots:
        - Variant Classification: Distribution of variant classifications
        - Variant Type: Distribution of variant types (SNP, INS, DEL, etc.)
        - SNV Class: Distribution of SNV classes (nucleotide changes like A>G, C>T, etc.)
        - Variants per Sample: Distribution of variants per sample and median (TMB)
        - Top Mutated Genes: Most frequently mutated genes

        Args:
            figsize: Figure size.
            title: Plot title.
            max_samples: Maximum number of samples to show in the variants per sample plot.
                        If None, all samples are shown.
            top_genes_count: Number of genes to show in the top mutated genes plot.
                        If there are fewer genes than this number, all will be shown.

        Returns:
            Matplotlib figure with the summary plot.
        """
        from .visualizations.summary import _create_summary_plot
        from .utils.data_processing import extract_variant_classifications, extract_variant_types

        # Preprocess data to ensure we have the necessary columns
        self.data = extract_variant_classifications(self.data, variant_column=VARIANT_CLASSIFICATION_COLUMN,
                                                    funcotation_column=FUNCOTATION_COLUMN)

        self.data = extract_variant_types(self.data, variant_column=VARIANT_TYPE_COLUMN,
                                          funcotation_column=FUNCOTATION_COLUMN)

        # Generate the summary plot
        fig = _create_summary_plot(self, figsize, title, max_samples, top_genes_count)

        plt.close(fig)
        return fig

    def variant_classification_plot(self, figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                                    title: str = "Variant Classification") -> plt.Figure:
        """
        Generate a horizontal bar plot showing the distribution of variant classifications.

        Args:
            figsize: Figure size.
            title: Plot title.

        Returns:
            Matplotlib figure with the variant classification plot.
        """
        from .visualizations.summary import _create_variant_classification_plot
        from .utils.data_processing import extract_variant_classifications

        # Preprocess data to ensure we have the necessary column
        self.data = extract_variant_classifications(self.data, variant_column="Variant_Classification",
                                                    funcotation_column="FUNCOTATION")

        # Create figure and axes
        fig, ax = plt.subplots(figsize=figsize)

        # Generate the plot, passing set_title=False to avoid duplicate title
        _create_variant_classification_plot(self, ax=ax, set_title=False)

        # Configure title
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')

        plt.tight_layout()

        plt.close(fig)
        return fig

    def variant_type_plot(self, figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                          title: str = "Variant Type") -> plt.Figure:
        """
        Generate a horizontal bar plot showing the distribution of variant types.

        Args:
            figsize: Figure size.
            title: Plot title.

        Returns:
            Matplotlib figure with the variant types plot.
        """
        from .visualizations.summary import _create_variant_type_plot
        from .utils.data_processing import extract_variant_types

        # Preprocess data to ensure we have the necessary column
        self.data = extract_variant_types(self.data, variant_column="Variant_Type", funcotation_column="FUNCOTATION")

        # Create figure and axes
        fig, ax = plt.subplots(figsize=figsize)

        # Generate the plot, passing set_title=False to avoid duplicate title
        _create_variant_type_plot(self, ax=ax, set_title=False)

        # Configure title
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')

        plt.tight_layout()

        plt.close(fig)
        return fig

    def snv_class_plot(self, figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE, title: str = "SNV Class",
                       ref_column: str = "REF", alt_column: str = "ALT") -> plt.Figure:
        """
        Generate a horizontal bar plot showing the distribution of SNV classes.

        Args:
            figsize: Figure size.
            title: Plot title.
            ref_column: Name of the column containing the reference allele.
            alt_column: Name of the column containing the alternative allele.

        Returns:
            Matplotlib figure with the SNV classes plot.
        """
        from .visualizations.summary import _create_snv_class_plot

        fig, ax = plt.subplots(figsize=figsize)
        _create_snv_class_plot(self, ref_column=ref_column, alt_column=alt_column, ax=ax, set_title=False
                               # Avoid duplicate title
                               )

        # Configure title
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')

        plt.tight_layout()

        plt.close(fig)
        return fig

    def variants_per_sample_plot(self, figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                                 title: str = "Variants per Sample", variant_column: str = "Variant_Classification",
                                 sample_column: str = "Tumor_Sample_Barcode",
                                 max_samples: Optional[int] = 200) -> plt.Figure:
        """
        Generate a stacked bar plot showing the number of variants per sample (TMB)
        and their composition by variant type.

        Args:
            figsize: Figure size.
            title: Plot title.
            variant_column: Name of the column containing the variant classification.
            sample_column: Name of the column containing the sample identifier,
                          or string used to identify sample columns if samples
                          are stored as columns.
            max_samples: Maximum number of samples to show. If None, all are shown.
                        If there are more samples than this number, only the first max_samples are shown.

        Returns:
            Matplotlib figure with the variants per sample plot.
        """
        from .visualizations.summary import _create_variants_per_sample_plot
        from .utils.data_processing import extract_variant_classifications

        # If variant_column is not in columns, try to normalize it
        if variant_column not in self.data.columns:
            # Check if there's a version with different capitalization
            column_lower = variant_column.lower()
            for col in self.data.columns:
                if col.lower() == column_lower:
                    variant_column = col
                    break

        # Ensure the variant classification column exists or is extracted
        self.data = extract_variant_classifications(self.data, variant_column=variant_column,
                                                    funcotation_column="FUNCOTATION")

        fig, ax = plt.subplots(figsize=figsize)
        _create_variants_per_sample_plot(self, variant_column=variant_column, sample_column=sample_column, ax=ax,
                                         set_title=False,  # Avoid duplicate title
                                         max_samples=max_samples  # Pass the configured sample limit
                                         )

        # Don't modify title if it contains the median
        if title and not title.startswith("Variants per Sample"):
            fig.suptitle(title, fontsize=16, fontweight='bold', y=1.02)
        elif title:
            fig.suptitle(title, fontsize=16, fontweight='bold')

        plt.tight_layout()

        plt.close(fig)
        return fig

    def variant_classification_summary_plot(self, figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                                            title: str = "Variant Classification Summary",
                                            variant_column: str = "Variant_Classification",
                                            sample_column: str = "Tumor_Sample_Barcode") -> plt.Figure:
        """
        Generate a box-and-whiskers plot (boxplot) that summarizes, for each variant classification,
        the distribution (among samples) of the number of detected alternative alleles.

        This plot shows the variability between samples for each type of variant classification,
        allowing identification of which ones present more differences between patients.

        Args:
            figsize: Figure size.
            title: Plot title.
            variant_column: Name of the column containing the variant classification.
            sample_column: Name of the column containing the sample identifier.
                           If it doesn't exist, samples are assumed to be columns (wide format).

        Returns:
            Matplotlib figure with the box-and-whiskers plot.
        """
        from .visualizations.summary import _create_variant_classification_summary_plot
        from .utils.data_processing import extract_variant_classifications

        # Ensure the variant classification column exists or is extracted
        self.data = extract_variant_classifications(self.data, variant_column=variant_column,
                                                    funcotation_column="FUNCOTATION")

        # Check if we're in wide format (samples as columns)
        is_wide_format = sample_column not in self.data.columns
        if is_wide_format:
            # Detect and show information about the format
            sample_cols = [col for col in self.data.columns if
                           col.startswith('TCGA-') or (isinstance(col, str) and col.count('-') >= 2)]
            if sample_cols:
                print(f"Detected wide format with {len(sample_cols)} possible sample columns.")

        fig, ax = plt.subplots(figsize=figsize)
        _create_variant_classification_summary_plot(self, variant_column=variant_column, sample_column=sample_column,
                                                    ax=ax, show_labels=True, set_title=False)

        # Configure title
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')

        plt.tight_layout()

        plt.close(fig)
        return fig

    def top_mutated_genes_plot(self, figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE, title: str = "Top Mutated Genes",
                               mode: str = MODE_VARIANTS, variant_column: str = VARIANT_CLASSIFICATION_COLUMN,
                               gene_column: str = GENE_COLUMN, sample_column: str = SAMPLE_COLUMN,
                               count: int = DEFAULT_TOP_GENES_COUNT) -> plt.Figure:
        """
        Generate a horizontal bar plot showing the most mutated genes and the distribution
        of variants according to their classification.

        Args:
            figsize: Figure size.
            title: Plot title.
            mode: Mutation counting mode: "variants" (counts total number of variants)
                  or "samples" (counts number of affected samples).
            variant_column: Name of the column containing the variant classification.
            gene_column: Name of the column containing the gene symbol.
            sample_column: Name of the column containing the sample identifier,
                          or prefix to identify sample columns if they are columns.
            count: Number of top genes to show.

        Returns:
            Matplotlib figure with the top mutated genes plot.

        Raises:
            ValueError: If 'count' is not a positive number or 'mode' is not a valid value.
        """
        from .visualizations.summary import _create_top_mutated_genes_plot
        from .utils.data_processing import extract_variant_classifications

        # Validate parameters
        if not isinstance(count, int):
            raise ValueError(f"The 'count' parameter must be an integer, received: {count}")
        if count <= 0:
            raise ValueError(f"The 'count' parameter must be a positive integer, received: {count}")

        # Check that mode is valid
        if mode not in VALID_PLOT_MODES:
            raise ValueError(f"Mode '{mode}' is not valid. Allowed values are: {', '.join(VALID_PLOT_MODES)}")

        # If variant_column is not in columns, try to normalize it
        if variant_column not in self.data.columns:
            # Check if there's a version with different capitalization
            column_lower = variant_column.lower()
            for col in self.data.columns:
                if col.lower() == column_lower:
                    variant_column = col
                    break

        # If gene_column is not in columns, try to normalize it
        if gene_column not in self.data.columns:
            # Check if there's a version with different capitalization
            column_lower = gene_column.lower()
            for col in self.data.columns:
                if col.lower() == column_lower:
                    gene_column = col
                    break

        # Ensure the variant classification column exists or is extracted
        self.data = extract_variant_classifications(self.data, variant_column=variant_column,
                                                    funcotation_column=FUNCOTATION_COLUMN)

        fig, ax = plt.subplots(figsize=figsize)
        _create_top_mutated_genes_plot(self, mode=mode, variant_column=variant_column, gene_column=gene_column,
                                       sample_column=sample_column, count=count, ax=ax, set_title=False
                                       # Avoid duplicate title
                                       )

        # Adjust custom title based on mode
        if title:
            if mode == "variants" and title == "Top Mutated Genes":
                fig.suptitle("Top mutated genes (variants)", fontsize=16, fontweight='bold', y=0.98)
            elif mode == "samples" and title == "Top Mutated Genes":
                fig.suptitle("Top mutated genes (samples)", fontsize=16, fontweight='bold', y=0.98)
            else:
                fig.suptitle(title, fontsize=16, fontweight='bold', y=0.98)

        # Use tight_layout with additional padding to improve margins
        plt.tight_layout(pad=1.2)

        # Adjust margins for more consistent appearance between modes
        # Increase left margin to prevent text from being cut off
        plt.subplots_adjust(left=0.15, right=0.9)

        plt.close(fig)
        return fig

    def oncoplot(self, figsize: Optional[Tuple[int, int]] = None, title: str = "Oncoplot",
                 gene_column: str = GENE_COLUMN, variant_column: str = VARIANT_CLASSIFICATION_COLUMN,
                 ref_column: str = REF_COLUMN, alt_column: str = ALT_COLUMN, top_genes_count: int = None,
                 max_samples: int = None) -> plt.Figure:
        """
        Generates an oncoplot showing mutation patterns in a heatmap.
        
        The oncoplot is a fundamental visualization in cancer genomics that shows
        mutation patterns across samples and genes in heatmap format.
        
        Features:
        - Automatic detection of sample columns (TCGA and .GT format)
        - Support for multiple genotype formats (A|G, A/G, etc.)
        - Multi_Hit detection for samples with multiple mutations
        - Standard color schemes for mutation types
        - Smart gene ordering by mutation frequency
        - Sample ordering by mutational burden
        
        Args:
            figsize: Figure size (width, height) in inches.
                    If None, uses DEFAULT_ONCOPLOT_FIGSIZE.
            title: Title for the visualization.
            gene_column: Name of the column containing gene symbols.
            variant_column: Name of the column containing variant classifications.
            ref_column: Name of the column containing reference alleles.
            alt_column: Name of the column containing alternative alleles.
            top_genes_count: Number of top mutated genes to show.
                           If None, uses DEFAULT_ONCOPLOT_TOP_GENES.
            max_samples: Maximum number of samples to show.
                        If None, uses DEFAULT_ONCOPLOT_MAX_SAMPLES.
            
        Returns:
            plt.Figure: matplotlib Figure object with the oncoplot.
            
        Raises:
            ValueError: If required columns are missing, no mutation data,
                       or problems with data format.
            
        Examples:
            Basic usage:
            >>> py_mut = PyMutation(data)
            >>> fig = py_mut.oncoplot()
            >>> fig.savefig('oncoplot.png')
            
            With custom parameters:
            >>> fig = py_mut.oncoplot(
            ...     title="TCGA Samples Oncoplot",
            ...     top_genes_count=20,
            ...     max_samples=100
            ... )
            
        Note:
            - The method automatically detects sample columns using common
              patterns like 'TCGA-*' and '*.GT'
            - Genes are ordered by mutation frequency (most mutated at top)
            - Samples are ordered by total mutational burden
            - Colors follow cancer genomics standards
            - Multi_Hit detection is handled automatically
        """
        # Validate input parameters
        if top_genes_count is None:
            top_genes_count = DEFAULT_ONCOPLOT_TOP_GENES
        if max_samples is None:
            max_samples = DEFAULT_ONCOPLOT_MAX_SAMPLES
        if figsize is None:
            figsize = DEFAULT_ONCOPLOT_FIGSIZE

        # Parameter validation
        if top_genes_count <= 0:
            raise ValueError("top_genes_count must be a positive integer")
        if max_samples <= 0:
            raise ValueError("max_samples must be a positive integer")
        if len(figsize) != 2 or any(x <= 0 for x in figsize):
            raise ValueError("figsize must be a tuple of two positive numbers")

        # Validate required columns
        required_columns = [gene_column, variant_column, ref_column, alt_column]
        missing_columns = [col for col in required_columns if col not in self.data.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")

        try:
            from .visualizations.oncoplot import _create_oncoplot_plot
            # Generate the oncoplot
            fig = _create_oncoplot_plot(py_mut=self, gene_column=gene_column, variant_column=variant_column,
                                        ref_column=ref_column, alt_column=alt_column, top_genes_count=top_genes_count,
                                        max_samples=max_samples,
                                        figsize=figsize, title=title)
            plt.close(fig)
            return fig

        except Exception as e:
            raise ValueError(f"Error generating oncoplot: {str(e)}")


