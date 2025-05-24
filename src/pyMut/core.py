import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Dict, Union, Optional, Tuple
from .utils.constants import (
    VARIANT_CLASSIFICATION_COLUMN, VARIANT_TYPE_COLUMN, SAMPLE_COLUMN, 
    GENE_COLUMN, REF_COLUMN, ALT_COLUMN, FUNCOTATION_COLUMN, 
    DEFAULT_SUMMARY_FIGSIZE, DEFAULT_PLOT_FIGSIZE, DEFAULT_PLOT_TITLE,
    DEFAULT_TOP_GENES_COUNT, MODE_VARIANTS, VALID_PLOT_MODES
)

class PyMutation:
    """
    Main class for visualizing genetic mutations from data in TSV format.
    
    This class provides methods for generating summary visualizations of genetic
    mutation data, showing general statistics such as distributions of variant
    types, classifications, and nucleotide changes.
    
    Attributes:
        data (pd.DataFrame): DataFrame containing mutation data.
    """
    
    def __init__(self, data: pd.DataFrame):
        """
        Initialize a PyMutation object with a pandas DataFrame.
        
        Args:
            data (pd.DataFrame): DataFrame containing mutation data.
        
        Raises:
            ValueError: If the DataFrame is empty or not a valid DataFrame.
        """
        if not isinstance(data, pd.DataFrame):
            raise ValueError("The 'data' parameter must be a pandas DataFrame.")
        
        if data.empty:
            raise ValueError("The provided DataFrame is empty. No data to analyze.")
        
        self.data = data
    
    def save_figure(self, figure: plt.Figure, filename: str, 
                   dpi: int = 300, bbox_inches: str = 'tight', **kwargs) -> None:
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
        
        print("✅ High-quality configuration activated for matplotlib")
        print("   • DPI: 300 (high resolution)")
        print("   • Margins: automatic (tight)")
        print("   • Format: optimized PNG")
        print("   ℹ️  Now all figures will be automatically saved in high quality")
    
    def summary_plot(self, 
                   figsize: Tuple[int, int] = DEFAULT_SUMMARY_FIGSIZE,
                   title: str = DEFAULT_PLOT_TITLE,
                   max_samples: Optional[int] = 200,
                   top_genes_count: int = DEFAULT_TOP_GENES_COUNT,
                   show_interactive: bool = False) -> plt.Figure:
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
            show_interactive: If True, display the visualization in interactive mode.
            
        Returns:
            Matplotlib figure with the summary plot.
        """
        from .visualizations.summary import create_summary_plot
        from .utils.data_processing import extract_variant_classifications, extract_variant_types
        
        # Preprocess data to ensure we have the necessary columns
        processed_data = extract_variant_classifications(
            self.data, 
            variant_column=VARIANT_CLASSIFICATION_COLUMN,
            funcotation_column=FUNCOTATION_COLUMN
        )
        
        processed_data = extract_variant_types(
            processed_data,
            variant_column=VARIANT_TYPE_COLUMN,
            funcotation_column=FUNCOTATION_COLUMN
        )
        
        # Generate the summary plot
        fig = create_summary_plot(processed_data, figsize, title, max_samples, top_genes_count)
        
        # If requested to show interactively
        if show_interactive:
            self._show_figure_interactive(fig)
        
        return fig
    
    def variant_classification_plot(self,
                                    figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                                    title: str = "Variant Classification",
                                    show_interactive: bool = False) -> plt.Figure:
        """
        Generate a horizontal bar plot showing the distribution of variant classifications.
        
        Args:
            figsize: Figure size.
            title: Plot title.
            show_interactive: If True, display the visualization in interactive mode.
            
        Returns:
            Matplotlib figure with the variant classification plot.
        """
        from .visualizations.summary import create_variant_classification_plot
        from .utils.data_processing import extract_variant_classifications
        
        # Preprocess data to ensure we have the necessary column
        processed_data = extract_variant_classifications(
            self.data, 
            variant_column="Variant_Classification",
            funcotation_column="FUNCOTATION"
        )
        
        # Create figure and axes
        fig, ax = plt.subplots(figsize=figsize)
        
        # Generate the plot, passing set_title=False to avoid duplicate title
        create_variant_classification_plot(processed_data, ax=ax, set_title=False)
        
        # Configure title
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        
        # If requested to show interactively
        if show_interactive:
            self._show_figure_interactive(fig)
        
        return fig
    
    def variant_type_plot(self,
                          figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                          title: str = "Variant Type",
                          show_interactive: bool = False) -> plt.Figure:
        """
        Generate a horizontal bar plot showing the distribution of variant types.
        
        Args:
            figsize: Figure size.
            title: Plot title.
            show_interactive: If True, display the visualization in interactive mode.
            
        Returns:
            Matplotlib figure with the variant types plot.
        """
        from .visualizations.summary import create_variant_type_plot
        from .utils.data_processing import extract_variant_types
        
        # Preprocess data to ensure we have the necessary column
        processed_data = extract_variant_types(
            self.data,
            variant_column="Variant_Type",
            funcotation_column="FUNCOTATION"
        )
        
        # Create figure and axes
        fig, ax = plt.subplots(figsize=figsize)
        
        # Generate the plot, passing set_title=False to avoid duplicate title
        create_variant_type_plot(processed_data, ax=ax, set_title=False)
        
        # Configure title
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        
        # If requested to show interactively
        if show_interactive:
            self._show_figure_interactive(fig)
        
        return fig
    
    def snv_class_plot(self,
                        figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                        title: str = "SNV Class",
                        ref_column: str = "REF",
                        alt_column: str = "ALT",
                        show_interactive: bool = False) -> plt.Figure:
        """
        Generate a horizontal bar plot showing the distribution of SNV classes.
        
        Args:
            figsize: Figure size.
            title: Plot title.
            ref_column: Name of the column containing the reference allele.
            alt_column: Name of the column containing the alternative allele.
            show_interactive: If True, display the visualization in interactive mode.
            
        Returns:
            Matplotlib figure with the SNV classes plot.
        """
        from .visualizations.summary import create_snv_class_plot
        
        fig, ax = plt.subplots(figsize=figsize)
        create_snv_class_plot(
            self.data, 
            ref_column=ref_column,
            alt_column=alt_column,
            ax=ax,
            set_title=False  # Avoid duplicate title
        )
        
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        # If requested to show interactively
        if show_interactive:
            self._show_figure_interactive(fig)
        
        return fig
    
    def variants_per_sample_plot(self,
                                 figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                                 title: str = "Variants per Sample",
                                 variant_column: str = "Variant_Classification",
                                 sample_column: str = "Tumor_Sample_Barcode",
                                 max_samples: Optional[int] = 200,
                                 show_interactive: bool = False) -> plt.Figure:
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
            show_interactive: If True, display the visualization in interactive mode.
            
        Returns:
            Matplotlib figure with the variants per sample plot.
        """
        from .visualizations.summary import create_variants_per_sample_plot
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
        processed_data = extract_variant_classifications(
            self.data, 
            variant_column=variant_column,
            funcotation_column="FUNCOTATION"
        )
        
        fig, ax = plt.subplots(figsize=figsize)
        create_variants_per_sample_plot(
            processed_data, 
            variant_column=variant_column,
            sample_column=sample_column,
            ax=ax,
            set_title=False,  # Avoid duplicate title
            max_samples=max_samples  # Pass the configured sample limit
        )
        
        # Don't modify title if it contains the median
        if title and not title.startswith("Variants per Sample"):
            fig.suptitle(title, fontsize=16, fontweight='bold', y=1.02)
        elif title:
            fig.suptitle(title, fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        
        # If requested to show interactively
        if show_interactive:
            self._show_figure_interactive(fig)
        
        return fig
    
    
    def variant_classification_summary_plot(self,
                                           figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                                           title: str = "Variant Classification Summary",
                                           variant_column: str = "Variant_Classification",
                                           sample_column: str = "Tumor_Sample_Barcode",
                                           show_interactive: bool = False) -> plt.Figure:
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
            show_interactive: If True, display the visualization in interactive mode.
            
        Returns:
            Matplotlib figure with the box-and-whiskers plot.
        """
        from .visualizations.summary import create_variant_classification_summary_plot
        from .utils.data_processing import extract_variant_classifications

        # Ensure the variant classification column exists or is extracted
        processed_data = extract_variant_classifications(
            self.data, 
            variant_column=variant_column,
            funcotation_column="FUNCOTATION"
        )
        
        # Check if we're in wide format (samples as columns)
        is_wide_format = sample_column not in processed_data.columns
        if is_wide_format:
            # Detect and show information about the format
            sample_cols = [col for col in processed_data.columns if col.startswith('TCGA-') or 
                           (isinstance(col, str) and col.count('-') >= 2)]
            if sample_cols:
                print(f"Detected wide format with {len(sample_cols)} possible sample columns.")
        
        fig, ax = plt.subplots(figsize=figsize)
        create_variant_classification_summary_plot(
            processed_data, 
            variant_column=variant_column,
            sample_column=sample_column,
            ax=ax,
            show_labels=True,  # Ensure it always shows labels when generated individually
            set_title=False  # Avoid duplicate title
        )
        
        # Configure title
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        
        # If requested to show interactively
        if show_interactive:
            self._show_figure_interactive(fig)
        
        return fig


    def top_mutated_genes_plot(self,
                              figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                              title: str = "Top Mutated Genes",
                              mode: str = MODE_VARIANTS, 
                              variant_column: str = VARIANT_CLASSIFICATION_COLUMN,
                              gene_column: str = GENE_COLUMN,
                              sample_column: str = SAMPLE_COLUMN,
                              count: int = DEFAULT_TOP_GENES_COUNT,
                              show_interactive: bool = False) -> plt.Figure:
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
            show_interactive: If True, display the visualization in interactive mode.
            
        Returns:
            Matplotlib figure with the top mutated genes plot.
        
        Raises:
            ValueError: If 'count' is not a positive number or 'mode' is not a valid value.
        """
        from .visualizations.summary import create_top_mutated_genes_plot
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
        processed_data = extract_variant_classifications(
            self.data, 
            variant_column=variant_column,
            funcotation_column=FUNCOTATION_COLUMN
        )
        
        fig, ax = plt.subplots(figsize=figsize)
        create_top_mutated_genes_plot(
            processed_data, 
            mode=mode,
            variant_column=variant_column,
            gene_column=gene_column,
            sample_column=sample_column,
            count=count,
            ax=ax,
            set_title=False  # Avoid duplicate title
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
        
        # If requested to show interactively
        if show_interactive:
            self._show_figure_interactive(fig)
        
        return fig
    
    def _show_figure_interactive(self, figure: plt.Figure) -> None:
        """
        Display a specific figure in interactive mode without affecting other figures.
        
        This private method is used internally to show only the specific figure
        in interactive mode when using show_interactive=True in visualization methods.
        
        Args:
            figure: The specific figure to display in interactive mode.
        """
        # Save current interactive mode state
        was_interactive = plt.isinteractive()
        
        # Variable to control when the window is closed
        window_closed = [False]  # Use list to make it mutable in inner function
        
        def on_close(event):
            """Callback executed when the window is closed."""
            window_closed[0] = True
        
        try:
            # Enable interactive mode temporarily only if it wasn't active
            if not was_interactive:
                plt.ion()
            
            # Connect the window close event
            figure.canvas.mpl_connect('close_event', on_close)
            
            # Show only this specific figure
            figure.show()
            
            # Force immediate figure render
            figure.canvas.draw()
            figure.canvas.flush_events()
            
            # Inform the user
            title_text = "Untitled"
            if figure._suptitle and figure._suptitle.get_text():
                title_text = figure._suptitle.get_text()
            
            print(f"Figure '{title_text}' displayed in interactive mode.")
            print("Close the window to continue with script execution.")
            
            # Wait until user closes the window
            # Use a loop that checks the window_closed variable
            import time
            while not window_closed[0] and plt.fignum_exists(figure.number):
                # Use only time.sleep() without plt.pause() to avoid side effects
                time.sleep(0.1)
                # Allow matplotlib to process events minimally
                figure.canvas.flush_events()
                    
        finally:
            # Restore original interactive mode state only if we changed it
            if not was_interactive:
                plt.ioff()  # Disable interactive mode if it wasn't active before

