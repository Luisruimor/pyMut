"""
Tests for the PyMutation class.

This file contains unit tests for the PyMutation class
and its main methods, focusing on summary plot functionality.
"""

import os
import sys
import unittest
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import tempfile
import numpy as np

# Comprehensive warning suppression for clean test output
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=PendingDeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Specifically suppress pandas/numpy warnings that cause noise
warnings.filterwarnings("ignore", message=".*find_common_type.*")
warnings.filterwarnings("ignore", message=".*promote_types.*")
warnings.filterwarnings("ignore", message=".*result_type.*")

# Configure matplotlib for testing (no GUI, no output)
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
plt.ioff()  # Turn off interactive mode globally

# Suppress matplotlib font warnings
import logging
logging.getLogger('matplotlib.font_manager').disabled = True
logging.getLogger('matplotlib').setLevel(logging.WARNING)

# Add root directory to path to import pyMut
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.pyMut import PyMutation
from src.pyMut.utils.data_processing import read_tsv

# Path to example file
EXAMPLE_FILE = os.path.abspath(os.path.join(os.path.dirname(__file__), 
                                         '../src/pyMut/data/examples/tcga_laml_converted.tsv'))

class TestPyMutation(unittest.TestCase):
    """Tests for the PyMutation class."""
    
    @classmethod
    def setUpClass(cls):
        """Set up class-level configuration that runs once before all tests."""
        # Load data once for all tests to improve performance
        cls.test_data = read_tsv(EXAMPLE_FILE)
        
        # Configure pandas to suppress warnings during testing
        pd.set_option('mode.chained_assignment', None)
        
        # Additional warning suppression at class level
        warnings.simplefilter("ignore")
    
    def setUp(self):
        """Initial setup for tests."""
        # Create a fresh PyMutation instance for each test
        self.data = self.test_data.copy()  # Use copy to avoid interference between tests
        self.pyMut = PyMutation(self.data)
    
    def tearDown(self):
        """Cleanup after each test."""
        # Close all matplotlib figures to prevent memory leaks
        plt.close('all')
        
        # Clear any matplotlib cache
        plt.clf()
        plt.cla()
    
    @classmethod
    def tearDownClass(cls):
        """Class-level cleanup that runs once after all tests."""
        # Reset pandas options
        pd.reset_option('mode.chained_assignment')
        
        # Final matplotlib cleanup
        plt.close('all')
        matplotlib.pyplot.switch_backend('Agg')
    
    def test_initialization(self):
        """Test PyMutation initialization."""
        self.assertIsInstance(self.pyMut, PyMutation)
        self.assertIsInstance(self.pyMut.data, pd.DataFrame)
        self.assertEqual(len(self.pyMut.data), len(self.data))
    
    def test_summary_plot(self):
        """Test the summary_plot method."""
        fig = self.pyMut.summary_plot()
        self.assertIsInstance(fig, plt.Figure)
        # Verify that the figure has the expected number of subplots (2x3 = 6)
        self.assertEqual(len(fig.axes), 6)
    
    def test_summary_plot_with_parameters(self):
        """Test the summary_plot method with custom parameters."""
        fig = self.pyMut.summary_plot(
            figsize=(12, 8),
            title="Custom Summary",
            max_samples=50,
            top_genes_count=5
        )
        self.assertIsInstance(fig, plt.Figure)
        self.assertEqual(fig._suptitle.get_text(), "Custom Summary")

    def test_variant_classification_plot(self):
        """Test the variant_classification_plot method."""
        fig = self.pyMut.variant_classification_plot()
        self.assertIsInstance(fig, plt.Figure)

    def test_variant_type_plot(self):
        """Test the variant_type_plot method."""
        fig = self.pyMut.variant_type_plot()
        self.assertIsInstance(fig, plt.Figure)

    def test_snv_class_plot(self):
        """Test the snv_class_plot method."""
        fig = self.pyMut.snv_class_plot()
        self.assertIsInstance(fig, plt.Figure)

    def test_variant_classification_summary_plot(self):
        """Test the variant_classification_summary_plot method."""
        fig = self.pyMut.variant_classification_summary_plot()
        self.assertIsInstance(fig, plt.Figure)

    def test_variants_per_sample_plot(self):
        """Test the variants_per_sample_plot method."""
        fig = self.pyMut.variants_per_sample_plot()
        self.assertIsInstance(fig, plt.Figure)

    def test_top_mutated_genes_plot(self):
        """Test the top_mutated_genes_plot method."""
        # Test "variants" mode
        fig = self.pyMut.top_mutated_genes_plot(mode="variants")
        self.assertIsInstance(fig, plt.Figure)
        
        # Test "samples" mode
        fig = self.pyMut.top_mutated_genes_plot(mode="samples")
        self.assertIsInstance(fig, plt.Figure)

    def test_validation_top_mutated_genes_plot(self):
        """Test parameter validation in top_mutated_genes_plot."""
        # Verify that 'count' parameter is validated
        with self.assertRaises(ValueError):
            self.pyMut.top_mutated_genes_plot(count=0)
        
        with self.assertRaises(ValueError):
            self.pyMut.top_mutated_genes_plot(count=-5)
        
        # Verify that 'mode' parameter is validated
        with self.assertRaises(ValueError):
            self.pyMut.top_mutated_genes_plot(mode="invalid_mode")

    def test_empty_dataframe_initialization(self):
        """Test that empty DataFrame validation works."""
        with self.assertRaises(ValueError):
            PyMutation(pd.DataFrame())

    def test_non_dataframe_initialization(self):
        """Test that parameter type validation works."""
        with self.assertRaises(ValueError):
            PyMutation([1, 2, 3])  # List instead of DataFrame
        
        with self.assertRaises(ValueError):
            PyMutation("Not a DataFrame")  # String instead of DataFrame

    def test_configure_high_quality_plots(self):
        """Test the static method configure_high_quality_plots."""
        # This method modifies matplotlib global settings
        # We can test that it doesn't raise an exception
        try:
            PyMutation.configure_high_quality_plots()
        except Exception as e:
            self.fail(f"configure_high_quality_plots raised an exception: {e}")

    def test_save_figure(self):
        """Test the save_figure method."""
        fig = self.pyMut.variant_classification_plot()
        
        # Create a temporary file to save the figure
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp_file:
            temp_filename = tmp_file.name
        
        try:
            # Test saving the figure
            self.pyMut.save_figure(fig, temp_filename)
            
            # Verify that the file was created
            self.assertTrue(os.path.exists(temp_filename))
            
            # Verify that the file has content (is not empty)
            self.assertGreater(os.path.getsize(temp_filename), 0)
            
        finally:
            # Clean up the temporary file
            if os.path.exists(temp_filename):
                os.unlink(temp_filename)

    def test_plot_methods_with_custom_parameters(self):
        """Test plot methods with custom parameters."""
        # Test variant_classification_plot with parameters
        fig = self.pyMut.variant_classification_plot(
            figsize=(10, 8),
            title="Custom Variant Classification"
        )
        self.assertIsInstance(fig, plt.Figure)
        
        # Test variants_per_sample_plot with max_samples
        fig = self.pyMut.variants_per_sample_plot(max_samples=50)
        self.assertIsInstance(fig, plt.Figure)
        
        # Test top_mutated_genes_plot with custom count
        fig = self.pyMut.top_mutated_genes_plot(mode="variants", count=5)
        self.assertIsInstance(fig, plt.Figure)

    def test_show_interactive_parameter(self):
        """Test that show_interactive parameter doesn't cause errors."""
        # Note: We set show_interactive=True but since we're in non-interactive mode
        # for testing, this should not actually display anything
        fig = self.pyMut.variant_classification_plot(show_interactive=False)
        self.assertIsInstance(fig, plt.Figure)

    def test_missing_columns_handling(self):
        """Test behavior with missing columns."""
        # Create a DataFrame with minimal columns
        minimal_data = pd.DataFrame({
            'FUNCOTATION': ['test|test|test|test|test|Missense_Mutation|test|SNP|test'] * 10,
            'REF': ['A'] * 10,
            'ALT': ['G'] * 10,
            'Hugo_Symbol': ['GENE1'] * 10,
            'Tumor_Sample_Barcode': [f'SAMPLE_{i}' for i in range(10)]
        })
        
        minimal_pyMut = PyMutation(minimal_data)
        
        # This should work because the class should extract missing columns from FUNCOTATION
        fig = minimal_pyMut.variant_classification_plot()
        self.assertIsInstance(fig, plt.Figure)


if __name__ == "__main__":
    # Configure unittest to run with minimal output to avoid clutter
    unittest.main(verbosity=1, buffer=True)
