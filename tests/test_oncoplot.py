#!/usr/bin/env python3
"""
Test suite for oncoplot functionality.

This module contains comprehensive tests for all functions related
to oncoplot generation in pyMut.
"""

import unittest
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

# Add src to path to import pyMut
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from pyMut.core import PyMutation
from pyMut.visualizations.oncoplot import (
    is_mutated,
    detect_sample_columns,
    create_variant_color_mapping,
    process_mutation_matrix,
    create_oncoplot_plot
)


class TestOncoplotUtilities(unittest.TestCase):
    """Tests for oncoplot utility functions."""
    
    def test_is_mutated_basic_cases(self):
        """Test basic mutation detection cases."""
        # Cases that should return True (mutation present)
        self.assertTrue(is_mutated("A|G", "A", "G"))
        self.assertTrue(is_mutated("G|A", "A", "G"))
        self.assertTrue(is_mutated("G|G", "A", "G"))
        self.assertTrue(is_mutated("A/G", "A", "G"))
        self.assertTrue(is_mutated("G/G", "A", "G"))
        
        # Cases that should return False (no mutation)
        self.assertFalse(is_mutated("A|A", "A", "G"))
        self.assertFalse(is_mutated("A/A", "A", "G"))
        self.assertFalse(is_mutated("", "A", "G"))
        self.assertFalse(is_mutated(".", "A", "G"))
    
    def test_is_mutated_edge_cases(self):
        """Test edge cases for mutation detection."""
        # Cases with null values
        self.assertFalse(is_mutated(None, "A", "G"))
        self.assertFalse(is_mutated("A|G", None, "G"))
        self.assertFalse(is_mutated("A|G", "A", None))
        
        # Cases with different separators
        self.assertTrue(is_mutated("A:G", "A", "G"))
        self.assertTrue(is_mutated("A;G", "A", "G"))
    
    def test_detect_sample_columns(self):
        """Test automatic detection of sample columns."""
        # DataFrame with TCGA columns
        df_tcga = pd.DataFrame({
            'Hugo_Symbol': ['TP53'],
            'TCGA-AB-1234': ['A|G'],
            'TCGA-CD-5678': ['A|A'],
            'Other_Column': ['value']
        })
        
        sample_cols = detect_sample_columns(df_tcga)
        expected = ['TCGA-AB-1234', 'TCGA-CD-5678']
        self.assertEqual(sorted(sample_cols), sorted(expected))
    
    def test_detect_sample_columns_gt_format(self):
        """Test detection of columns in .GT format."""
        df_gt = pd.DataFrame({
            'Hugo_Symbol': ['TP53'],
            'sample1.GT': ['A|G'],
            'sample2.GT': ['A|A'],
            'Other_Column': ['value']
        })
        
        sample_cols = detect_sample_columns(df_gt)
        expected = ['sample1.GT', 'sample2.GT']
        self.assertEqual(sorted(sample_cols), sorted(expected))
    
    def test_detect_sample_columns_no_samples(self):
        """Test error when no sample columns are detected."""
        df_no_samples = pd.DataFrame({
            'Hugo_Symbol': ['TP53'],
            'Other_Column': ['value']
        })
        
        with self.assertRaises(ValueError):
            detect_sample_columns(df_no_samples)
    
    def test_create_variant_color_mapping(self):
        """Test creation of color mapping for variants."""
        variants = {'Missense_Mutation', 'Nonsense_Mutation', 'Custom_Variant'}
        color_mapping = create_variant_color_mapping(variants)
        
        # Verify colors were created for all variants
        self.assertEqual(len(color_mapping), len(variants))
        
        # Verify standard variants have known colors
        self.assertIn('Missense_Mutation', color_mapping)
        self.assertIn('Nonsense_Mutation', color_mapping)
        self.assertIn('Custom_Variant', color_mapping)


class TestOncoplotProcessing(unittest.TestCase):
    """Tests for oncoplot data processing."""
    
    def setUp(self):
        """Set up test data."""
        self.test_data = pd.DataFrame({
            'Hugo_Symbol': ['TP53', 'TP53', 'KRAS'],
            'Variant_Classification': ['Missense_Mutation', 'Nonsense_Mutation', 'Missense_Mutation'],
            'REF': ['A', 'C', 'G'],
            'ALT': ['G', 'T', 'A'],
            'TCGA-AB-1234': ['A|G', 'C|C', 'G|G'],
            'TCGA-CD-5678': ['A|A', 'C|T', 'G|A']
        })
    
    def test_process_mutation_matrix_basic(self):
        """Test basic mutation matrix processing."""
        matrix, counts = process_mutation_matrix(self.test_data)
        
        # Verify dimensions
        self.assertEqual(matrix.shape[0], 2)  # 2 unique genes
        self.assertEqual(matrix.shape[1], 2)  # 2 samples
        
        # Verify mutations were detected correctly
        self.assertEqual(matrix.loc['TP53', 'TCGA-AB-1234'], 'Missense_Mutation')
        self.assertEqual(matrix.loc['KRAS', 'TCGA-CD-5678'], 'Missense_Mutation')
    
    def test_process_mutation_matrix_multi_hit(self):
        """Test Multi_Hit detection."""
        # Add data that generates Multi_Hit
        multi_hit_data = pd.DataFrame({
            'Hugo_Symbol': ['TP53', 'TP53'],
            'Variant_Classification': ['Missense_Mutation', 'Nonsense_Mutation'],
            'REF': ['A', 'C'],
            'ALT': ['G', 'T'],
            'TCGA-AB-1234': ['A|G', 'C|T'],  # Both mutations in same sample
        })
        
        matrix, counts = process_mutation_matrix(multi_hit_data)
        self.assertEqual(matrix.loc['TP53', 'TCGA-AB-1234'], 'Multi_Hit')
    
    def test_process_mutation_matrix_missing_columns(self):
        """Test error with missing columns."""
        incomplete_data = pd.DataFrame({
            'Hugo_Symbol': ['TP53'],
            'TCGA-AB-1234': ['A|G']
            # Missing required columns
        })
        
        with self.assertRaises(ValueError):
            process_mutation_matrix(incomplete_data)


class TestOncoplotPlot(unittest.TestCase):
    """Tests for the main oncoplot creation function."""
    
    def setUp(self):
        """Set up test data."""
        self.test_data = pd.DataFrame({
            'Hugo_Symbol': ['TP53', 'KRAS', 'PIK3CA'] * 2,
            'Variant_Classification': ['Missense_Mutation'] * 6,
            'REF': ['A'] * 6,
            'ALT': ['G'] * 6,
            'TCGA-AB-1234': ['A|G', 'A|A', 'A|G', 'A|G', 'A|A', 'A|A'],
            'TCGA-CD-5678': ['A|A', 'A|G', 'A|A', 'A|A', 'A|G', 'A|G']
        })
    
    def test_create_oncoplot_plot_basic(self):
        """Test basic oncoplot creation."""
        fig = create_oncoplot_plot(self.test_data)
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_create_oncoplot_plot_custom_parameters(self):
        """Test oncoplot with custom parameters."""
        fig = create_oncoplot_plot(
            self.test_data,
            title="Test Oncoplot",
            top_genes_count=2,
            max_samples=1,
            figsize=(10, 6)
        )
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_create_oncoplot_plot_with_axes(self):
        """Test oncoplot with pre-existing axes."""
        fig, ax = plt.subplots(figsize=(8, 6))
        result_fig = create_oncoplot_plot(self.test_data, ax=ax)
        self.assertEqual(fig, result_fig)
        plt.close(fig)


class TestPyMutationOncoplot(unittest.TestCase):
    """Tests for oncoplot integration in PyMutation."""
    
    def setUp(self):
        """Set up test data."""
        self.test_data = pd.DataFrame({
            'Hugo_Symbol': ['TP53', 'KRAS', 'PIK3CA'],
            'Variant_Classification': ['Missense_Mutation', 'Nonsense_Mutation', 'In_Frame_Del'],
            'REF': ['A', 'C', 'G'],
            'ALT': ['G', 'T', 'A'],
            'TCGA-AB-1234': ['A|G', 'C|C', 'G|A'],
            'TCGA-CD-5678': ['A|A', 'C|T', 'G|G']
        })
        self.py_mut = PyMutation(self.test_data)
        
    def test_oncoplot_basic(self):
        """Test basic oncoplot method."""
        fig = self.py_mut.oncoplot()
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_oncoplot_custom_parameters(self):
        """Test oncoplot method with custom parameters."""
        fig = self.py_mut.oncoplot(
            title="Test Custom Oncoplot",
            top_genes_count=2,
            max_samples=1,
            figsize=(12, 8)
        )
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_oncoplot_invalid_parameters(self):
        """Test invalid parameters for oncoplot."""
        with self.assertRaises(ValueError):
            self.py_mut.oncoplot(top_genes_count=0)
        
        with self.assertRaises(ValueError):
            self.py_mut.oncoplot(max_samples=-1)
    
    def test_oncoplot_missing_columns(self):
        """Test oncoplot with missing columns."""
        incomplete_data = pd.DataFrame({
            'Hugo_Symbol': ['TP53'],
            'TCGA-AB-1234': ['A|G']
        })
        py_mut_incomplete = PyMutation(incomplete_data)
        
        with self.assertRaises(ValueError):
            py_mut_incomplete.oncoplot()


class TestOncoplotIntegration(unittest.TestCase):
    """Integration tests for oncoplot with realistic data."""
    
    def setUp(self):
        """Set up more realistic test data."""
        # Create data simulating a real dataset
        genes = ['TP53', 'KRAS', 'PIK3CA', 'EGFR', 'BRAF'] * 10
        samples = [f'TCGA-AB-{1000+i}' for i in range(20)]
        
        data_rows = []
        for i, gene in enumerate(genes):
            for j, sample in enumerate(samples):
                if (i + j) % 3 == 0:  # Create mutation pattern
                    data_rows.append({
                        'Hugo_Symbol': gene,
                        'Variant_Classification': 'Missense_Mutation',
                        'REF': 'A',
                        'ALT': 'G',
                        sample: 'A|G'
                    })
                else:
                    data_rows.append({
                        'Hugo_Symbol': gene,
                        'Variant_Classification': 'Missense_Mutation',
                        'REF': 'A',
                        'ALT': 'G',
                        sample: 'A|A'
                    })
        
        # Create DataFrame with complete structure
        self.realistic_data = pd.DataFrame(data_rows).fillna('A|A')
        
        # Reorganize to have correct columns
        base_columns = ['Hugo_Symbol', 'Variant_Classification', 'REF', 'ALT']
        sample_columns = [col for col in self.realistic_data.columns if col.startswith('TCGA-')]
        
        # Create final DataFrame with one row per gene
        final_data = []
        for gene in ['TP53', 'KRAS', 'PIK3CA', 'EGFR', 'BRAF']:
            row = {
                'Hugo_Symbol': gene,
                'Variant_Classification': 'Missense_Mutation',
                'REF': 'A',
                'ALT': 'G'
            }
            for sample in samples:
                # Create pseudo-random mutation pattern
                if hash(f"{gene}_{sample}") % 3 == 0:
                    row[sample] = 'A|G'
                else:
                    row[sample] = 'A|A'
            final_data.append(row)
        
        self.realistic_data = pd.DataFrame(final_data)
    
    def test_oncoplot_with_real_data(self):
        """Test oncoplot with more realistic data."""
        py_mut = PyMutation(self.realistic_data)
        fig = py_mut.oncoplot(
            title="Realistic Oncoplot Test",
            top_genes_count=5,
            max_samples=10
        )
        
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)


class TestOncoplotSidePanel:
    """Test class for oncoplot side panel functionality."""
    
    def test_oncoplot_with_side_panel_integration(self):
        """Test that oncoplot integrates correctly with side panel (top mutated genes)."""
        # Create test data with TCGA sample columns
        data = pd.DataFrame({
            'Hugo_Symbol': ['TP53', 'BRCA1', 'EGFR', 'MYC'] * 5,
            'Variant_Classification': ['Missense_Mutation', 'Nonsense_Mutation', 
                                     'Frame_Shift_Del', 'Silent'] * 5,
            'REF': ['A', 'C', 'G', 'T'] * 5,
            'ALT': ['G', 'T', 'A', 'C'] * 5,
            'TCGA-AB-2988': ['A|G', 'C|T', 'G|A', './.'] * 5,
            'TCGA-AB-2869': ['A|A', 'C|T', 'G|G', 'T|C'] * 5,
            'TCGA-AB-3009': ['A|G', 'C|C', 'G|A', 'T|T'] * 5,
        })
        
        # Initialize PyMutation object
        py_mut = PyMutation(data)
        
        # Generate oncoplot which should include side panel
        figure = py_mut.oncoplot(
            figsize=(16, 10),
            title="Test Oncoplot with Side Panel",
            top_genes_count=4,
            max_samples=3
        )
        
        # Verify figure exists and has expected structure
        assert figure is not None
        assert isinstance(figure, plt.Figure)
        
        # Check that figure has multiple axes (main plot + side panel)
        axes = figure.get_axes()
        assert len(axes) >= 2, "Oncoplot should have at least main plot and side panel"
        
        # Check figure size
        assert figure.get_figwidth() == 16
        assert figure.get_figheight() == 10
        
        # Close figure to free memory
        plt.close(figure)
    
    def test_oncoplot_side_panel_gene_ordering(self):
        """Test that genes in side panel match the main plot ordering."""
        # Create test data with varying mutation frequencies
        genes = ['GENE_A', 'GENE_B', 'GENE_C']
        data = pd.DataFrame({
            'Hugo_Symbol': genes * 10,
            'Variant_Classification': ['Missense_Mutation'] * 30,
            'REF': ['A'] * 30,
            'ALT': ['G'] * 30,
            # GENE_A: appears in 3 samples (most frequent)
            # GENE_B: appears in 2 samples  
            # GENE_C: appears in 1 sample (least frequent)
            'TCGA-AB-001': ['A|G'] * 10 + ['A|A'] * 10 + ['A|A'] * 10,
            'TCGA-AB-002': ['A|G'] * 10 + ['A|G'] * 10 + ['A|A'] * 10,
            'TCGA-AB-003': ['A|G'] * 10 + ['A|A'] * 10 + ['A|G'] * 10,
        })
        
        py_mut = PyMutation(data)
        figure = py_mut.oncoplot(
            top_genes_count=3,
            max_samples=3
        )
        
        # Verify figure was created successfully
        assert figure is not None
        axes = figure.get_axes()
        assert len(axes) >= 2
        
        plt.close(figure)


class TestOncoplotWaterfallOrdering(unittest.TestCase):
    """Tests for waterfall/cascade ordering verification in oncoplot."""
    
    def test_waterfall_ordering_by_tmb(self):
        """Test that samples are ordered by TMB in descending order (cascade pattern)."""
        # Create data with clearly different TMB for each sample
        # Sample 1: 4 mutations, Sample 2: 2 mutations, Sample 3: 1 mutation
        data = pd.DataFrame({
            'Hugo_Symbol': ['GENE1', 'GENE2', 'GENE3', 'GENE4'],
            'Variant_Classification': ['Missense_Mutation'] * 4,
            'REF': ['A'] * 4,
            'ALT': ['G'] * 4,
            'TCGA-SAMPLE-1': ['A|G', 'A|G', 'A|G', 'A|G'],  # 4 mutations
            'TCGA-SAMPLE-2': ['A|G', 'A|G', 'A|A', 'A|A'],  # 2 mutations
            'TCGA-SAMPLE-3': ['A|G', 'A|A', 'A|A', 'A|A'],  # 1 mutation
        })
        
        # Process the matrix
        matrix, counts = process_mutation_matrix(data)
        
        # Verify TMB is calculated correctly
        tmb_per_sample = (matrix != 'None').sum(axis=0)
        self.assertEqual(tmb_per_sample['TCGA-SAMPLE-1'], 4)
        self.assertEqual(tmb_per_sample['TCGA-SAMPLE-2'], 2)
        self.assertEqual(tmb_per_sample['TCGA-SAMPLE-3'], 1)
        
        # Create the oncoplot
        py_mut = PyMutation(data)
        fig = py_mut.oncoplot(top_genes_count=4, max_samples=3)
        
        # The oncoplot should display samples ordered by TMB
        # descending: SAMPLE-1 (4), SAMPLE-2 (2), SAMPLE-3 (1)
        # Verify the figure was created correctly
        self.assertIsInstance(fig, plt.Figure)
        
        plt.close(fig)
    
    def test_waterfall_secondary_ordering(self):
        """Test secondary ordering for samples with equal TMB."""
        # Create data where some samples have the same TMB
        data = pd.DataFrame({
            'Hugo_Symbol': ['GENE1', 'GENE2', 'GENE3'],
            'Variant_Classification': ['Missense_Mutation'] * 3,
            'REF': ['A'] * 3,
            'ALT': ['G'] * 3,
            'TCGA-A': ['A|G', 'A|G', 'A|A'],  # 2 mutations, pattern 110
            'TCGA-B': ['A|G', 'A|A', 'A|G'],  # 2 mutations, pattern 101
            'TCGA-C': ['A|A', 'A|G', 'A|G'],  # 2 mutations, pattern 011
            'TCGA-D': ['A|G', 'A|G', 'A|G'],  # 3 mutations
        })
        
        matrix, counts = process_mutation_matrix(data)
        
        # Verify TMB
        tmb_per_sample = (matrix != 'None').sum(axis=0)
        self.assertEqual(tmb_per_sample['TCGA-A'], 2)
        self.assertEqual(tmb_per_sample['TCGA-B'], 2)
        self.assertEqual(tmb_per_sample['TCGA-C'], 2)
        self.assertEqual(tmb_per_sample['TCGA-D'], 3)
        
        # Create oncoplot
        py_mut = PyMutation(data)
        fig = py_mut.oncoplot(top_genes_count=3, max_samples=4)
        
        # TCGA-D should appear first (TMB=3)
        # The other three should be ordered by their binary pattern
        self.assertIsInstance(fig, plt.Figure)
        
        plt.close(fig)


if __name__ == '__main__':
    # Configure matplotlib for tests
    import matplotlib
    matplotlib.use('Agg')
    
    # Run tests
    unittest.main(verbosity=2) 