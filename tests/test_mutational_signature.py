"""
Test module for mutational signature analysis functionality.
"""

import pytest
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyMut import PyMutation
from pyMut.visualizations.mutational_signature import (
    get_trinucleotide_contexts,
    get_96_category,
    extract_mutation_matrix,
    perform_nmf,
    create_signature_profile_plot,
    create_cosine_similarity_heatmap,
    create_signature_contribution_heatmap,
    create_signature_contribution_barplot,
    create_signature_donut_plot,
    create_mutational_signature_analysis
)


@pytest.fixture
def sample_mutation_data():
    """Create sample mutation data for testing."""
    np.random.seed(42)
    
    # Generate synthetic mutation data
    n_mutations = 100
    samples = [f"Sample_{i:03d}" for i in range(1, 11)]
    genes = ['TP53', 'KRAS', 'EGFR', 'BRAF', 'PIK3CA']
    nucleotides = ['A', 'C', 'G', 'T']
    
    mutations = []
    for _ in range(n_mutations):
        ref = np.random.choice(nucleotides)
        alt = np.random.choice([n for n in nucleotides if n != ref])
        upstream = np.random.choice(nucleotides)
        downstream = np.random.choice(nucleotides)
        
        mutations.append({
            'Hugo_Symbol': np.random.choice(genes),
            'Variant_Classification': 'Missense_Mutation',
            'Variant_Type': 'SNP',
            'REF': ref,
            'ALT': alt,
            'Reference_Allele': ref,
            'Tumor_Seq_Allele2': alt,
            'Tumor_Sample_Barcode': np.random.choice(samples),
            'Reference_Context': f"{upstream}{ref}{downstream}"
        })
    
    return pd.DataFrame(mutations)


@pytest.fixture
def py_mut_object(sample_mutation_data):
    """Create PyMutation object with sample data."""
    from pyMut.core import MutationMetadata
    metadata = MutationMetadata(
        source_format="test",
        file_path="test_data.tsv",
        filters=[],
        fasta="test.fasta"
    )
    return PyMutation(sample_mutation_data, metadata)


class TestMutationalSignatureHelpers:
    """Test helper functions for mutational signature analysis."""
    
    def test_get_trinucleotide_contexts(self):
        """Test trinucleotide context generation."""
        contexts = get_trinucleotide_contexts()
        assert len(contexts) == 96
        assert contexts[0] == "A[C>A]A"
        assert contexts[-1] == "T[T>G]T"
        
        # Check all substitution types are present
        for sub in ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']:
            sub_contexts = [c for c in contexts if sub in c]
            assert len(sub_contexts) == 16
    
    def test_get_96_category_valid(self):
        """Test valid category assignment."""
        # Test C>T mutation with ACG context
        assert get_96_category('C', 'T', 'ACG') == 34  # C>T, ACG context A_G
        
        # Test reverse complement conversion
        assert get_96_category('G', 'A', 'CGT') == get_96_category('C', 'T', 'ACG')
        
        # Test all substitution types
        assert get_96_category('C', 'A', 'ACA') is not None
        assert get_96_category('T', 'G', 'TTT') is not None
    
    def test_get_96_category_invalid(self):
        """Test invalid inputs."""
        assert get_96_category('C', 'T', 'A') is None  # Too short
        assert get_96_category('X', 'Y', 'ACG') is None  # Invalid nucleotides
        assert get_96_category('C', 'C', 'ACG') is None  # Same ref/alt
    
    def test_extract_mutation_matrix(self, sample_mutation_data):
        """Test mutation matrix extraction."""
        matrix, samples = extract_mutation_matrix(
            sample_mutation_data,
            context_column='Reference_Context'
        )
        
        assert matrix.shape[0] == len(samples)
        assert matrix.shape[1] == 96
        assert matrix.sum() > 0  # Should have some mutations
        assert np.all(matrix >= 0)  # No negative counts


class TestNMF:
    """Test NMF functionality."""
    
    def test_perform_nmf(self):
        """Test NMF decomposition."""
        # Create synthetic matrix
        matrix = np.random.poisson(5, size=(10, 96))
        
        W, H = perform_nmf(matrix, n_signatures=3)
        
        assert W.shape == (10, 3)
        assert H.shape == (3, 96)
        assert np.all(W >= 0)
        assert np.all(H >= 0)
        
        # Check reconstruction
        reconstruction = np.dot(W, H)
        assert reconstruction.shape == matrix.shape


class TestVisualizationFunctions:
    """Test individual visualization functions."""
    
    def test_create_signature_profile_plot(self):
        """Test signature profile plot creation."""
        # Create synthetic signatures
        signatures = np.random.random((3, 96))
        signatures = signatures / signatures.sum(axis=1, keepdims=True)
        
        fig = create_signature_profile_plot(signatures)
        
        assert isinstance(fig, plt.Figure)
        assert len(fig.axes) >= 3  # One per signature
        plt.close(fig)
    
    def test_create_cosine_similarity_heatmap(self):
        """Test cosine similarity heatmap."""
        signatures = np.random.random((3, 96))
        
        fig = create_cosine_similarity_heatmap(signatures)
        
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_create_contribution_heatmap(self):
        """Test contribution heatmap."""
        W = np.random.random((20, 3))
        samples = [f"Sample_{i}" for i in range(20)]
        
        fig = create_signature_contribution_heatmap(W, samples)
        
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_create_contribution_barplot(self):
        """Test contribution barplot."""
        W = np.random.random((20, 3))
        samples = [f"Sample_{i}" for i in range(20)]
        
        fig = create_signature_contribution_barplot(W, samples)
        
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_create_donut_plot(self):
        """Test donut plot creation."""
        W = np.random.random((20, 3))
        
        fig = create_signature_donut_plot(W)
        
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


class TestPyMutationIntegration:
    """Test PyMutation class integration."""
    
    def test_mutational_signature_plot_basic(self, py_mut_object):
        """Test basic mutational signature plot."""
        fig = py_mut_object.mutational_signature_plot(n_signatures=2)
        
        assert isinstance(fig, plt.Figure)
        assert fig.texts[0].get_text() == 'Mutational Signature Analysis'
        plt.close(fig)
    
    def test_mutational_signature_plot_custom_params(self, py_mut_object):
        """Test with custom parameters."""
        fig = py_mut_object.mutational_signature_plot(
            n_signatures=3,
            context_column='Reference_Context',
            figsize=(15, 20),
            title="Custom Title"
        )
        
        assert isinstance(fig, plt.Figure)
        assert fig.texts[0].get_text() == 'Custom Title'
        plt.close(fig)
    
    def test_mutational_signature_plot_no_context(self, sample_mutation_data):
        """Test when context column is missing."""
        # Remove context column
        data_no_context = sample_mutation_data.drop('Reference_Context', axis=1)
        from pyMut.core import MutationMetadata
        metadata = MutationMetadata(
            source_format="test",
            file_path="test_data.tsv",
            filters=[],
            fasta="test.fasta"
        )
        py_mut = PyMutation(data_no_context, metadata)
        
        # Should still work (might generate synthetic contexts or handle gracefully)
        fig = py_mut.mutational_signature_plot(n_signatures=2)
        
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_mutational_signature_plot_empty_data(self):
        """Test with empty data."""
        empty_data = pd.DataFrame({
            'Hugo_Symbol': [],
            'Variant_Classification': [],
            'REF': [],
            'ALT': [],
            'Tumor_Sample_Barcode': []
        })
        from pyMut.core import MutationMetadata
        metadata = MutationMetadata(
            source_format="test",
            file_path="test_data.tsv",
            filters=[],
            fasta="test.fasta"
        )
        py_mut = PyMutation(empty_data, metadata)
        
        fig = py_mut.mutational_signature_plot(n_signatures=3)
        
        assert isinstance(fig, plt.Figure)
        # Should show "No mutations found" message
        plt.close(fig)
    
    def test_mutational_signature_plot_with_cosmic(self, py_mut_object):
        """Test with COSMIC signatures."""
        # Create mock COSMIC signatures
        cosmic_sigs = np.random.random((96, 10))
        cosmic_df = pd.DataFrame(
            cosmic_sigs,
            columns=[f'SBS{i+1}' for i in range(10)]
        )
        
        fig = py_mut_object.mutational_signature_plot(
            n_signatures=3,
            cosmic_signatures=cosmic_df
        )
        
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_single_sample(self):
        """Test with single sample."""
        data = pd.DataFrame({
            'Hugo_Symbol': ['TP53'] * 10,
            'Variant_Classification': ['Missense_Mutation'] * 10,
            'REF': ['C'] * 10,
            'ALT': ['T'] * 10,
            'Tumor_Sample_Barcode': ['Sample_001'] * 10,
            'Reference_Context': ['ACG'] * 10
        })
        from pyMut.core import MutationMetadata
        metadata = MutationMetadata(
            source_format="test",
            file_path="test_data.tsv",
            filters=[],
            fasta="test.fasta"
        )
        py_mut = PyMutation(data, metadata)
        
        fig = py_mut.mutational_signature_plot(n_signatures=1)
        
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_invalid_nucleotides(self):
        """Test with invalid nucleotides."""
        data = pd.DataFrame({
            'Hugo_Symbol': ['TP53'],
            'Variant_Classification': ['Missense_Mutation'],
            'REF': ['N'],  # Invalid
            'ALT': ['X'],  # Invalid
            'Tumor_Sample_Barcode': ['Sample_001'],
            'Reference_Context': ['NNN']
        })
        from pyMut.core import MutationMetadata
        metadata = MutationMetadata(
            source_format="test",
            file_path="test_data.tsv",
            filters=[],
            fasta="test.fasta"
        )
        py_mut = PyMutation(data, metadata)
        
        # Should handle gracefully
        fig = py_mut.mutational_signature_plot(n_signatures=2)
        
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_large_signature_number(self, py_mut_object):
        """Test with large number of signatures."""
        # Request more signatures than samples
        fig = py_mut_object.mutational_signature_plot(n_signatures=20)
        
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_mixed_variant_types(self):
        """Test with mixed variant types (not just SNPs)."""
        data = pd.DataFrame({
            'Hugo_Symbol': ['TP53'] * 4,
            'Variant_Classification': ['Missense_Mutation'] * 4,
            'Variant_Type': ['SNP', 'INS', 'DEL', 'SNP'],
            'REF': ['C', 'A', 'ATG', 'G'],
            'ALT': ['T', 'AGG', 'A', 'A'],
            'Tumor_Sample_Barcode': ['Sample_001'] * 4,
            'Reference_Context': ['ACG', 'TAT', 'ATGC', 'CGT']
        })
        from pyMut.core import MutationMetadata
        metadata = MutationMetadata(
            source_format="test",
            file_path="test_data.tsv",
            filters=[],
            fasta="test.fasta"
        )
        py_mut = PyMutation(data, metadata)
        
        # Should filter to only SNPs
        fig = py_mut.mutational_signature_plot(n_signatures=2)
        
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


# Run tests if executed directly
if __name__ == "__main__":
    pytest.main([__file__, "-v"]) 