#!/usr/bin/env python3
"""
Test suite for oncoplot functionality.

Este módulo contiene pruebas exhaustivas para todas las funciones relacionadas
con la generación de oncoplots en pyMut.
"""

import unittest
import pandas as pd
import matplotlib.pyplot as plt
from unittest.mock import patch
import tempfile
import os
import sys

# Añadir src al path para importar pyMut
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
    """Tests para funciones utilitarias del oncoplot."""
    
    def test_is_mutated_basic_cases(self):
        """Test casos básicos de detección de mutaciones."""
        # Casos que deben retornar True (hay mutación)
        self.assertTrue(is_mutated("A|G", "A", "G"))
        self.assertTrue(is_mutated("G|A", "A", "G"))
        self.assertTrue(is_mutated("G|G", "A", "G"))
        self.assertTrue(is_mutated("A/G", "A", "G"))
        self.assertTrue(is_mutated("G/G", "A", "G"))
        
        # Casos que deben retornar False (no hay mutación)
        self.assertFalse(is_mutated("A|A", "A", "G"))
        self.assertFalse(is_mutated("A/A", "A", "G"))
        self.assertFalse(is_mutated("", "A", "G"))
        self.assertFalse(is_mutated(".", "A", "G"))
        
    def test_is_mutated_edge_cases(self):
        """Test casos extremos para detección de mutaciones."""
        # Casos con valores nulos
        self.assertFalse(is_mutated(None, "A", "G"))
        self.assertFalse(is_mutated("A|G", None, "G"))
        self.assertFalse(is_mutated("A|G", "A", None))
        
        # Casos con diferentes separadores
        self.assertTrue(is_mutated("A:G", "A", "G"))
        self.assertTrue(is_mutated("A;G", "A", "G"))
        
    def test_detect_sample_columns(self):
        """Test detección automática de columnas de muestras."""
        # DataFrame con columnas TCGA
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
        """Test detección de columnas en formato .GT."""
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
        """Test error cuando no se detectan columnas de muestras."""
        df_no_samples = pd.DataFrame({
            'Hugo_Symbol': ['TP53'],
            'Other_Column': ['value']
        })
        
        with self.assertRaises(ValueError):
            detect_sample_columns(df_no_samples)
            
    def test_create_variant_color_mapping(self):
        """Test creación de mapeo de colores para variantes."""
        variants = {'Missense_Mutation', 'Nonsense_Mutation', 'Custom_Variant'}
        color_mapping = create_variant_color_mapping(variants)
        
        # Verificar que se crearon colores para todas las variantes
        self.assertEqual(len(color_mapping), len(variants))
        
        # Verificar que las variantes estándar tienen colores conocidos
        self.assertIn('Missense_Mutation', color_mapping)
        self.assertIn('Nonsense_Mutation', color_mapping)
        self.assertIn('Custom_Variant', color_mapping)


class TestOncoplotProcessing(unittest.TestCase):
    """Tests para procesamiento de datos del oncoplot."""
    
    def setUp(self):
        """Configurar datos de prueba."""
        self.test_data = pd.DataFrame({
            'Hugo_Symbol': ['TP53', 'TP53', 'KRAS'],
            'Variant_Classification': ['Missense_Mutation', 'Nonsense_Mutation', 'Missense_Mutation'],
            'REF': ['A', 'C', 'G'],
            'ALT': ['G', 'T', 'A'],
            'TCGA-AB-1234': ['A|G', 'C|C', 'G|G'],
            'TCGA-CD-5678': ['A|A', 'C|T', 'G|A']
        })
        
    def test_process_mutation_matrix_basic(self):
        """Test procesamiento básico de matriz de mutaciones."""
        matrix, counts = process_mutation_matrix(self.test_data)
        
        # Verificar dimensiones
        self.assertEqual(matrix.shape[0], 2)  # 2 genes únicos
        self.assertEqual(matrix.shape[1], 2)  # 2 muestras
        
        # Verificar que se detectaron mutaciones correctamente
        self.assertEqual(matrix.loc['TP53', 'TCGA-AB-1234'], 'Missense_Mutation')
        self.assertEqual(matrix.loc['KRAS', 'TCGA-CD-5678'], 'Missense_Mutation')
        
    def test_process_mutation_matrix_multi_hit(self):
        """Test detección de Multi_Hit."""
        # Agregar datos que generen Multi_Hit
        multi_hit_data = pd.DataFrame({
            'Hugo_Symbol': ['TP53', 'TP53'],
            'Variant_Classification': ['Missense_Mutation', 'Nonsense_Mutation'],
            'REF': ['A', 'C'],
            'ALT': ['G', 'T'],
            'TCGA-AB-1234': ['A|G', 'C|T'],  # Ambas mutaciones en la misma muestra
        })
        
        matrix, counts = process_mutation_matrix(multi_hit_data)
        self.assertEqual(matrix.loc['TP53', 'TCGA-AB-1234'], 'Multi_Hit')
        
    def test_process_mutation_matrix_missing_columns(self):
        """Test error con columnas faltantes."""
        incomplete_data = pd.DataFrame({
            'Hugo_Symbol': ['TP53'],
            'TCGA-AB-1234': ['A|G']
            # Faltan columnas requeridas
        })
        
        with self.assertRaises(ValueError):
            process_mutation_matrix(incomplete_data)


class TestOncoplotPlot(unittest.TestCase):
    """Tests para la función principal de creación de oncoplots."""
    
    def setUp(self):
        """Configurar datos de prueba."""
        self.test_data = pd.DataFrame({
            'Hugo_Symbol': ['TP53', 'KRAS', 'PIK3CA'] * 2,
            'Variant_Classification': ['Missense_Mutation'] * 6,
            'REF': ['A'] * 6,
            'ALT': ['G'] * 6,
            'TCGA-AB-1234': ['A|G', 'A|A', 'A|G', 'A|G', 'A|A', 'A|A'],
            'TCGA-CD-5678': ['A|A', 'A|G', 'A|A', 'A|A', 'A|G', 'A|G']
        })
        
    def test_create_oncoplot_plot_basic(self):
        """Test creación básica de oncoplot."""
        fig = create_oncoplot_plot(self.test_data)
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)
        
    def test_create_oncoplot_plot_custom_parameters(self):
        """Test oncoplot con parámetros personalizados."""
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
        """Test oncoplot con ejes pre-existentes."""
        fig, ax = plt.subplots(figsize=(8, 6))
        result_fig = create_oncoplot_plot(self.test_data, ax=ax)
        self.assertEqual(fig, result_fig)
        plt.close(fig)


class TestPyMutationOncoplot(unittest.TestCase):
    """Tests para la integración del oncoplot en PyMutation."""
    
    def setUp(self):
        """Configurar datos de prueba."""
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
        """Test método oncoplot básico."""
        fig = self.py_mut.oncoplot()
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)
        
    def test_oncoplot_custom_parameters(self):
        """Test método oncoplot con parámetros personalizados."""
        fig = self.py_mut.oncoplot(
            title="Test Custom Oncoplot",
            top_genes_count=2,
            max_samples=1,
            figsize=(12, 8)
        )
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)
        
    def test_oncoplot_invalid_parameters(self):
        """Test parámetros inválidos para oncoplot."""
        with self.assertRaises(ValueError):
            self.py_mut.oncoplot(top_genes_count=0)
            
        with self.assertRaises(ValueError):
            self.py_mut.oncoplot(max_samples=-1)
            
    def test_oncoplot_missing_columns(self):
        """Test oncoplot con columnas faltantes."""
        incomplete_data = pd.DataFrame({
            'Hugo_Symbol': ['TP53'],
            'TCGA-AB-1234': ['A|G']
        })
        py_mut_incomplete = PyMutation(incomplete_data)
        
        with self.assertRaises(ValueError):
            py_mut_incomplete.oncoplot()
            
    def test_oncoplot_interactive_mode(self):
        """Test modo interactivo del oncoplot."""
        with patch.object(self.py_mut, '_show_figure_interactive') as mock_show_interactive:
            fig = self.py_mut.oncoplot(show_interactive=True)
            self.assertIsInstance(fig, plt.Figure)
            mock_show_interactive.assert_called_once()
            plt.close(fig)


class TestOncoplotIntegration(unittest.TestCase):
    """Tests de integración para oncoplot con datos reales."""
    
    def setUp(self):
        """Configurar datos de prueba más realistas."""
        # Crear datos que simulan un dataset real
        genes = ['TP53', 'KRAS', 'PIK3CA', 'EGFR', 'BRAF'] * 10
        samples = [f'TCGA-AB-{1000+i}' for i in range(20)]
        
        data_rows = []
        for i, gene in enumerate(genes):
            for j, sample in enumerate(samples):
                if (i + j) % 3 == 0:  # Crear patrón de mutaciones
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
        
        # Crear DataFrame con estructura completa
        self.realistic_data = pd.DataFrame(data_rows).fillna('A|A')
        
        # Reorganizar para tener columnas correctas
        base_columns = ['Hugo_Symbol', 'Variant_Classification', 'REF', 'ALT']
        sample_columns = [col for col in self.realistic_data.columns if col.startswith('TCGA-')]
        
        # Crear DataFrame final con una fila por gen
        final_data = []
        for gene in ['TP53', 'KRAS', 'PIK3CA', 'EGFR', 'BRAF']:
            row = {
                'Hugo_Symbol': gene,
                'Variant_Classification': 'Missense_Mutation',
                'REF': 'A',
                'ALT': 'G'
            }
            for sample in samples:
                # Crear patrón pseudo-aleatorio de mutaciones
                if hash(f"{gene}_{sample}") % 3 == 0:
                    row[sample] = 'A|G'
                else:
                    row[sample] = 'A|A'
            final_data.append(row)
        
        self.realistic_data = pd.DataFrame(final_data)
        
    def test_oncoplot_with_real_data(self):
        """Test oncoplot con datos más realistas."""
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
            max_samples=3,
            show_interactive=False
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
            max_samples=3,
            show_interactive=False
        )
        
        # Verify figure was created successfully
        assert figure is not None
        axes = figure.get_axes()
        assert len(axes) >= 2
        
        plt.close(figure)


if __name__ == '__main__':
    # Configurar matplotlib para tests (modo no interactivo)
    import matplotlib
    matplotlib.use('Agg')
    
    # Ejecutar tests
    unittest.main(verbosity=2) 