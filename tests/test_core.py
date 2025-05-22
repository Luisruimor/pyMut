"""
Pruebas para la clase PyMutation.

Este archivo contiene pruebas unitarias para la clase PyMutation
y sus métodos principales, enfocándose en la funcionalidad de gráficos de resumen.
"""

import os
import sys
import unittest
import pandas as pd
import matplotlib.pyplot as plt

# Añadir el directorio raíz al path para importar pyMut
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.pyMut import PyMutation
from src.pyMut.utils.data_processing import read_tsv

# Ruta al archivo de ejemplo
EXAMPLE_FILE = os.path.abspath(os.path.join(os.path.dirname(__file__), 
                                         '../src/pyMut/data/examples/tcga_laml_converted.tsv'))

class TestPyMutation(unittest.TestCase):
    """Pruebas para la clase PyMutation."""
    
    def setUp(self):
        """Configuración inicial para las pruebas."""
        self.data = read_tsv(EXAMPLE_FILE)
        self.pyMut = PyMutation(self.data)
    
    def test_initialization(self):
        """Prueba la inicialización de PyMutation."""
        self.assertIsInstance(self.pyMut, PyMutation)
        self.assertIsInstance(self.pyMut.data, pd.DataFrame)
        self.assertEqual(len(self.pyMut.data), len(self.data))
    
    def test_summary_plot(self):
        """Prueba el método summary_plot."""
        fig = self.pyMut.summary_plot()
        self.assertIsInstance(fig, plt.Figure)

    def test_variant_classification_plot(self):
        """Prueba el método variant_classification_plot."""
        fig = self.pyMut.variant_classification_plot()
        self.assertIsInstance(fig, plt.Figure)

    def test_variant_type_plot(self):
        """Prueba el método variant_type_plot."""
        fig = self.pyMut.variant_type_plot()
        self.assertIsInstance(fig, plt.Figure)

    def test_snv_class_plot(self):
        """Prueba el método snv_class_plot."""
        fig = self.pyMut.snv_class_plot()
        self.assertIsInstance(fig, plt.Figure)

    def test_variant_classification_summary_plot(self):
        """Prueba el método variant_classification_summary_plot."""
        fig = self.pyMut.variant_classification_summary_plot()
        self.assertIsInstance(fig, plt.Figure)

    def test_variants_per_sample_plot(self):
        """Prueba el método variants_per_sample_plot."""
        fig = self.pyMut.variants_per_sample_plot()
        self.assertIsInstance(fig, plt.Figure)

    def test_top_mutated_genes_plot(self):
        """Prueba el método top_mutated_genes_plot."""
        # Probar el modo "variants"
        fig = self.pyMut.top_mutated_genes_plot(mode="variants")
        self.assertIsInstance(fig, plt.Figure)
        
        # Probar el modo "samples"
        fig = self.pyMut.top_mutated_genes_plot(mode="samples")
        self.assertIsInstance(fig, plt.Figure)

    def test_validation_top_mutated_genes_plot(self):
        """Prueba la validación de parámetros en top_mutated_genes_plot."""
        # Verificar que se valide el parámetro 'count'
        with self.assertRaises(ValueError):
            self.pyMut.top_mutated_genes_plot(count=0)
        
        with self.assertRaises(ValueError):
            self.pyMut.top_mutated_genes_plot(count=-5)
        
        # Verificar que se valide el parámetro 'mode'
        with self.assertRaises(ValueError):
            self.pyMut.top_mutated_genes_plot(mode="invalid_mode")

    def test_empty_dataframe_initialization(self):
        """Prueba que se valide que el DataFrame no esté vacío."""
        with self.assertRaises(ValueError):
            PyMutation(pd.DataFrame())

    def test_non_dataframe_initialization(self):
        """Prueba que se valide que el parámetro sea un DataFrame."""
        with self.assertRaises(ValueError):
            PyMutation([1, 2, 3])  # Lista en lugar de DataFrame
        
        with self.assertRaises(ValueError):
            PyMutation("Not a DataFrame")  # String en lugar de DataFrame



if __name__ == "__main__":
    unittest.main()
