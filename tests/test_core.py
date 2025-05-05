#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

if __name__ == "__main__":
    unittest.main()
