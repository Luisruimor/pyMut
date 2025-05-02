#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Ejemplo básico de uso de la biblioteca pyMut.

Este script muestra cómo utilizar la biblioteca pyMut para generar
visualizaciones de resumen a partir de datos de mutaciones genéticas.
"""

import os
import sys
import matplotlib.pyplot as plt
import pandas as pd

# Añadir el directorio raíz al path para importar pyMut
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from src.pyMut import PyMutation
from src.pyMut.utils.data_processing import read_tsv

# Ruta al archivo de ejemplo
EXAMPLE_FILE = os.path.abspath(os.path.join(os.path.dirname(__file__), 
                                          '../../src/pyMut/data/examples/sample_mutations.tsv'))

def main():
    # Lectura de datos
    print(f"Leyendo datos de {EXAMPLE_FILE}")
    data = read_tsv(EXAMPLE_FILE)
    
    # Creación de un objeto PyMutation
    print("Creando objeto PyMutation")
    pyMut = PyMutation(data)
    
    # Generación de un gráfico de resumen
    print("Generando gráfico de resumen")
    fig = pyMut.summary_plot(title="Resumen de mutaciones de ejemplo")
    fig.savefig("summary_example.png")
    
    print("Ejemplo generado correctamente")

if __name__ == "__main__":
    main() 