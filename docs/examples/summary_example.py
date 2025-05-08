#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Ejemplo de uso del gráfico de resumen de mutaciones.

Este script muestra cómo utilizar la biblioteca pyMut para generar
un gráfico de resumen a partir de datos de mutaciones genéticas.

Para ejecutar este ejemplo:
1. Asegúrate de estar en la raíz del proyecto pyMut
2. Ejecuta: python docs/examples/summary_example.py
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
                                          '../../src/pyMut/data/examples/tcga_laml_converted.tsv'))

def main():
    # Lectura de datos
    print(f"Leyendo datos de {EXAMPLE_FILE}")
    try:
        data = read_tsv(EXAMPLE_FILE)
        print(f"Se leyeron {len(data)} filas de datos.")
        
        # Crear una columna Variant_Classification si no existe
        if 'Variant_Classification' not in data.columns:
            if 'Variant_Classification' in data.columns:
                # Renombrar la columna para evitar confusiones
                data.rename(columns={'Variant_Classification': 'Variant_Classification'}, inplace=True)
            else:
                # Si no existe, usar la que ya tiene
                print("Usando la columna existente de clasificación de variantes.")
    
        # Creación de un objeto PyMutation
        print("Creando objeto PyMutation")
        pyMut = PyMutation(data)
        
        # Generación de un gráfico de resumen
        print("Generando gráfico de resumen")
        fig = pyMut.summary_plot(title="Resumen de mutaciones en el ejemplo")
        
        # Guardar y mostrar el gráfico
        output_file = "summary_example.png"
        fig.savefig(output_file)
        print(f"Gráfico de resumen guardado en: {output_file}")
        
        # También podemos mostrar el gráfico interactivamente
        # plt.show()  # Descomenta esta línea si quieres ver el gráfico interactivamente
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main() 