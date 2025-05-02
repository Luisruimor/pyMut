#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Ejemplo sencillo de uso de pyMut para visualizaciones de resumen.

Este script muestra cómo utilizar la biblioteca pyMut para generar
visualizaciones de resumen a partir de un archivo TSV.

- Gráfico de Resumen: Incluye múltiples visualizaciones:
  * Variant Classification: Distribución de clasificaciones de variantes
  * Variant Type: Distribución de tipos de variantes (SNP, INS, DEL, etc.)
  * SNV Class: Distribución de clases de SNV (cambios nucleotídicos como A>G, C>T, etc.)
"""

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import os

# Importar pyMut
from pyMut import PyMutation

# Obtener la ruta absoluta al directorio del proyecto
project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Ruta al archivo de ejemplo
EXAMPLE_FILE = os.path.join(project_dir, "src", "pyMut", "data", "examples", "sample_mutations.tsv")

# Verificar si existe el archivo
if not os.path.exists(EXAMPLE_FILE):
    print(f"Error: El archivo {EXAMPLE_FILE} no existe.")
    exit(1)

# Leer los datos
data = pd.read_csv(EXAMPLE_FILE, sep="\t")
print(f"Se leyeron {len(data)} filas de datos.")

# Crear un objeto PyMutation
pyMut = PyMutation(data)

# Generar un gráfico de resumen (incluye Variant Classification, Variant Type y SNV Class)
print("Generando gráfico de resumen...")
fig = pyMut.summary_plot(title="Resumen de mutaciones")
fig.savefig("summary_example.png")
print("Gráfico de resumen guardado en summary_example.png")

# También podemos mostrar el gráfico interactivamente
plt.show() 