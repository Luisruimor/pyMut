"""
Ejemplo sencillo de visualización individual de clases de SNV.

Este ejemplo muestra cómo generar un gráfico de clases de SNV
de forma individual utilizando la API de pyMut.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Importar pyMut
from pyMut import PyMutation

def main():
    # Obtener la ruta absoluta al directorio del proyecto
    project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    
    # Ruta al archivo de ejemplo
    tsv_file_path = os.path.join(project_dir, "src", "pyMut", "data", "examples", "tcga_laml_converted.tsv")
    
    # Verificar si existe el archivo
    if not os.path.exists(tsv_file_path):
        print(f"Error: El archivo {tsv_file_path} no existe.")
        exit(1)
    
    # Leer los datos
    print(f"Leyendo datos de {tsv_file_path}...")
    data = pd.read_csv(tsv_file_path, sep="\t")
    print(f"Se leyeron {len(data)} filas de datos.")
    
    # Crear un objeto PyMutation
    pyMut = PyMutation(data)
    
    # Generar gráfico de clases de SNV
    print("Generando gráfico de clases de SNV...")
    fig = pyMut.snv_class_plot(
        figsize=(10, 6),
        title="Clases de SNV",
        ref_column="REF",
        alt_column="ALT"
    )
    
    # Guardar el gráfico
    output_file = "snv_class_individual_example.png"
    fig.savefig(output_file)
    print(f"Gráfico guardado en: {output_file}")
    
    # Mostrar el gráfico interactivamente
    plt.show()

if __name__ == "__main__":
    main() 