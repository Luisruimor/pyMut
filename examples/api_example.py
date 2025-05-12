"""
Ejemplo de uso de la API de pyMut.

Este ejemplo muestra cómo utilizar la API simplificada de pyMut para generar
tanto gráficos de resumen completos como visualizaciones individuales.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

# Importar pyMut
from pyMut import PyMutation

def main():
    """Función principal para demostrar el uso de la API de pyMut."""
    
    # Obtener la ruta absoluta al directorio del proyecto
    project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    
    # Ruta al archivo de ejemplo
    tsv_file = os.path.join(project_dir, "src", "pyMut", "data", "examples", "tcga_laml_converted.tsv")
    
    # Verificar si existe el archivo
    if not os.path.exists(tsv_file):
        print(f"Error: El archivo {tsv_file} no existe.")
        exit(1)
    
    # Leer los datos
    print(f"Leyendo datos de {tsv_file}...")
    data = pd.read_csv(tsv_file, sep="\t")
    print(f"Se leyeron {len(data)} filas de datos.")
    
    # Crear un objeto PyMutation
    print("Creando objeto PyMutation...")
    pyMut = PyMutation(data)
    
    # 1. Generar el gráfico de resumen completo
    print("\n1. Generando gráfico de resumen completo...")
    summary_fig = pyMut.summary_plot(title="Resumen de Mutaciones")
    summary_fig.savefig("summary_api_example.png")
    print("Gráfico de resumen guardado en: summary_api_example.png")
    
    # 2. Generar visualizaciones individuales
    print("\n2. Generando visualizaciones individuales...")
    
    # 2.1 Clasificación de variantes
    print("\n2.1 Generando gráfico de clasificación de variantes...")
    vc_fig = pyMut.variant_classification_plot(title="Clasificación de Variantes")
    vc_fig.savefig("variant_classification_api_example.png")
    print("Gráfico guardado en: variant_classification_api_example.png")
    
    # 2.2 Tipos de variantes
    print("\n2.2 Generando gráfico de tipos de variantes...")
    vt_fig = pyMut.variant_type_plot(title="Tipos de Variantes")
    vt_fig.savefig("variant_type_api_example.png")
    print("Gráfico guardado en: variant_type_api_example.png")
    
    # 2.3 Clases de SNV
    print("\n2.3 Generando gráfico de clases de SNV...")
    snv_fig = pyMut.snv_class_plot(title="Clases de SNV")
    snv_fig.savefig("snv_class_api_example.png")
    print("Gráfico guardado en: snv_class_api_example.png")
    
    print("\nTodos los gráficos han sido generados correctamente.")
    
    # Mostrar un gráfico como ejemplo (descomenta para mostrar interactivamente)
    # plt.figure(summary_fig.number)
    # plt.show()

if __name__ == "__main__":
    main() 