"""
Ejemplo mejorado para visualización de mutaciones con pyMut.

Este ejemplo muestra cómo generar visualizaciones de alta calidad de datos de mutaciones
genéticas utilizando pyMut, con opciones para ver tanto el resumen completo como
visualizaciones individuales.
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

# Añadir el directorio src al path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Importar PyMutation
from src.pyMut import PyMutation

def main():
    """
    Script principal para generar visualizaciones de mutaciones genéticas.
    """
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
    mutation_data = pd.read_csv(tsv_file, sep='\t')
    print(f"Datos cargados: {mutation_data.shape[0]} filas, {mutation_data.shape[1]} columnas")
    
    # Crear objeto PyMutation
    print("Creando objeto PyMutation...")
    py_mut = PyMutation(mutation_data)
    
    # 1. Generar el gráfico de resumen completo
    print("\n1. Generando gráfico de resumen completo...")
    summary_fig = py_mut.summary_plot(
        title="Resumen de Mutaciones", 
        show_interactive=False
    )
    
    # Guardar la figura con alta resolución
    output_path = os.path.join(os.path.dirname(__file__), 'summary_plot.png')
    summary_fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Gráfico de resumen guardado en: {output_path}")
    
    # 2. Generar visualizaciones individuales
    print("\n2. Generando visualizaciones individuales...")
    
    # 2.1 Variantes por muestra (TMB) - Esta se muestra por defecto
    print("\n2.1 Generando gráfico de variantes por muestra (TMB)...")
    tmb_fig = py_mut.variants_per_sample_plot(
        title="Carga Mutacional por Muestra (TMB)",
        figsize=(12, 6),  # Figura más ancha
        show_interactive=False
    )
    tmb_output_path = os.path.join(os.path.dirname(__file__), 'variants_per_sample.png')
    tmb_fig.savefig(tmb_output_path, dpi=300, bbox_inches='tight')
    print(f"Gráfico guardado en: {tmb_output_path}")
    
    # 2.2 Variant Classification Summary (Nuevo)
    print("\n2.2 Generando gráfico de resumen de clasificación de variantes (Boxplot)...")
    vcs_fig = py_mut.variant_classification_summary_plot(
        title="Resumen de Clasificación de Variantes por Muestra",
        figsize=(12, 6),
        show_interactive=False
    )
    vcs_output_path = os.path.join(os.path.dirname(__file__), 'variant_classification_summary.png')
    vcs_fig.savefig(vcs_output_path, dpi=300, bbox_inches='tight')
    print(f"Gráfico guardado en: {vcs_output_path}")
    
    # Descomenta las siguientes secciones para generar otros gráficos individuales:
    
    # # 2.3 Clasificación de variantes
    # print("\n2.3 Generando gráfico de clasificación de variantes...")
    # vc_fig = py_mut.variant_classification_plot(
    #     title="Clasificación de Variantes",
    #     figsize=(12, 6),  # Figura más ancha
    #     show_interactive=False
    # )
    # vc_output_path = os.path.join(os.path.dirname(__file__), 'variant_classification.png')
    # vc_fig.savefig(vc_output_path, dpi=300, bbox_inches='tight')
    # print(f"Gráfico guardado en: {vc_output_path}")
    
    # # 2.4 Tipos de variantes
    # print("\n2.4 Generando gráfico de tipos de variantes...")
    # vt_fig = py_mut.variant_type_plot(
    #     title="Tipos de Variantes",
    #     figsize=(12, 6),  # Figura más ancha
    #     show_interactive=False
    # )
    # vt_output_path = os.path.join(os.path.dirname(__file__), 'variant_type.png')
    # vt_fig.savefig(vt_output_path, dpi=300, bbox_inches='tight')
    # print(f"Gráfico guardado en: {vt_output_path}")
    
    # # 2.5 Clases de SNV
    # print("\n2.5 Generando gráfico de clases de SNV...")
    # snv_fig = py_mut.snv_class_plot(
    #     title="Clases de SNV",
    #     figsize=(12, 6),  # Figura más ancha
    #     show_interactive=False
    # )
    # snv_output_path = os.path.join(os.path.dirname(__file__), 'snv_class.png')
    # snv_fig.savefig(snv_output_path, dpi=300, bbox_inches='tight')
    # print(f"Gráfico guardado en: {snv_output_path}")
    
    print("\nEl ejemplo de generación de gráficos ha finalizado.")
    print("Descomenta las secciones correspondientes en el código para generar otras visualizaciones individuales.")
    
    # Mostrar las visualizaciones interactivamente
    print("\nMostrando visualizaciones interactivas...")
    plt.figure(summary_fig.number)
    plt.figure(tmb_fig.number)
    plt.figure(vcs_fig.number)
    plt.show()

if __name__ == "__main__":
    main() 