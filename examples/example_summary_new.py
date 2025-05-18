"""
Ejemplo básico para visualización de mutaciones con pyMut.

Este ejemplo muestra cómo generar visualizaciones de datos de mutaciones genéticas
utilizando pyMut, mostrando tanto el resumen completo como visualizaciones individuales
de la forma más simple posible.
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
    mutation_data = pd.read_csv(tsv_file, sep='\t')
    
    # Crear objeto PyMutation
    print("Creando objeto PyMutation...")
    py_mut = PyMutation(mutation_data)
    
    # 1. Generar el gráfico de resumen completo
    print("\n1. Generando gráfico de resumen completo...")
    summary_fig = py_mut.summary_plot()
    # Parámetros opcionales comentados:
    # summary_fig = py_mut.summary_plot(
         # title="Ejemplo del plot summary",     # Título personalizado
         # figsize=(16, 12),                     # Tamaño de la figura
         # max_samples=50,                       # Máximo de muestras a mostrar
         # top_genes_count=15,                   # Número de genes en Top Mutated Genes (por defecto 10)
         # show_interactive=True                 # Mostrar interactivamente
    # )
    
    # Guardar la figura
    output_path = os.path.join(os.path.dirname(__file__), 'summary_plot.png')
    summary_fig.savefig(output_path)
    # Opciones adicionales para guardar:
    # summary_fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Gráfico de resumen guardado en: {output_path}")
    
    # 2. Variant Classification Plot
    print("\n2. Generando gráfico de clasificación de variantes...")
    vc_fig = py_mut.variant_classification_plot()
    # Parámetros opcionales comentados:
    # vc_fig = py_mut.variant_classification_plot(
    #     figsize=(12, 6),                      # Tamaño de la figura
    #     title="Variant Classification",       # Título personalizado
    #     show_interactive=True                 # Mostrar interactivamente
    # )
    
    vc_output_path = os.path.join(os.path.dirname(__file__), 'variant_classification.png')
    vc_fig.savefig(vc_output_path)
    print(f"Gráfico guardado en: {vc_output_path}")
    
    # 3. Variant Type Plot
    print("\n3. Generando gráfico de tipos de variantes...")
    vt_fig = py_mut.variant_type_plot()
    # Parámetros opcionales comentados:
    # vt_fig = py_mut.variant_type_plot(
    #      figsize=(12, 6),                      # Tamaño de la figura
    #      title="Variant Type",                 # Título personalizado
    #      show_interactive=True                 # Mostrar interactivamente
    #  )
    
    vt_output_path = os.path.join(os.path.dirname(__file__), 'variant_type.png')
    vt_fig.savefig(vt_output_path)
    print(f"Gráfico guardado en: {vt_output_path}")
    
    # 4. SNV Class Plot
    print("\n4. Generando gráfico de clases de SNV...")
    snv_fig = py_mut.snv_class_plot()
    # Parámetros opcionales comentados:
    # snv_fig = py_mut.snv_class_plot(
    #     figsize=(12, 6),                      # Tamaño de la figura
    #     title="SNV Class",                    # Título personalizado
    #     ref_column="REF",                     # Columna con alelo de referencia
    #     alt_column="ALT",                     # Columna con alelo alternativo
    #     show_interactive=True                 # Mostrar interactivamente
    # )
    
    snv_output_path = os.path.join(os.path.dirname(__file__), 'snv_class.png')
    snv_fig.savefig(snv_output_path)
    print(f"Gráfico guardado en: {snv_output_path}")
    
    # 5. Variants per Sample Plot (TMB)
    print("\n5. Generando gráfico de variantes por muestra (TMB)...")
    tmb_fig = py_mut.variants_per_sample_plot()
    # Parámetros opcionales comentados:
    # tmb_fig = py_mut.variants_per_sample_plot(
         # figsize=(12, 6),                         # Tamaño de la figura
         # max_samples=50,                          # Máximo de muestras a mostrar
         # title="Variants per Sample",             # Título personalizado
         # variant_column="Variant_Classification", # Columna con clasificación de variante
         # sample_column="Tumor_Sample_Barcode",    # Columna con ID de muestra
         # show_interactive=True                    # Mostrar interactivamente
    # )
    
    tmb_output_path = os.path.join(os.path.dirname(__file__), 'variants_per_sample.png')
    tmb_fig.savefig(tmb_output_path)
    print(f"Gráfico guardado en: {tmb_output_path}")
    
    # 6. Variant Classification Summary Plot (Boxplot)
    print("\n6. Generando gráfico de resumen de clasificación de variantes (Boxplot)...")
    vcs_fig = py_mut.variant_classification_summary_plot()
    # Parámetros opcionales comentados:
    # vcs_fig = py_mut.variant_classification_summary_plot(
    #     figsize=(12, 6),                      # Tamaño de la figura
    #     title="Variant Classification Summary", # Título personalizado
    #     variant_column="Variant_Classification", # Columna con clasificación de variante
    #     sample_column="Tumor_Sample_Barcode", # Columna con ID de muestra
    #     show_interactive=True                 # Mostrar interactivamente
    # )
    
    vcs_output_path = os.path.join(os.path.dirname(__file__), 'variant_classification_summary.png')
    vcs_fig.savefig(vcs_output_path)
    print(f"Gráfico guardado en: {vcs_output_path}")
    
    # 7a. Top Mutated Genes Plot (modo variants)
    print("\n7a. Generando gráfico de genes más mutados (modo variants)...")
    tmg_variants_fig = py_mut.top_mutated_genes_plot(mode="variants")
    # Parámetros opcionales comentados:
    # tmg_variants_fig = py_mut.top_mutated_genes_plot(
    #     mode="variants",                      # Modo de conteo (obligatorio)
    #     figsize=(12, 6),                      # Tamaño de la figura
    #     title="Top Mutated Genes",            # Título personalizado
    #     variant_column="Variant_Classification", # Columna con clasificación de variante
    #     gene_column="Hugo_Symbol",            # Columna con símbolo del gen
    #     sample_column="Tumor_Sample_Barcode", # Columna con ID de muestra
    #     count=10,                             # Número de genes a mostrar (por defecto 10, muestra todos si hay menos)
    #     show_interactive=True                 # Mostrar interactivamente
    # )
    
    tmg_variants_path = os.path.join(os.path.dirname(__file__), 'top_mutated_genes_variants.png')
    tmg_variants_fig.savefig(tmg_variants_path)
    print(f"Gráfico guardado en: {tmg_variants_path}")
    
    # 7b. Top Mutated Genes Plot (modo samples)
    print("\n7b. Generando gráfico de genes más mutados (modo samples)...")
    tmg_samples_fig = py_mut.top_mutated_genes_plot(mode="samples")
    # Parámetros opcionales comentados:
    # tmg_samples_fig = py_mut.top_mutated_genes_plot(
    #     mode="samples",                       # Modo de conteo (obligatorio)
    #     figsize=(12, 6),                      # Tamaño de la figura
    #     title="Top Mutated Genes",            # Título personalizado
    #     variant_column="Variant_Classification", # Columna con clasificación de variante
    #     gene_column="Hugo_Symbol",            # Columna con símbolo del gen
    #     sample_column="Tumor_Sample_Barcode", # Columna con ID de muestra
    #     count=15,                             # Número de genes a mostrar
    #     show_interactive=True                 # Mostrar interactivamente
    # )
    
    tmg_samples_path = os.path.join(os.path.dirname(__file__), 'top_mutated_genes_samples.png')
    tmg_samples_fig.savefig(tmg_samples_path)
    print(f"Gráfico guardado en: {tmg_samples_path}")
    
    print("\nEl ejemplo de generación de gráficos ha finalizado.")
    
    # Para mostrar las visualizaciones interactivamente, descomenta las siguientes líneas:
    # plt.ion()  # Activar modo interactivo
    # plt.show()

if __name__ == "__main__":
    main()
