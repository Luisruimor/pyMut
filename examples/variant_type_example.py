#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Ejemplo de visualización de tipos de variantes.

Este script muestra cómo utilizar la visualización de tipos de variantes
de forma independiente, sin necesidad de utilizar toda la biblioteca pyMut.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
from pathlib import Path

def extract_variant_type(funcotation_str):
    """
    Recibe una cadena del campo FUNCOTATION y extrae el valor del tipo de variante,
    que se asume es el octavo campo (índice 7) de la cadena separada por "|".
    """
    # Asegurarse de trabajar con una cadena
    funcotation_str = str(funcotation_str).strip()
    
    # Eliminar un posible prefijo "FUNCOTATION="
    if funcotation_str.startswith("FUNCOTATION="):
        funcotation_str = funcotation_str[len("FUNCOTATION="):].strip()
        
    # Quitar los corchetes inicial y final, si están presentes
    if funcotation_str.startswith('[') and funcotation_str.endswith(']'):
        funcotation_str = funcotation_str[1:-1]
    
    # Separar la cadena por "|" y extraer el octavo elemento (índice 7)
    fields = funcotation_str.split('|')
    if len(fields) > 7:
        return fields[7].strip()
    else:
        return None

def plot_variant_type(tsv_file,
                      variant_column="Variant_Type",
                      funcotation_column="FUNCOTATION",
                      output_file="variant_type.png"):
    """
    Lee un archivo TSV con anotaciones de variantes y genera un diagrama de barras horizontal
    mostrando el recuento por cada tipo de variante.
    
    Se espera que la información de tipo de variante esté en la columna 'variant_column'.
    Si no existe, se intenta extraer desde la columna 'funcotation_column' usando la función extract_variant_type.
    """
    try:
        # Inicializar un diccionario para acumular los recuentos
        variant_counts = {}

        # Leer el archivo en fragmentos (chunks) para archivos muy grandes
        for chunk in pd.read_csv(tsv_file, sep="\t", comment='#', engine='python', chunksize=1000):
            # Verificar si la columna de tipo de variante existe en el fragmento
            if variant_column not in chunk.columns:
                if funcotation_column in chunk.columns:
                    print(f"La columna '{variant_column}' no se encontró. Se extraerá desde '{funcotation_column}'.")
                    # Crear la columna aplicando la función de extracción
                    chunk[variant_column] = chunk[funcotation_column].apply(
                        lambda x: extract_variant_type(x) if pd.notnull(x) else None
                    )
                else:
                    raise ValueError(f"Ninguna de las columnas '{variant_column}' ni '{funcotation_column}' se encontró en el archivo.")
            
            # Eliminar filas donde no se pudo extraer el tipo de variante
            chunk = chunk[chunk[variant_column].notnull()]
            
            # Contar los tipos de variantes encontrados en este fragmento
            for variant in chunk[variant_column]:
                variant_counts[variant] = variant_counts.get(variant, 0) + 1

        # Ordenar los recuentos de menor a mayor (para graficar de abajo hacia arriba)
        variant_counts = dict(sorted(variant_counts.items(), key=lambda item: item[1]))

        # Crear la gráfica de barras horizontal
        plt.figure(figsize=(6, 5))
        # Colores específicos (si el número de categorías excede la lista, se utilizará un colormap)
        colors = ['#D3C4E7', '#FFFACD', '#87CEEB']
        if len(variant_counts) > len(colors):
            cmap = plt.colormaps['tab20']
            colors = cmap(range(len(variant_counts)))
        else:
            colors = colors[:len(variant_counts)]
            
        bars = plt.barh(list(variant_counts.keys()), list(variant_counts.values()), color=colors)
        
        # Configuración del título y etiquetas de los ejes
        plt.title("Variant Type", fontsize=14, loc="center")
        plt.xlabel("")  # Se elimina el título del eje X
        plt.ylabel("")  # Se elimina el título del eje Y
        
        # Configuración de las etiquetas del eje Y en cursiva
        plt.yticks(fontsize=12, fontstyle='italic')
        
        # Añadir etiquetas numéricas al final de cada barra
        for bar in bars:
            plt.text(bar.get_width() + 50,
                     bar.get_y() + bar.get_height() / 2,
                     f'{int(bar.get_width())}', va='center', fontsize=10)
        
        # Ajustar márgenes y eliminar algunos bordes para un estilo más limpio
        plt.subplots_adjust(left=0.25, right=0.95, top=0.9, bottom=0.1)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['left'].set_visible(False)
        
        # Guardar la gráfica en el archivo de salida
        plt.tight_layout()
        plt.savefig(output_file)
        print(f"Gráfica guardada en: {output_file}")
        
        # Mostrar la gráfica interactivamente
        plt.show()
    
    except Exception as e:
        print(f"Error: {e}")

# Ejemplo de uso con el archivo de muestra incluido en pyMut
if __name__ == "__main__":
    # Obtener la ruta absoluta al directorio del proyecto
    project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    
    # Ruta al archivo de ejemplo
    tsv_file_path = os.path.join(project_dir, "src", "pyMut", "data", "examples", "sample_mutations.tsv")
    
    # Verificar si existe el archivo
    if not os.path.exists(tsv_file_path):
        print(f"Error: El archivo {tsv_file_path} no existe.")
        exit(1)
        
    plot_variant_type(
        tsv_file_path,
        variant_column="Variant_Type",
        funcotation_column="FUNCOTATION",
        output_file="variant_type_example.png"
    ) 