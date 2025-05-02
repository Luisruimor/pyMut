#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Ejemplo de visualización de clases de SNV (Single Nucleotide Variant).

Este script muestra cómo utilizar la visualización de clases de SNV
de forma independiente, sin necesidad de utilizar toda la biblioteca pyMut.
"""

import pandas as pd
import matplotlib.pyplot as plt
import re
import os
from pathlib import Path

def extract_genome_change(funcotation_str):
    """
    Recibe una cadena del campo FUNCOTATION y extrae el valor del cambio genómico,
    que se asume es el duodécimo campo (índice 11) de la cadena separada por "|".
    """
    # Convertir a cadena y quitar espacios en blanco
    funcotation_str = str(funcotation_str).strip()
    
    # Eliminar el posible prefijo "FUNCOTATION="
    if funcotation_str.startswith("FUNCOTATION="):
        funcotation_str = funcotation_str[len("FUNCOTATION="):].strip()
        
    # Quitar los corchetes de apertura y cierre si están presentes
    if funcotation_str.startswith('[') and funcotation_str.endswith(']'):
        funcotation_str = funcotation_str[1:-1]
    
    # Separar la cadena por "|" y extraer el campo correspondiente (índice 11)
    fields = funcotation_str.split('|')
    if len(fields) > 11:
        return fields[11].strip()
    else:
        return None

def plot_snv_class(tsv_file,
                   snv_column="Genome_Change",
                   ref_allele_column="Reference_Allele",
                   tumor_allele_column="Tumor_Allele",
                   funcotation_column="FUNCOTATION",
                   output_file="snv_class_example.png"):
    """
    Lee un archivo TSV con anotaciones de variantes y genera un diagrama de barras horizontal
    mostrando el recuento por cada clase de SNV (Single Nucleotide Variant).
    
    Se espera que la información de cambio genómico esté en la columna 'snv_column'.
    Si no existe, se intenta extraer desde 'funcotation_column' usando la función extract_genome_change.
    Si tampoco existe 'funcotation_column', se intenta generar a partir de 'ref_allele_column' y 'tumor_allele_column'.
    
    A partir del valor en snv_column se extrae la SNV Class aplicando la expresión regular:
       ([A-Z]>[A-Z])$
    que busca una transformación tipo "G>A" al final de la cadena.
    """
    try:
        # Leer todos los datos de una vez (en un caso real con archivos grandes, usar chunksize)
        df = pd.read_csv(tsv_file, sep="\t", comment='#', engine='python')
        
        # Verificar qué estrategia usar para obtener la SNV Class
        if snv_column in df.columns:
            # Usar la columna de cambio genómico existente
            print(f"Usando la columna '{snv_column}' para extraer la SNV Class.")
            df['SNV_Class'] = df[snv_column].str.extract(r'([A-Z]>[A-Z])$')
        elif funcotation_column in df.columns:
            # Extraer el cambio genómico desde FUNCOTATION
            print(f"Extrayendo el cambio genómico desde '{funcotation_column}'.")
            df[snv_column] = df[funcotation_column].apply(
                lambda x: extract_genome_change(x) if pd.notnull(x) else None
            )
            df['SNV_Class'] = df[snv_column].str.extract(r'([A-Z]>[A-Z])$')
        elif ref_allele_column in df.columns and tumor_allele_column in df.columns:
            # Generar la SNV Class directamente desde los alelos
            print(f"Generando la SNV Class desde '{ref_allele_column}' y '{tumor_allele_column}'.")
            df['SNV_Class'] = df[ref_allele_column] + '>' + df[tumor_allele_column]
            # Filtrar para mantener solo los cambios de nucléotido único (un solo carácter a cada lado)
            df = df[df['SNV_Class'].str.match(r'^[A-Z]>[A-Z]$')]
        else:
            raise ValueError("No se encontraron columnas necesarias para generar la SNV Class.")
        
        # Eliminar filas sin información extraída
        df = df[df['SNV_Class'].notnull()]
        
        if len(df) == 0:
            raise ValueError("No se pudieron extraer clases de SNV de los datos.")
        
        # Contar las clases de SNV
        snv_counts = df['SNV_Class'].value_counts().to_dict()
        
        # Ordenar los recuentos de menor a mayor para que la gráfica se muestre de abajo hacia arriba
        snv_counts = dict(sorted(snv_counts.items(), key=lambda item: item[1]))
        
        # Crear la gráfica de barras horizontal
        plt.figure(figsize=(8, 6))
        # Definir colores específicos para cada clase
        colors = ['#FF8C00', '#9ACD32', '#FFD700', '#FF4500', '#4169E1', '#1E90FF']
        # Si hay más categorías que colores definidos, se puede optar por un colormap
        if len(snv_counts) > len(colors):
            cmap = plt.colormaps['tab20']
            colors = cmap(range(len(snv_counts)))
        else:
            colors = colors[:len(snv_counts)]
            
        bars = plt.barh(list(snv_counts.keys()), list(snv_counts.values()), color=colors)
        
        # Configurar título y eliminar etiqueta del eje Y
        plt.title("SNV Class", fontsize=14, loc="center")
        plt.ylabel("")
        
        # Añadir etiquetas numéricas en cada barra (en cursiva)
        for bar, label in zip(bars, list(snv_counts.values())):
            plt.text(bar.get_width() + 10,
                     bar.get_y() + bar.get_height() / 2,
                     f'{int(label)}', va='center', fontsize=10, fontstyle='italic')
        
        # Ajustar márgenes y quitar recuadros innecesarios
        plt.subplots_adjust(left=0.3, right=0.95, top=0.9, bottom=0.1)
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
        
    plot_snv_class(
        tsv_file_path,
        snv_column="Genome_Change",
        ref_allele_column="Reference_Allele",
        tumor_allele_column="Tumor_Allele",
        funcotation_column="FUNCOTATION",
        output_file="snv_class_example.png"
    ) 