"""
Módulo para el procesamiento de datos de mutaciones.

Este módulo contiene funciones para procesar y transformar datos de mutaciones
genéticas para su posterior visualización en gráficos de resumen.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Union, Optional, Tuple
import os


def extract_variant_classification(funcotation_str):
    """
    Recibe una cadena del campo FUNCOTATION y extrae el valor de la clasificación de variante,
    que se asume es el sexto campo (índice 5) de la cadena separada por "|".
    
    Args:
        funcotation_str: Cadena del campo FUNCOTATION.
        
    Returns:
        Clasificación de la variante extraída o None si no se puede extraer.
    """
    # Asegurarse de trabajar con una cadena (en caso de que sea otro tipo)
    funcotation_str = str(funcotation_str).strip()
    
    # Si la cadena incluye un prefijo como "FUNCOTATION=", eliminarlo
    if funcotation_str.startswith("FUNCOTATION="):
        funcotation_str = funcotation_str[len("FUNCOTATION="):].strip()
        
    # Quitar los corchetes inicial y final, si están presentes
    if funcotation_str.startswith('[') and funcotation_str.endswith(']'):
        funcotation_str = funcotation_str[1:-1]
    
    # Separar la cadena por "|" y extraer el sexto elemento (índice 5)
    fields = funcotation_str.split('|')
    if len(fields) > 5:
        return fields[5].strip()
    else:
        return None


def extract_variant_classifications(data: pd.DataFrame, 
                                  variant_column: str = "Variant_Classification",
                                  funcotation_column: str = "FUNCOTATION") -> pd.DataFrame:
    """
    Extrae la clasificación de variante desde el campo FUNCOTATION si no existe una columna
    específica para ello.
    
    Args:
        data: DataFrame con los datos de mutaciones.
        variant_column: Nombre de la columna que debe contener la clasificación de variante.
        funcotation_column: Nombre de la columna que contiene la cadena FUNCOTATION.
        
    Returns:
        DataFrame con la columna de clasificación de variante añadida o sin cambios si ya existía.
    """
    # Crear una copia para no modificar el original
    result = data.copy()
    
    # Verificar si la columna de clasificación ya existe
    if variant_column not in result.columns:
        # Si existe la columna FUNCOTATION, extraer de ahí
        if funcotation_column in result.columns:
            print(f"La columna '{variant_column}' no se encontró. Se extraerá desde '{funcotation_column}'.")
            # Crear la columna aplicando la función de extracción
            result[variant_column] = result[funcotation_column].apply(
                lambda x: extract_variant_classification(x) if pd.notnull(x) else None
            )
            
            # Verificar si se obtuvieron datos válidos
            if result[variant_column].isna().all():
                print(f"No se pudieron extraer datos de clasificación de variante desde '{funcotation_column}'.")
                result[variant_column] = "Unknown"
        else:
            print(f"Ninguna de las columnas '{variant_column}' ni '{funcotation_column}' se encontró en el DataFrame.")
            # Crear una columna con valor desconocido
            result[variant_column] = "Unknown"
            
    # Rellenar los valores NaN con "Unknown"
    result[variant_column] = result[variant_column].fillna("Unknown")
    
    return result


def extract_variant_type(funcotation_str):
    """
    Recibe una cadena del campo FUNCOTATION y extrae el valor del tipo de variante,
    que se asume es el octavo campo (índice 7) de la cadena separada por "|".
    
    Args:
        funcotation_str: Cadena del campo FUNCOTATION.
        
    Returns:
        Tipo de variante extraído o None si no se puede extraer.
    """
    # Asegurarse de trabajar con una cadena (en caso de que sea otro tipo)
    funcotation_str = str(funcotation_str).strip()
    
    # Si la cadena incluye un prefijo como "FUNCOTATION=", eliminarlo
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


def extract_variant_types(data: pd.DataFrame, 
                       variant_column: str = "Variant_Type",
                       funcotation_column: str = "FUNCOTATION") -> pd.DataFrame:
    """
    Extrae el tipo de variante desde el campo FUNCOTATION si no existe una columna
    específica para ello.
    
    Args:
        data: DataFrame con los datos de mutaciones.
        variant_column: Nombre de la columna que debe contener el tipo de variante.
        funcotation_column: Nombre de la columna que contiene la cadena FUNCOTATION.
        
    Returns:
        DataFrame con la columna de tipo de variante añadida o sin cambios si ya existía.
    """
    # Crear una copia para no modificar el original
    result = data.copy()
    
    # Verificar si la columna de tipo de variante ya existe
    if variant_column not in result.columns:
        # Si existe la columna FUNCOTATION, extraer de ahí
        if funcotation_column in result.columns:
            print(f"La columna '{variant_column}' no se encontró. Se extraerá desde '{funcotation_column}'.")
            # Crear la columna aplicando la función de extracción
            result[variant_column] = result[funcotation_column].apply(
                lambda x: extract_variant_type(x) if pd.notnull(x) else None
            )
            
            # Verificar si se obtuvieron datos válidos
            if result[variant_column].isna().all():
                print(f"No se pudieron extraer datos de tipos de variante desde '{funcotation_column}'.")
                result[variant_column] = "Unknown"
        else:
            print(f"Ninguna de las columnas '{variant_column}' ni '{funcotation_column}' se encontró en el DataFrame.")
            # Crear una columna con valor desconocido
            result[variant_column] = "Unknown"
    
    # Rellenar los valores NaN con "Unknown"
    result[variant_column] = result[variant_column].fillna("Unknown")
    
    return result


def extract_genome_change(funcotation_str):
    """
    Recibe una cadena del campo FUNCOTATION y extrae el valor del cambio genómico,
    que se asume es el duodécimo campo (índice 11) de la cadena separada por "|".
    
    Args:
        funcotation_str: Cadena del campo FUNCOTATION.
        
    Returns:
        Cambio genómico extraído o None si no se puede extraer.
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


def extract_genome_changes(data: pd.DataFrame, 
                        genome_change_column: str = "Genome_Change",
                        funcotation_column: str = "FUNCOTATION") -> pd.DataFrame:
    """
    Extrae el cambio genómico desde el campo FUNCOTATION si no existe una columna
    específica para ello.
    
    Args:
        data: DataFrame con los datos de mutaciones.
        genome_change_column: Nombre de la columna que debe contener el cambio genómico.
        funcotation_column: Nombre de la columna que contiene la cadena FUNCOTATION.
        
    Returns:
        DataFrame con la columna de cambio genómico añadida o sin cambios si ya existía.
    """
    # Crear una copia para no modificar el original
    result = data.copy()
    
    # Verificar si la columna de cambio genómico ya existe
    if genome_change_column not in result.columns:
        # Si existe la columna FUNCOTATION, extraer de ahí
        if funcotation_column in result.columns:
            print(f"La columna '{genome_change_column}' no se encontró. Se extraerá desde '{funcotation_column}'.")
            # Crear la columna aplicando la función de extracción
            result[genome_change_column] = result[funcotation_column].apply(
                lambda x: extract_genome_change(x) if pd.notnull(x) else None
            )
            
            # Verificar si se obtuvieron datos válidos
            if result[genome_change_column].isna().all():
                print(f"No se pudieron extraer datos de cambio genómico desde '{funcotation_column}'.")
                result[genome_change_column] = "Unknown"
        else:
            print(f"Ninguna de las columnas '{genome_change_column}' ni '{funcotation_column}' se encontró en el DataFrame.")
            # Crear una columna con valor desconocido
            result[genome_change_column] = "Unknown"
    
    # Rellenar los valores NaN con "Unknown"
    result[genome_change_column] = result[genome_change_column].fillna("Unknown")
    
    return result


def read_tsv(file_path: str) -> pd.DataFrame:
    """
    Lee un archivo TSV y devuelve un DataFrame con los datos de mutaciones.
    
    Args:
        file_path: Ruta al archivo TSV.
        
    Returns:
        DataFrame con los datos de mutaciones.
    
    Raises:
        FileNotFoundError: Si el archivo no existe.
        pd.errors.EmptyDataError: Si el archivo está vacío.
        pd.errors.ParserError: Si hay problemas al analizar el archivo.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"El archivo '{file_path}' no existe.")
        
    try:
        # Primero intentamos leer sin comentarios
        data = pd.read_csv(file_path, sep='\t')
    except (pd.errors.ParserError, pd.errors.EmptyDataError) as e:
        # Si falla por un error de análisis, intentamos con el parámetro comment
        try:
            data = pd.read_csv(file_path, sep='\t', comment='#', engine='python')
        except Exception as err:
            raise ValueError(f"No se pudo leer el archivo: {str(err)}")
    except Exception as e:
        raise ValueError(f"Error al leer el archivo: {str(e)}")
        
    return data