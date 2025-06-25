from typing import re
import pandas as pd

def formatear_rs(cadena):
    """
    Removes the 'rs' prefix from rsID codes in a pipe-separated string.
    
    Takes a string containing rsID codes separated by '|' and removes
    the 'rs' prefix from each code, returning only the numeric part.
    
    Parameters
    ----------
    cadena : str
        Pipe-separated string containing rsID codes (e.g., "rs123|rs456").
    
    Returns
    -------
    str
        Pipe-separated string with numeric parts only (e.g., "123|456").
    
    Examples
    --------
    >>> formatear_rs("rs123|rs456")
    '123|456'
    >>> formatear_rs("rs789")
    '789'
    >>> formatear_rs("123|rs456")  # Mixed format
    '123|456'
    """
    # Split the string by '|' to get each code
    codigos = cadena.split('|')
    # Remove the "rs" prefix from each code
    codigos_solo_numeros = [codigo[2:] if codigo.startswith("rs") else codigo for codigo in codigos]

    return '|'.join(codigos_solo_numeros)

def formatear_chr(string: str):
    """
    Formats chromosome identifiers to a standard format.
    
    Converts chromosome identifiers to the standard 'chr' format:
    - Converts "23" to "X"
    - Converts "24" to "Y" 
    - Adds "chr" prefix if not already present
    - Leaves existing "chr" prefixed identifiers unchanged
    
    Parameters
    ----------
    string : str
        Chromosome identifier to format.
    
    Returns
    -------
    str
        Standardized chromosome identifier with 'chr' prefix or 'X'/'Y' for sex chromosomes.
    
    Examples
    --------
    >>> formatear_chr("1")
    'chr1'
    >>> formatear_chr("23")
    'X'
    >>> formatear_chr("24")
    'Y'
    >>> formatear_chr("chr5")
    'chr5'
    """
    if string == "23":
        return "X"
    elif string == "24":
        return "Y"
    elif string.startswith("chr"):
        return string
    else:
        return "chr" + string

def normalize_variant_classification(df: pd.DataFrame) -> pd.DataFrame:
    """
    Converts values in any 'Variant Classification' column to uppercase,
    regardless of Gencode prefix, version, or capitalization.

    Examples of matched column names:
        - Gencode_43_variantClassification
        - gencode_34_variantclassification
        - variant_classification
        - Variant_Classification
        - gencode_99_VariantClassification

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.

    Returns
    -------
    pd.DataFrame
        The same DataFrame with matched columns normalized to uppercase.
        The same reference is returned to allow method chaining if desired.
    """
    # Regular expression:
    #  - ^                  : start of string
    #  - (gencode_\d+_)?    : optional prefix 'gencode_<num>_' (case insensitive)
    #  - variant[_]?classification : body of the name (allows 'variantclassification' or with '_')
    #  - $                  : end of string
    patron = re.compile(r'^(gencode_\d+_)?variant[_]?classification$', flags=re.IGNORECASE)

    # Find columns that match the pattern
    columnas_objetivo = [col for col in df.columns if patron.match(col)]

    # Convert values to uppercase for each found column
    for col in columnas_objetivo:
        # Only makes sense for object type columns (strings)
        if pd.api.types.is_string_dtype(df[col]):
            df[col] = df[col].str.upper()

    return df