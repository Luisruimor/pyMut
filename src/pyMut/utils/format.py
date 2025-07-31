import re

import pandas as pd


def format_rs(string: str) -> str:
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
    >>> format_rs("rs123|rs456")
    '123|456'
    >>> format_rs("rs789")
    '789'
    >>> format_rs("123|rs456")  # Mixed format
    '123|456'
    """
    # Split the string by '|' to get each code
    codigos = string.split('|')
    # Remove the "rs" prefix from each code
    codigos_solo_numeros = [codigo[2:] if codigo.startswith("rs") else codigo for codigo in codigos]

    return '|'.join(codigos_solo_numeros)


def format_chr(string: str):
    """
    Formats chromosome identifiers to a standard format.

    Converts chromosome identifiers to the standard 'chr' format:
    - Converts "23" to "X"
    - Converts "24" to "Y"
    - Converts "25" to "MT"
    - Converts "chr23" to "chrX"
    - Converts "chr24" to "chrY"
    - Converts "chr25" to "chrM"
    - Adds "chr" prefix if not already present
    - Leaves existing "chr" prefixed identifiers unchanged

    Parameters
    ----------
    string : str
        Chromosome identifier to format.

    Returns
    -------
    str
        Standardized chromosome identifier with 'chr' prefix or 'X'/'Y'/'MT' for special chromosomes.

    Examples
    --------
    >>> format_chr("1")
    'chr1'
    >>> format_chr("23")
    'X'
    >>> format_chr("24")
    'Y'
    >>> format_chr("25")
    'MT'
    >>> format_chr("chr23")
    'chrX'
    >>> format_chr("chr24")
    'chrY'
    >>> format_chr("chr25")
    'chrM'
    >>> format_chr("chr5")
    'chr5'
    """
    # Handle numeric chromosome identifiers
    if string == "23":
        return "X"
    elif string == "24":
        return "Y"
    elif string == "25":
        return "MT"
    # Handle chr-prefixed chromosome identifiers
    elif string == "chr23":
        return "chrX"
    elif string == "chr24":
        return "chrY"
    elif string == "chr25":
        return "chrM"
    elif string.startswith("chr"):
        return string
    else:
        return "chr" + string


def reverse_format_chr(string: str) -> str:
    """
    Removes the 'chr' prefix from chromosome identifiers.

    Converts chromosome identifiers by removing the 'chr' prefix:
    - Converts "chr1" to "1"
    - Converts "chrX" to "X"
    - Converts "chrY" to "Y"
    - Converts "chrM" to "MT"
    - Converts "chrMT" to "MT"
    - Leaves identifiers without 'chr' prefix unchanged

    Parameters
    ----------
    string : str
        Chromosome identifier to format.

    Returns
    -------
    str
        Chromosome identifier without 'chr' prefix.

    Examples
    --------
    >>> reverse_format_chr("chr1")
    '1'
    >>> reverse_format_chr("chrX")
    'X'
    >>> reverse_format_chr("chrY")
    'Y'
    >>> reverse_format_chr("chrM")
    'MT'
    >>> reverse_format_chr("chrMT")
    'MT'
    >>> reverse_format_chr("X")
    'X'
    >>> reverse_format_chr("1")
    '1'
    """
    if string.startswith("chr"):
        # Remove the "chr" prefix
        result = string[3:]
        # Handle special case: chrM -> MT
        if result == "M":
            return "MT"
        return result
    else:
        return string


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
    pattern = re.compile(r'^(gencode_\d+_)?variant[_]?classification$', flags=re.IGNORECASE)

    # Find columns that match the pattern
    variant_col = [col for col in df.columns if pattern.match(col)]

    # Convert values to uppercase for each found column
    for col in variant_col:
        # Only makes sense for object type columns (strings)
        if pd.api.types.is_string_dtype(df[col]):
            df[col] = df[col].str.upper()

    return df
