"""
Module for processing mutation data.

This module contains functions for processing and transforming genetic
mutation data for subsequent visualization in summary charts.
"""

import os

import pandas as pd


def extract_variant_classification(funcotation_str):
    """
    Receives a string from the FUNCOTATION field and extracts the variant classification value,
    which is assumed to be the sixth field (index 5) of the pipe-separated string.
    
    Args:
        funcotation_str: String from the FUNCOTATION field.
        
    Returns:
        Extracted variant classification or None if it cannot be extracted.
    """
    # Ensure we're working with a string (in case it's another type)
    funcotation_str = str(funcotation_str).strip()

    # If the string includes a prefix like "FUNCOTATION=", remove it
    if funcotation_str.startswith("FUNCOTATION="):
        funcotation_str = funcotation_str[len("FUNCOTATION="):].strip()

    # Remove initial and final brackets, if present
    if funcotation_str.startswith('[') and funcotation_str.endswith(']'):
        funcotation_str = funcotation_str[1:-1]

    # Split the string by "|" and extract the sixth element (index 5)
    fields = funcotation_str.split('|')
    if len(fields) > 5:
        return fields[5].strip()
    else:
        return None


def extract_variant_classifications(data: pd.DataFrame,
                                    variant_column: str = "Variant_Classification",
                                    funcotation_column: str = "FUNCOTATION") -> pd.DataFrame:
    """
    Extract variant classification from the FUNCOTATION field if there is no specific
    column for it.
    
    Args:
        data: DataFrame with mutation data.
        variant_column: Name of the column that should contain the variant classification.
        funcotation_column: Name of the column containing the FUNCOTATION string.
        
    Returns:
        DataFrame with the variant classification column added or unchanged if it already existed.
    """
    # Create a copy to not modify the original
    result = data.copy()

    # Check if the classification column already exists
    if variant_column not in result.columns:
        # If the FUNCOTATION column exists, extract from there
        if funcotation_column in result.columns:
            print(f"Column '{variant_column}' not found. Extracting from '{funcotation_column}'.")
            # Create the column by applying the extraction function
            result[variant_column] = result[funcotation_column].apply(
                lambda x: extract_variant_classification(x) if pd.notnull(x) else None
            )

            # Check if valid data was obtained
            if result[variant_column].isna().all():
                print(f"Failed to extract variant classification data from '{funcotation_column}'.")
                result[variant_column] = "Unknown"
        else:
            print(f"Neither '{variant_column}' nor '{funcotation_column}' columns were found in the DataFrame.")
            # Create a column with unknown value
            result[variant_column] = "Unknown"

    # Fill NaN values with "Unknown"
    result[variant_column] = result[variant_column].fillna("Unknown")

    return result


def extract_variant_type(funcotation_str):
    """
    Receives a string from the FUNCOTATION field and extracts the variant type value,
    which is assumed to be the eighth field (index 7) of the pipe-separated string.
    
    Args:
        funcotation_str: String from the FUNCOTATION field.
        
    Returns:
        Extracted variant type or None if it cannot be extracted.
    """
    # Ensure we're working with a string (in case it's another type)
    funcotation_str = str(funcotation_str).strip()

    # If the string includes a prefix like "FUNCOTATION=", remove it
    if funcotation_str.startswith("FUNCOTATION="):
        funcotation_str = funcotation_str[len("FUNCOTATION="):].strip()

    # Remove initial and final brackets, if present
    if funcotation_str.startswith('[') and funcotation_str.endswith(']'):
        funcotation_str = funcotation_str[1:-1]

    # Split the string by "|" and extract the eighth element (index 7)
    fields = funcotation_str.split('|')
    if len(fields) > 7:
        return fields[7].strip()
    else:
        return None


def extract_variant_types(data: pd.DataFrame,
                          variant_column: str = "Variant_Type",
                          funcotation_column: str = "FUNCOTATION") -> pd.DataFrame:
    """
    Extract variant type from the FUNCOTATION field if there is no specific
    column for it.
    
    Args:
        data: DataFrame with mutation data.
        variant_column: Name of the column that should contain the variant type.
        funcotation_column: Name of the column containing the FUNCOTATION string.
        
    Returns:
        DataFrame with the variant type column added or unchanged if it already existed.
    """
    # Create a copy to not modify the original
    result = data.copy()

    # Check if the variant type column already exists
    if variant_column not in result.columns:
        # If the FUNCOTATION column exists, extract from there
        if funcotation_column in result.columns:
            print(f"Column '{variant_column}' not found. Extracting from '{funcotation_column}'.")
            # Create the column by applying the extraction function
            result[variant_column] = result[funcotation_column].apply(
                lambda x: extract_variant_type(x) if pd.notnull(x) else None
            )

            # Check if valid data was obtained
            if result[variant_column].isna().all():
                print(f"Failed to extract variant type data from '{funcotation_column}'.")
                result[variant_column] = "Unknown"
        else:
            print(f"Neither '{variant_column}' nor '{funcotation_column}' columns were found in the DataFrame.")
            # Create a column with unknown value
            result[variant_column] = "Unknown"

    # Fill NaN values with "Unknown"
    result[variant_column] = result[variant_column].fillna("Unknown")

    return result


def extract_genome_change(funcotation_str):
    """
    Receives a string from the FUNCOTATION field and extracts the genomic change value,
    which is assumed to be the twelfth field (index 11) of the pipe-separated string.
    
    Args:
        funcotation_str: String from the FUNCOTATION field.
        
    Returns:
        Extracted genomic change or None if it cannot be extracted.
    """
    # Convert to string and remove whitespace
    funcotation_str = str(funcotation_str).strip()

    # Remove possible "FUNCOTATION=" prefix
    if funcotation_str.startswith("FUNCOTATION="):
        funcotation_str = funcotation_str[len("FUNCOTATION="):].strip()

    # Remove opening and closing brackets if present
    if funcotation_str.startswith('[') and funcotation_str.endswith(']'):
        funcotation_str = funcotation_str[1:-1]

    # Split the string by "|" and extract the corresponding field (index 11)
    fields = funcotation_str.split('|')
    if len(fields) > 11:
        return fields[11].strip()
    else:
        return None


def extract_genome_changes(data: pd.DataFrame,
                           genome_change_column: str = "Genome_Change",
                           funcotation_column: str = "FUNCOTATION") -> pd.DataFrame:
    """
    Extract genomic change from the FUNCOTATION field if there is no specific
    column for it.
    
    Args:
        data: DataFrame with mutation data.
        genome_change_column: Name of the column that should contain the genomic change.
        funcotation_column: Name of the column containing the FUNCOTATION string.
        
    Returns:
        DataFrame with the genomic change column added or unchanged if it already existed.
    """
    # Create a copy to not modify the original
    result = data.copy()

    # Check if the genomic change column already exists
    if genome_change_column not in result.columns:
        # If the FUNCOTATION column exists, extract from there
        if funcotation_column in result.columns:
            print(f"Column '{genome_change_column}' not found. Extracting from '{funcotation_column}'.")
            # Create the column by applying the extraction function
            result[genome_change_column] = result[funcotation_column].apply(
                lambda x: extract_genome_change(x) if pd.notnull(x) else None
            )

            # Check if valid data was obtained
            if result[genome_change_column].isna().all():
                print(f"Failed to extract genome change data from '{funcotation_column}'.")
                result[genome_change_column] = "Unknown"
        else:
            print(f"Neither '{genome_change_column}' nor '{funcotation_column}' columns were found in the DataFrame.")
            # Create a column with unknown value
            result[genome_change_column] = "Unknown"

    # Fill NaN values with "Unknown"
    result[genome_change_column] = result[genome_change_column].fillna("Unknown")

    return result


def read_tsv(file_path: str) -> pd.DataFrame:
    """
    Read a TSV file and return a DataFrame with mutation data.
    
    Args:
        file_path: Path to the TSV file.
        
    Returns:
        DataFrame with mutation data.
    
    Raises:
        FileNotFoundError: If the file does not exist.
        pd.errors.EmptyDataError: If the file is empty.
        pd.errors.ParserError: If there are problems parsing the file.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File '{file_path}' does not exist.")

    try:
        # First try to read without comments
        data = pd.read_csv(file_path, sep='\t')
    except (pd.errors.ParserError, pd.errors.EmptyDataError):
        # If it fails due to a parsing error, try with the comment parameter
        try:
            data = pd.read_csv(file_path, sep='\t', comment='#', engine='python')
        except Exception as err:
            raise ValueError(f"Could not read file: {str(err)}")
    except Exception as e:
        raise ValueError(f"Error reading file: {str(e)}")

    return data
