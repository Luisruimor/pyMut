"""
Combination module for PyMutation objects.

This module provides functionality to combine two PyMutation objects.
"""

import logging

import pandas as pd

from .core import PyMutation, MutationMetadata

# Logger configuration
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if not logger.handlers:
    logging.basicConfig(format="%(asctime)s | %(levelname)s | %(name)s | %(message)s", level=logging.INFO, )


def combine_pymutations(pymut1: PyMutation, pymut2: PyMutation) -> PyMutation:
    """
    Combine two PyMutation instances into a single PyMutation instance.

    This function combines two PyMutation instances, ensuring:
    - Both instances have the same assembly
    - Samples are combined without duplicates
    - No duplicate variants (identified by CHROM, POS, REF, ALT)
    - Annotations are combined properly (same column name = same column, different name = new column)

    Parameters
    ----------
    pymut1 : PyMutation
        First PyMutation instance
    pymut2 : PyMutation
        Second PyMutation instance

    Returns
    -------
    PyMutation
        A new PyMutation instance that combines the two input instances

    Raises
    ------
    ValueError
        If the two instances have different assemblies
    """
    logger.info(
        f"Starting combination of PyMutation instances: {pymut1.metadata.file_path} and {pymut2.metadata.file_path}")
    if pymut1.metadata.assembly != pymut2.metadata.assembly:
        logger.error(
            f"Cannot combine PyMutation instances with different assemblies: {pymut1.metadata.assembly} and {pymut2.metadata.assembly}")
        raise ValueError(f"Cannot combine PyMutation instances with different assemblies: "
                         f"{pymut1.metadata.assembly} and {pymut2.metadata.assembly}")
    logger.info(f"Assembly check passed: both instances have assembly {pymut1.metadata.assembly}")

    combined_samples = list(set(pymut1.samples + pymut2.samples))
    logger.info(
        f"Combined samples: {len(pymut1.samples)} from first instance + {len(pymut2.samples)} from second instance = {len(combined_samples)} unique samples")

    # Identify variant columns for deduplication
    variant_id_columns = ['CHROM', 'POS', 'REF', 'ALT']

    # Create a copy of the DataFrames to avoid modifying the originals
    df1 = pymut1.data.copy()
    df2 = pymut2.data.copy()

    df1['_variant_id'] = df1.apply(lambda row: f"{row['CHROM']}_{row['POS']}_{row['REF']}_{row['ALT']}", axis=1)
    df2['_variant_id'] = df2.apply(lambda row: f"{row['CHROM']}_{row['POS']}_{row['REF']}_{row['ALT']}", axis=1)
    logger.info(
        f"Created unique variant identifiers for {len(df1)} variants in first instance and {len(df2)} variants in second instance")

    # Get the set of columns from each DataFrame
    cols1 = set(df1.columns)
    cols2 = set(df2.columns)

    common_cols = cols1.intersection(cols2)
    unique_cols1 = cols1 - cols2
    unique_cols2 = cols2 - cols1
    logger.info(
        f"Column analysis: {len(common_cols)} common columns, {len(unique_cols1)} unique to first instance, {len(unique_cols2)} unique to second instance")

    # Ensure variant ID columns are in common columns
    for col in variant_id_columns:
        if col in common_cols:
            continue
        if col in unique_cols1:
            unique_cols1.remove(col)
            common_cols.add(col)
        elif col in unique_cols2:
            unique_cols2.remove(col)
            common_cols.add(col)

    # Add _variant_id to common columns
    common_cols.add('_variant_id')

    all_variant_ids = set(df1['_variant_id']).union(set(df2['_variant_id']))
    logger.info(f"Found {len(all_variant_ids)} unique variants across both instances")

    # Initialize the combined DataFrame with all unique variant IDs
    combined_df = pd.DataFrame({'_variant_id': list(all_variant_ids)})

    # Create a dictionary to store all column data
    column_data = {}

    for col in common_cols:
        if col == '_variant_id':
            continue

        # Create a mapping from variant ID to value for each DataFrame
        map1 = dict(zip(df1['_variant_id'], df1[col]))
        map2 = dict(zip(df2['_variant_id'], df2[col]))

        # Combine the values, preferring non-null values
        combined_values = []
        for variant_id in combined_df['_variant_id']:
            val1 = map1.get(variant_id)
            val2 = map2.get(variant_id)

            if pd.isna(val1) and pd.isna(val2):
                combined_values.append(None)
            elif pd.isna(val1):
                combined_values.append(val2)
            elif pd.isna(val2):
                combined_values.append(val1)
            else:
                # Both values are non-null, prefer the first one
                combined_values.append(val1)

        column_data[col] = combined_values

    for col in unique_cols1:
        map1 = dict(zip(df1['_variant_id'], df1[col]))
        column_data[col] = combined_df['_variant_id'].map(map1).tolist()

    for col in unique_cols2:
        map2 = dict(zip(df2['_variant_id'], df2[col]))
        column_data[col] = combined_df['_variant_id'].map(map2).tolist()

    # Create a DataFrame from all column data and concatenate with the variant_id DataFrame
    columns_df = pd.DataFrame(column_data, index=combined_df.index)
    combined_df = pd.concat([combined_df, columns_df], axis=1)
    logger.info(
        f"Processed {len(common_cols)} common columns, {len(unique_cols1)} columns unique to first instance, and {len(unique_cols2)} columns unique to second instance")

    # Remove the temporary variant ID column
    combined_df = combined_df.drop(columns=['_variant_id'])

    std_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']

    # Ensure all standard columns exist, add them if they don't
    for col in std_columns:
        if col not in combined_df.columns:
            combined_df[col] = None

    # Get all columns that are not standard columns and not sample columns
    annotation_columns = [col for col in combined_df.columns if col not in std_columns and col not in combined_samples]

    # Reorder the columns
    ordered_columns = std_columns + combined_samples + annotation_columns

    # Filter to only include columns that actually exist in the DataFrame
    ordered_columns = [col for col in ordered_columns if col in combined_df.columns]

    # Reorder the DataFrame
    combined_df = combined_df[ordered_columns]
    logger.info(
        f"Reordered columns: {len(std_columns)} standard columns + {len(combined_samples)} sample columns + {len(annotation_columns)} annotation columns")

    combined_metadata = MutationMetadata(
        source_format=f"COMBINED {pymut1.metadata.source_format} + {pymut2.metadata.source_format}",
        file_path=f"{pymut1.metadata.file_path} + {pymut2.metadata.file_path}", filters=["."],
        assembly=pymut1.metadata.assembly,
        notes=f"Combined from {pymut1.metadata.file_path} and {pymut2.metadata.file_path}")
    logger.info(f"Created new metadata with assembly {combined_metadata.assembly}")

    logger.info(
        f"Combination complete: created PyMutation instance with {len(combined_df)} variants and {len(combined_df.columns)} columns")
    return PyMutation(combined_df, combined_metadata, combined_samples)
