import json
import logging
from pathlib import Path
from typing import Union, List, Dict, Optional, Tuple

import pandas as pd

from ..utils.fields import FIELDS

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if not logger.handlers:
    logging.basicConfig(
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        level=logging.INFO,
    )

_rna_cancer_consensus_cache: Optional[Dict[str, Dict[str, float]]] = None


def _load_rna_cancer_consensus(file_path: Union[str, Path] = None) -> Dict[str, Dict[str, float]]:
    """
    Auxiliary method to load RNA cancer consensus data with caching.

    This method loads the rna_cancer_consensus.json file and caches the result
    to avoid reloading the same data multiple times.
    """
    global _rna_cancer_consensus_cache

    if _rna_cancer_consensus_cache is not None:
        return _rna_cancer_consensus_cache

    if file_path is None:
        current_dir = Path(__file__).parent
        default_path = current_dir.parent / "data" / "rna_cancer_consensus.json"
        file_path = default_path

    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"RNA cancer consensus file not found: {file_path}")

    try:
        with open(file_path, 'r') as f:
            data = json.load(f)

        _rna_cancer_consensus_cache = data
        return data

    except json.JSONDecodeError as e:
        raise json.JSONDecodeError(f"Error parsing JSON file {file_path}: {e}", e.doc, e.pos)


def get_gene_symbol(row: pd.Series) -> Optional[str]:
    """
    Get gene symbol from a row by iterating through Hugo_Symbol field list from fields.py.
    
    Returns the first non-empty value that is not empty, dot, or dash.
    If the value contains "&" (like "GENE1&GENE2"), returns only the first gene.
    
    Parameters
    ----------
    row : pd.Series
        A row from mutation data
        
    Returns
    -------
    Optional[str]
        The first valid gene symbol found, or None if no valid symbol is found
    """

    hugo_symbol_fields = FIELDS.get("Hugo_Symbol", [])

    for field in hugo_symbol_fields:
        if field in row.index and pd.notna(row[field]):
            value = str(row[field]).strip()
            # Skip empty values, dots, and dashes
            if value and value not in ["", ".", "-"]:
                # If contains "&", return only the first gene
                if "&" in value:
                    value = value.split("&")[0].strip()
                if value and value not in ["", ".", "-"]:
                    return value

    return None


def tissue_expression(data: Union[str, pd.Series], tissue: List[Union[str, float]]) -> bool:
    """
    Check if a gene is sufficiently expressed in a specific tissue/cancer based on the provided threshold.

    This auxiliary function checks if the gene where a variant falls is sufficiently "turned on" (expressed)
    in a specific tissue/cancer according to the provided threshold. The RNA cancer consensus data is loaded
    automatically with caching to avoid reloading the same data multiple times.

    Parameters
    ----------
    data : Union[str, pd.Series]
        Can be either:
        - A gene symbol (Hugo_Symbol) as a string
        - A row from mutation data that contains Hugo_Symbol or its synonyms
    tissue : List[Union[str, float]]
        A list with two elements: [tissue_code, threshold]
        Example: ['BLCA', 5] where 'BLCA' is Bladder Cancer and 5 is the expression threshold

    Returns
    -------
    bool
        True if the gene is sufficiently expressed (above threshold) in the specified tissue,
        False otherwise

    Raises
    ------
    ValueError
        If tissue parameter doesn't have exactly 2 elements or if tissue code is not found
    KeyError
        If gene symbol is not found in the expression data
    TypeError
        If data parameter is not a string or pandas Series
    FileNotFoundError
        If the RNA cancer consensus file is not found
    json.JSONDecodeError
        If the RNA cancer consensus file is malformed

    Examples
    --------
    >>> from pyMut.filters.tissue_expression import tissue_expression
    >>> 
    >>> # Using gene symbol directly (data loaded automatically with caching)
    >>> tissue_expression("TSPAN6", ["BLCA", 5])
    True

    >>> # Using a row from mutation data
    >>> import pandas as pd
    >>> row = pd.Series({'Hugo_Symbol': 'TSPAN6', 'Chromosome': 'X'})
    >>> tissue_expression(row, ["BRCA", 10])
    False

    Note
    ----
    The RNA cancer consensus data is loaded automatically from rna_cancer_consensus.json
    and cached for subsequent calls to avoid reloading the same data multiple times.
    """
    if not isinstance(tissue, list) or len(tissue) != 2:
        raise ValueError("tissue parameter must be a list with exactly 2 elements: [tissue_code, threshold]")

    tissue_code, threshold = tissue

    if not isinstance(tissue_code, str):
        raise ValueError("tissue_code must be a string")

    if not isinstance(threshold, (int, float)):
        raise ValueError("threshold must be a number")

    all_dict = _load_rna_cancer_consensus()
    gene_symbol = None

    if isinstance(data, str):
        gene_symbol = data

    elif isinstance(data, pd.Series):
        gene_symbol = get_gene_symbol(data)

        if gene_symbol is None:
            raise KeyError(
                "Could not find gene symbol in the provided data row. Expected columns: Hugo_Symbol, Gene_Symbol, Gene, etc.")

    else:
        raise TypeError("data parameter must be either a string (gene symbol) or a pandas Series (data row)")

    if gene_symbol not in all_dict:
        return False

    gene_expression_data = all_dict[gene_symbol]

    if tissue_code not in gene_expression_data:
        return False

    expression_value = gene_expression_data[tissue_code]
    is_expressed = expression_value >= threshold

    return is_expressed


class TissueExpressionMixin:
    """
    Mixin class providing tissue expression filtering functionality for PyMutation objects.
    
    This mixin adds the ability to filter by tissue expression levels,
    following the same architectural pattern as other mixins in the project.
    """

    def filter_by_tissue_expression(self, tissues: List[Tuple[str, float]], keep_expressed: bool = True) -> 'PyMutation':
        """
        Filter PyMutation data based on tissue expression for one or more tissues.

        This method filters self.data based on whether genes are sufficiently expressed
        in the specified tissues according to their respective thresholds. It can filter
        for genes that are expressed (default) or not expressed in any of the specified tissues.

        Parameters
        ----------
        tissues : List[Tuple[str, float]]
            List of tuples where each tuple contains:
            - tissue_code (str): TCGA tissue/cancer code (e.g., 'BLCA', 'BRCA', 'LUAD')
            - threshold (float): Expression threshold for that tissue
            Example: [('BLCA', 5), ('BRCA', 3), ('LUAD', 4)]
        keep_expressed : bool, default True
            If True, keeps rows where genes are expressed in at least one of the specified tissues.
            If False, keeps rows where genes are NOT expressed in any of the specified tissues.

        Returns
        -------
        PyMutation
            A new PyMutation object with filtered data

        Raises
        ------
        ValueError
            If tissues list is empty or contains invalid tissue specifications
        KeyError
            If required gene symbol columns are not found in the data

        Examples
        --------
        >>> # Filter for genes expressed in bladder cancer (threshold 5) or breast cancer (threshold 3)
        >>> filtered_mut = py_mut.filter_by_tissue_expression([('BLCA', 5), ('BRCA', 3)])

        >>> # Filter for genes NOT expressed in lung adenocarcinoma (threshold 4)
        >>> not_expressed_mut = py_mut.filter_by_tissue_expression([('LUAD', 4)], keep_expressed=False)

        >>> # Filter for genes expressed in multiple tissues with different thresholds
        >>> multi_tissue_mut = py_mut.filter_by_tissue_expression([
        ...     ('BLCA', 5), ('BRCA', 3), ('LUAD', 4), ('COAD', 6)
        ... ])
        """
        if not isinstance(tissues, list) or len(tissues) == 0:
            raise ValueError("tissues parameter must be a non-empty list of tuples")

        for i, tissue_spec in enumerate(tissues):
            if not isinstance(tissue_spec, tuple) or len(tissue_spec) != 2:
                raise ValueError(f"Each tissue specification must be a tuple with 2 elements (tissue_code, threshold). "
                                 f"Invalid specification at index {i}: {tissue_spec}")

            tissue_code, threshold = tissue_spec
            if not isinstance(tissue_code, str):
                raise ValueError(f"Tissue code must be a string. Invalid at index {i}: {tissue_code}")

            if not isinstance(threshold, (int, float)):
                raise ValueError(f"Threshold must be a number. Invalid at index {i}: {threshold}")

        results_data = []
        filtered_data = self.data.copy()
        expression_results = []

        for idx, row in filtered_data.iterrows():
            is_expressed_in_any_tissue = False

            gene_symbol = get_gene_symbol(row)

            tissue_results = {}
            for tissue_code, threshold in tissues:
                try:
                    is_expressed = tissue_expression(row, [tissue_code, threshold])
                    tissue_results[f"{tissue_code}_expressed"] = is_expressed
                    tissue_results[f"{tissue_code}_threshold"] = threshold
                    if is_expressed:
                        is_expressed_in_any_tissue = True
                        # Don't break here - we want to check all tissues for the results dataframe
                except (KeyError, ValueError):
                    tissue_results[f"{tissue_code}_expressed"] = False
                    tissue_results[f"{tissue_code}_threshold"] = threshold
                    continue

            expression_results.append(is_expressed_in_any_tissue)

            result_row = {
                'Index': idx,
                'Gene_Symbol': gene_symbol,
                'Expressed_in_Any_Tissue': is_expressed_in_any_tissue
            }
            result_row.update(tissue_results)
            results_data.append(result_row)

        results_df = pd.DataFrame(results_data)

        if keep_expressed:
            mask = pd.Series(expression_results, index=filtered_data.index)
            filtered_data = filtered_data[mask]
        else:
            mask = pd.Series([not result for result in expression_results], index=filtered_data.index)
            filtered_data = filtered_data[mask]

        from ..core import PyMutation

        new_metadata = None
        if hasattr(self, 'metadata') and self.metadata is not None:
            import copy
            new_metadata = copy.deepcopy(self.metadata)
            tissue_filter_desc = f"tissue_expression_filter({tissues}, keep_expressed={keep_expressed})"
            if hasattr(new_metadata, 'filters'):
                new_metadata.filters.append(tissue_filter_desc)
            else:
                new_metadata.filters = [tissue_filter_desc]

        new_samples = getattr(self, 'samples', [])
        new_pymutation = PyMutation(filtered_data, metadata=new_metadata, samples=new_samples)
        new_pymutation.tissue_expression_results = results_df

        return new_pymutation
