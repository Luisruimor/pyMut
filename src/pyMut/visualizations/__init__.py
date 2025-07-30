"""
Visualizations module for pyMut.

This module contains functions and classes for generating visualizations
of genetic mutation data.
"""

from .summary import (
    _create_variant_classification_plot,
    _create_variant_type_plot,
    _create_snv_class_plot,
    _create_variants_per_sample_plot,
    _create_variant_classification_summary_plot,
    _create_top_mutated_genes_plot,
    _create_summary_plot,
)
from .oncoplot import _create_oncoplot_plot

# Note: These are internal functions used by PyMutation methods
# They are not exposed in __all__ as they are not meant for direct use
# All visualization functionality is accessed through PyMutation instance methods
__all__ = [
    # Currently no public functions - all visualization is done through PyMutation methods
    # This may change in future versions if standalone visualization functions are added
]
