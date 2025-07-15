"""
Visualizations module for pyMut.

This module contains functions and classes for generating visualizations
of genetic mutation data.
"""
# src/pyMut/visualizations/__init__.py
from .summary import (
    create_variant_classification_plot,
    create_variant_type_plot,
    create_snv_class_plot,
    create_variant_classification_summary_plot,
    create_variants_per_sample_plot,
    create_top_mutated_genes_plot,
    create_summary_plot
)

from .oncoplot import (
    create_oncoplot_plot
)

from .mutational_signature import (
    create_mutational_signature_analysis,
    extract_mutation_matrix,
    perform_nmf,
    create_signature_profile_plot,
    create_cosine_similarity_heatmap,
    create_signature_contribution_heatmap,
    create_signature_contribution_barplot,
    create_signature_donut_plot
)

__all__ = [
    # Summary visualizations
    "create_variant_classification_plot",
    "create_variant_type_plot",
    "create_snv_class_plot",
    "create_variant_classification_summary_plot",
    "create_variants_per_sample_plot",
    "create_top_mutated_genes_plot",
    "create_summary_plot",
    # Oncoplot
    "create_oncoplot_plot",
    # Mutational signature analysis
    "create_mutational_signature_analysis",
    "extract_mutation_matrix",
    "perform_nmf",
    "create_signature_profile_plot",
    "create_cosine_similarity_heatmap",
    "create_signature_contribution_heatmap",
    "create_signature_contribution_barplot",
    "create_signature_donut_plot"
]
