"""
Utilities module for pyMut.

This module provides comprehensive utility functions for processing, formatting,
and manipulating genetic mutation data. It includes tools for data formatting,
field manipulation, data processing, and VEP annotation integration.

Key features:
- Data formatting utilities (chromosome, rsID, variant classification)
- Field manipulation and canonicalization functions
- Data processing and extraction utilities
- VEP annotation merging capabilities
- Constants and configuration values used throughout the library
- Database utilities and data subset creation tools
"""

# Import modules without executing their module-level code immediately
from . import constants
from . import format
from . import fields
from . import data_processing
from . import merge_vep_annotation
from . import database
from . import create_subset

# Make specific functions and constants available at package level
# Constants
from .constants import (
    VARIANT_CLASSIFICATION_COLUMN,
    VARIANT_TYPE_COLUMN,
    SAMPLE_COLUMN,
    GENE_COLUMN,
    REF_COLUMN,
    ALT_COLUMN,
    FUNCOTATION_COLUMN,
    GENOME_CHANGE_COLUMN,
    DEFAULT_UNKNOWN_VALUE,
    DEFAULT_SUMMARY_FIGSIZE,
    DEFAULT_PLOT_FIGSIZE,
    DEFAULT_WATERFALL_FIGSIZE,
    DEFAULT_PLOT_TITLE,
    DEFAULT_TOP_GENES_COUNT,
    DEFAULT_ONCOPLOT_TOP_GENES,
    DEFAULT_ONCOPLOT_MAX_SAMPLES,
    DEFAULT_ONCOPLOT_FIGSIZE,
    MODE_VARIANTS,
    MODE_SAMPLES,
    VALID_PLOT_MODES,
)

# Format utilities
from .format import (
    format_rs,
    format_chr,
    reverse_format_chr,
    normalize_variant_classification,
)

# Field manipulation utilities
from .fields import (
    canonical_name,
    find_alias,
    col,
    canonicalize_columns,
)

# Data processing utilities
from .data_processing import (
    extract_variant_classification,
    extract_variant_classifications,
    extract_variant_type,
    extract_variant_types,
    extract_genome_change,
    extract_genome_changes,
    read_tsv,
)

# VEP annotation utilities
from .merge_vep_annotation import (
    merge_maf_with_vep_annotations,
)

__all__ = [
    # Constants
    'VARIANT_CLASSIFICATION_COLUMN',
    'VARIANT_TYPE_COLUMN',
    'SAMPLE_COLUMN',
    'GENE_COLUMN',
    'REF_COLUMN',
    'ALT_COLUMN',
    'FUNCOTATION_COLUMN',
    'GENOME_CHANGE_COLUMN',
    'DEFAULT_UNKNOWN_VALUE',
    'DEFAULT_SUMMARY_FIGSIZE',
    'DEFAULT_PLOT_FIGSIZE',
    'DEFAULT_WATERFALL_FIGSIZE',
    'DEFAULT_PLOT_TITLE',
    'DEFAULT_TOP_GENES_COUNT',
    'DEFAULT_ONCOPLOT_TOP_GENES',
    'DEFAULT_ONCOPLOT_MAX_SAMPLES',
    'DEFAULT_ONCOPLOT_FIGSIZE',
    'MODE_VARIANTS',
    'MODE_SAMPLES',
    'VALID_PLOT_MODES',
    
    # Format utilities
    'format_rs',
    'format_chr',
    'reverse_format_chr',
    'normalize_variant_classification',
    
    # Field manipulation utilities
    'canonical_name',
    'find_alias',
    'col',
    'canonicalize_columns',
    
    # Data processing utilities
    'extract_variant_classification',
    'extract_variant_classifications',
    'extract_variant_type',
    'extract_variant_types',
    'extract_genome_change',
    'extract_genome_changes',
    'read_tsv',
    
    # VEP annotation utilities
    'merge_maf_with_vep_annotations',
    
    # Modules for advanced usage
    'constants',
    'format',
    'fields',
    'data_processing',
    'merge_vep_annotation',
    'database',
    'create_subset',
]
