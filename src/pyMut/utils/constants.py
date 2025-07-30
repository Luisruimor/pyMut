"""
Constants used in the pyMut library.

This module contains the constants used throughout the library
to maintain consistency and facilitate maintenance.
"""

# Common column names
VARIANT_CLASSIFICATION_COLUMN = "Variant_Classification"
VARIANT_TYPE_COLUMN = "Variant_Type"
SAMPLE_COLUMN = "Tumor_Sample_Barcode"
GENE_COLUMN = "Hugo_Symbol"
REF_COLUMN = "REF"
ALT_COLUMN = "ALT"
FUNCOTATION_COLUMN = "FUNCOTATION"
GENOME_CHANGE_COLUMN = "Genome_Change"

# Default values
DEFAULT_UNKNOWN_VALUE = "Unknown"

# Visualization parameters
DEFAULT_SUMMARY_FIGSIZE = (14, 8)
DEFAULT_PLOT_FIGSIZE = (12, 6)
DEFAULT_WATERFALL_FIGSIZE = (16, 10)
DEFAULT_PLOT_TITLE = "Mutation Summary"
DEFAULT_TOP_GENES_COUNT = 10

# Oncoplot specific parameters
DEFAULT_ONCOPLOT_TOP_GENES = 10
DEFAULT_ONCOPLOT_MAX_SAMPLES = 180
DEFAULT_ONCOPLOT_FIGSIZE = (16, 10)

# Visualization modes
MODE_VARIANTS = "variants"
MODE_SAMPLES = "samples"
VALID_PLOT_MODES = [MODE_VARIANTS, MODE_SAMPLES]

# Data formats
FORMAT_PIPE_SEPARATED = "pipe_separated"
FORMAT_SLASH_SEPARATED = "slash_separated"
FORMAT_OTHER = "other"

# Patterns for detecting sample columns
TCGA_SAMPLE_PREFIX = "TCGA-"
SAMPLE_ID_MIN_HYPHENS = 2
