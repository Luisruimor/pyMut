"""
Constantes utilizadas en la librería pyMut.

Este módulo contiene las constantes que se utilizan en varios lugares de la librería,
para mantener consistencia y facilitar el mantenimiento.
"""

# Nombres de columnas comunes
VARIANT_CLASSIFICATION_COLUMN = "Variant_Classification"
VARIANT_TYPE_COLUMN = "Variant_Type"
SAMPLE_COLUMN = "Tumor_Sample_Barcode"
GENE_COLUMN = "Hugo_Symbol"
REF_COLUMN = "REF"
ALT_COLUMN = "ALT"
FUNCOTATION_COLUMN = "FUNCOTATION"
GENOME_CHANGE_COLUMN = "Genome_Change"

# Valores por defecto
DEFAULT_UNKNOWN_VALUE = "Unknown"

# Parámetros de visualización
DEFAULT_SUMMARY_FIGSIZE = (16, 12)
DEFAULT_PLOT_FIGSIZE = (10, 6)
DEFAULT_PLOT_TITLE = "Resumen de Mutaciones"
DEFAULT_TOP_GENES_COUNT = 10

# Modos de visualización
MODE_VARIANTS = "variants"
MODE_SAMPLES = "samples"
VALID_PLOT_MODES = [MODE_VARIANTS, MODE_SAMPLES]

# Formatos de datos
FORMAT_PIPE_SEPARATED = "pipe_separated"
FORMAT_SLASH_SEPARATED = "slash_separated"
FORMAT_OTHER = "other"

# Patrones para detectar columnas de muestra
TCGA_SAMPLE_PREFIX = "TCGA-"
SAMPLE_ID_MIN_HYPHENS = 2