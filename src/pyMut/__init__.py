"""
pyMut: A comprehensive Python library for mutation analysis and visualization.

This library provides tools for:
- Reading and writing mutation data (MAF, VCF formats)
- Mutation burden analysis and TMB calculation
- Mutational signature analysis
- Pfam domain annotation
- VEP annotation integration
- Filtering and genomic range operations
- Comprehensive visualization tools
- Statistical analysis and SMG detection
"""

# Core classes and main functionality
from .core import PyMutation, MutationMetadata
from .input import read_maf, read_vcf
from .combination import combine_pymutations
from .version import __version__, get_version

# Import modules that add methods to PyMutation class via mixins
# This ensures all mixin methods are available on PyMutation instances
from .analysis import mutation_burden
from .analysis import pfam_annotation
from .analysis import mutational_signature
from .analysis import smg_detection
from . import output
from . import filters

# Key utility functions for annotation
from .annotate import wrap_maf_vep_annotate_protein, wrap_vcf_vep_annotate_unified

# Important utility functions
from .utils import merge_maf_with_vep_annotations

# Key analysis functions that can be used independently
from .analysis import log_tmb_summary

# Configure logging for the library
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

# Public API - everything users should be able to import directly
__all__ = [
    # Core functionality
    'PyMutation',
    'MutationMetadata',
    
    # Data I/O
    'read_maf',
    'read_vcf',
    'combine_pymutations',
    
    # Annotation functions
    'wrap_maf_vep_annotate_protein',
    'wrap_vcf_vep_annotate_unified',
    
    # Utility functions
    'merge_maf_with_vep_annotations',
    'log_tmb_summary',
    
    # Version info
    '__version__',
    'get_version',
]

# Library metadata
__author__ = "Luis Ruiz Moreno"
__email__ = "luisruimore@gmail.com"
__description__ = "A comprehensive Python library for mutation analysis and visualization"
__url__ = "https://github.com/luisruimore/pyMut"
