"""
Analysis module for pyMut.

Contains mutation burden analysis, pfam annotation, mutational signature, 
and SMG detection functionality.
"""

# Import modules without executing their module-level code immediately
from . import mutation_burden
from . import pfam_annotation
from . import mutational_signature
from . import smg_detection

# Make specific classes and functions available at package level
from .mutation_burden import MutationBurdenMixin, log_tmb_summary
from .pfam_annotation import PfamAnnotationMixin
from .mutational_signature import MutationalSignatureMixin
from .smg_detection import SmgDetectionMixin

__all__ = [
    # Mixins
    'MutationBurdenMixin',
    'PfamAnnotationMixin', 
    'MutationalSignatureMixin',
    'SmgDetectionMixin',
    
    # Utility functions
    'log_tmb_summary',
    
    # Modules for advanced usage
    'mutation_burden',
    'pfam_annotation',
    'mutational_signature',
    'smg_detection',
]
