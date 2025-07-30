"""
Annotation module for pyMut.

This module provides comprehensive annotation functionality for genetic variants,
including VEP (Variant Effect Predictor) annotation, COSMIC Cancer Gene Census
annotation, and OncoKB actionable mutation annotation.

Key features:
- VEP annotation for protein, gene, and variant classification
- COSMIC Cancer Gene Census integration for cancer-related gene annotation
- OncoKB integration for actionable mutation annotation
- Batch processing capabilities for large datasets
- Support for both MAF and VCF formats
"""

# Import modules without executing their module-level code immediately
from . import vep_annotate
from . import cosmic_cancer_annotate
from . import actionable_mutation

# Make specific classes and functions available at package level
from .vep_annotate import wrap_maf_vep_annotate_protein, wrap_vcf_vep_annotate_unified
from .cosmic_cancer_annotate import CancerAnnotateMixin, _maf_COSMIC_OncoKB_annotation_aux
from .actionable_mutation import ActionableMutationMixin

__all__ = [
    # Mixins for PyMutation class
    'CancerAnnotateMixin',
    'ActionableMutationMixin',
    
    # VEP annotation functions
    'wrap_maf_vep_annotate_protein',
    'wrap_vcf_vep_annotate_unified',
    
    # Utility functions for advanced usage
    '_maf_COSMIC_OncoKB_annotation_aux',
    
    # Modules for advanced usage
    'vep_annotate',
    'cosmic_cancer_annotate',
    'actionable_mutation',
]
