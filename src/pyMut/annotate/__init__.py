"""
VEP Annotation Module

This module provides functionality for annotating variants.
"""

from .vep_annotate import wrap_maf_vep_annotate_protein, wrap_vcf_vep_annotate_protein

__all__ = ['wrap_maf_vep_annotate_protein', 'wrap_vcf_vep_annotate_protein']
