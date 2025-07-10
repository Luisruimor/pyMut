"""
Utilities module for pyMut.

This module contains utility functions for processing and manipulating
genetic mutation data.
"""

from .merge_vep_annotation import (
    merge_maf_with_vep_annotations,
)

__all__ = [
    'merge_maf_with_vep_annotations',
]
