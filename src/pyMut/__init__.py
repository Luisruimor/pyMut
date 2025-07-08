from .core import PyMutation
from .input import read_maf, read_vcf
from .version import __version__

# Import modules that add methods to PyMutation class
from . import mutation_burden
from . import output

from .pfam_annotation import (
    connect_db, VariantAnnotator, annotate_pfam, pfam_domains, build_embedded_db
)

__all__ = [
    'PyMutation', 'read_maf', 'read_vcf', 
    'connect_db', 'VariantAnnotator', 'annotate_pfam', 'pfam_domains', 'build_embedded_db'
]
