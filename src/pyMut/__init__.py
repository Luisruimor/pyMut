from .core import PyMutation
from .input import read_maf, read_vcf
from .version import __version__

# Import modules that add methods to PyMutation class
from .analysis import mutation_burden
from . import output
from .analysis import pfam_annotation


__all__ = [
    'PyMutation', 'read_maf', 'read_vcf',
]
