from .core import PyMutation
from .input import read_maf, read_vcf
from .version import __version__

# Import modules that add methods to PyMutation class
from .analysis import mutation_burden
from . import output
from .analysis import pfam_annotation
from .analysis import mutational_signature

# Methods are now automatically available through Mixin architecture


__all__ = [
    'PyMutation', 'read_maf', 'read_vcf',
]
