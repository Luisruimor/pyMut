from .core import PyMutation
from .input import read_maf, read_vcf
from .version import __version__

# Import modules that add methods to PyMutation class
from .analysis import mutation_burden
from . import output
from .analysis import pfam_annotation
from .analysis import mutational_signature

# Add methods to PyMutation class
mutational_signature.add_trinucleotide_method_to_pymutation()


__all__ = [
    'PyMutation', 'read_maf', 'read_vcf',
]
