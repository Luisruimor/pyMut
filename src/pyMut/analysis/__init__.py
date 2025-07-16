# Analysis module for pyMut
# Contains mutation burden analysis, pfam annotation, and mutational signature functionality

# Import modules without executing their module-level code immediately
from . import mutation_burden
from . import pfam_annotation
from . import mutational_signature

# Make specific functions available at package level
from .mutation_burden import calculate_tmb_analysis, log_tmb_summary
from .pfam_annotation import annotate_pfam, pfam_domains
from .mutational_signature import trinucleotideMatrix, add_trinucleotide_method_to_pymutation
