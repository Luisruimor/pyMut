# Analysis module for pyMut
# Contains mutation burden analysis and pfam annotation functionality

# Import modules without executing their module-level code immediately
from . import mutation_burden
from . import pfam_annotation

# Make specific functions available at package level
from .mutation_burden import calculate_tmb_analysis, log_tmb_summary
from .pfam_annotation import annotate_pfam, pfam_domains
