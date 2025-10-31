"""
ViroForge: Synthetic virome data generator

A comprehensive mock metavirome data generator for testing and validating
virome analysis pipelines.
"""

__version__ = "0.1.0"
__author__ = "Scott Handley Lab"
__email__ = "scott.handley@wustl.edu"

# Available modules
# Phase 1: Core community & contamination (80% complete)
# from .core import community, contamination

# Phase 2: VLP enrichment & amplification bias (complete)
from . import enrichment
from . import amplification

# Future modules
# from .simulators import illumina
# from .artifacts import platform_specific
# from .utils import genome_sampler, abundance, metrics
