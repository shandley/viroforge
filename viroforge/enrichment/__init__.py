"""
ViroForge Enrichment Module

VLP (Virus-Like Particle) enrichment modeling for realistic virome simulation.

Main classes:
- VLPEnrichment: Main enrichment model
- VLPProtocol: Pre-defined protocol configurations
- VirionSizeEstimator: Estimate virion physical properties
- FiltrationCurve: Size-based retention modeling

Usage:
    from viroforge.enrichment import VLPEnrichment, VLPProtocol

    # Create VLP enrichment with tangential flow filtration
    vlp = VLPEnrichment(protocol=VLPProtocol.tangential_flow_standard())

    # Apply to genome abundances
    enriched_abundances, stats = vlp.apply_enrichment(genomes, abundances)
"""

from .vlp import (
    VLPEnrichment,
    VLPProtocol,
    VLPProtocolConfig,
    VirionSizeEstimator,
    VirionSizeEstimate,
    FiltrationCurve
)

__all__ = [
    'VLPEnrichment',
    'VLPProtocol',
    'VLPProtocolConfig',
    'VirionSizeEstimator',
    'VirionSizeEstimate',
    'FiltrationCurve'
]
