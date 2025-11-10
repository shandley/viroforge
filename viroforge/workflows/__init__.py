"""
ViroForge Workflows Module

This module provides workflow classes for different virome sequencing preparations:
- RNAViromeWorkflow: RNA virus sequencing with RT and rRNA depletion
- Additional workflows planned for future phases
"""

from .rna_virome import (
    RNAViromeWorkflow,
    ReverseTranscription,
    RiboDepletion,
    RNADegradation
)

__all__ = [
    'RNAViromeWorkflow',
    'ReverseTranscription',
    'RiboDepletion',
    'RNADegradation'
]
