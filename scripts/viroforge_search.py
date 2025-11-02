#!/usr/bin/env python3
"""
ViroForge Advanced Genome Search

Advanced search tool with complex filtering, statistics, and export capabilities.

Usage:
    # Complex search
    viroforge-search --family Siphoviridae --length-min 30000 --gc-range 0.4-0.6 --show-stats

    # Search and export
    viroforge-search --realm Riboviria --genome-type ssRNA --export results.fasta

    # Compare two searches
    viroforge-search --compare --query1 "family:Siphoviridae" --query2 "family:Myoviridae"

Author: ViroForge Development Team
Date: 2025-11-01
