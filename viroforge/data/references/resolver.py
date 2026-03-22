"""
Reference sequence resolver for ViroForge contamination modeling.

Locates reference FASTA files using a priority chain:
1. User-supplied path (explicit argument)
2. Environment variable (VIROFORGE_*)
3. Bundled references (shipped with package)
4. None (caller falls back to synthetic generation)
"""

import logging
import os
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# Directory containing bundled reference FASTA files
_REFERENCES_DIR = Path(__file__).parent


def _resolve(
    user_path: Optional[Path],
    env_var: str,
    bundled_name: str,
    label: str,
) -> Optional[Path]:
    """Generic resolver with priority chain.

    Args:
        user_path: Explicit path provided by user (highest priority).
        env_var: Environment variable name to check.
        bundled_name: Filename of bundled reference in this package.
        label: Human-readable label for log messages.

    Returns:
        Path to reference FASTA, or None if nothing found.
    """
    # Priority 1: User-supplied path
    if user_path is not None:
        p = Path(user_path)
        if p.exists():
            logger.info(f"Using user-supplied {label}: {p}")
            return p
        else:
            logger.warning(f"User-supplied {label} not found: {p}")

    # Priority 2: Environment variable
    env_val = os.environ.get(env_var)
    if env_val:
        p = Path(env_val)
        if p.exists():
            logger.info(f"Using {label} from {env_var}: {p}")
            return p
        else:
            logger.warning(f"{env_var} set but file not found: {p}")

    # Priority 3: Bundled reference
    bundled = _REFERENCES_DIR / bundled_name
    if bundled.exists():
        logger.debug(f"Using bundled {label}: {bundled}")
        return bundled

    # Priority 4: No reference available
    logger.info(
        f"No {label} reference found; contamination will use synthetic sequences"
    )
    return None


def get_phix_path(user_path: Optional[Path] = None) -> Optional[Path]:
    """Locate PhiX174 reference genome."""
    return _resolve(user_path, "VIROFORGE_PHIX_GENOME", "phix174.fasta", "PhiX174")


def get_rrna_path(user_path: Optional[Path] = None) -> Optional[Path]:
    """Locate rRNA reference database."""
    return _resolve(
        user_path, "VIROFORGE_RRNA_DB", "rrna_representatives.fasta", "rRNA database"
    )


def get_host_fragments_path(user_path: Optional[Path] = None) -> Optional[Path]:
    """Locate host genome fragments for contamination modeling.

    For full host genomes (e.g., GRCh38), use get_host_genome_path() or
    set VIROFORGE_HOST_GENOME.
    """
    return _resolve(
        user_path,
        "VIROFORGE_HOST_FRAGMENTS",
        "host_fragments.fasta",
        "host fragments",
    )


def get_host_genome_path(user_path: Optional[Path] = None) -> Optional[Path]:
    """Locate full host genome (user-supplied only, not bundled).

    Full genomes are too large to bundle. Users can provide a path to
    GRCh38, T2T-CHM13, or other host genome FASTA.
    """
    if user_path is not None:
        p = Path(user_path)
        if p.exists():
            logger.info(f"Using user-supplied host genome: {p}")
            return p
        else:
            logger.warning(f"User-supplied host genome not found: {p}")

    env_val = os.environ.get("VIROFORGE_HOST_GENOME")
    if env_val:
        p = Path(env_val)
        if p.exists():
            logger.info(f"Using host genome from VIROFORGE_HOST_GENOME: {p}")
            return p
        else:
            logger.warning(f"VIROFORGE_HOST_GENOME set but file not found: {p}")

    return None


def get_adapter_path(user_path: Optional[Path] = None) -> Optional[Path]:
    """Locate Illumina adapter sequences."""
    return _resolve(
        user_path, "VIROFORGE_ADAPTERS", "adapters.fasta", "adapter sequences"
    )


def has_bundled_references() -> dict[str, bool]:
    """Check which bundled reference files are available."""
    return {
        "phix174": (_REFERENCES_DIR / "phix174.fasta").exists(),
        "rrna": (_REFERENCES_DIR / "rrna_representatives.fasta").exists(),
        "host_fragments": (_REFERENCES_DIR / "host_fragments.fasta").exists(),
        "adapters": (_REFERENCES_DIR / "adapters.fasta").exists(),
    }
