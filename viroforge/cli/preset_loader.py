#!/usr/bin/env python3
"""
Preset Loader

Load and manage ViroForge configuration presets.

Author: ViroForge Development Team
Date: 2025-11-10
"""

import yaml
from pathlib import Path
from typing import Dict, List, Optional


def get_preset_directories() -> List[Path]:
    """Get list of preset directories (built-in and user)."""
    directories = []

    # Built-in presets
    builtin_dir = Path(__file__).parent.parent / "presets"
    if builtin_dir.exists():
        directories.append(builtin_dir)

    # User presets (~/.viroforge/presets/)
    user_dir = Path.home() / ".viroforge" / "presets"
    if user_dir.exists():
        directories.append(user_dir)

    return directories


def load_preset(name: str) -> Optional[Dict]:
    """
    Load a preset by name.

    Parameters
    ----------
    name : str
        Preset name (without .yaml extension)

    Returns
    -------
    Optional[Dict]
        Preset configuration, or None if not found
    """
    for preset_dir in get_preset_directories():
        preset_file = preset_dir / f"{name}.yaml"
        if preset_file.exists():
            with open(preset_file) as f:
                return yaml.safe_load(f)

    return None


def list_all_presets() -> Dict[str, List[Dict]]:
    """
    List all available presets organized by category.

    Returns
    -------
    Dict[str, List[Dict]]
        Presets grouped by category
    """
    presets_by_category = {}
    seen_names = set()

    for preset_dir in get_preset_directories():
        if not preset_dir.exists():
            continue

        for preset_file in preset_dir.glob("*.yaml"):
            # Skip if we've already seen this preset (user overrides built-in)
            if preset_file.stem in seen_names:
                continue

            try:
                with open(preset_file) as f:
                    preset = yaml.safe_load(f)

                # Add to appropriate category
                category = preset.get('category', 'other')
                if category not in presets_by_category:
                    presets_by_category[category] = []

                preset['preset_name'] = preset_file.stem
                preset['source'] = 'user' if 'viroforge/presets' not in str(preset_dir) else 'built-in'

                presets_by_category[category].append(preset)
                seen_names.add(preset_file.stem)

            except Exception as e:
                print(f"Warning: Could not load preset {preset_file}: {e}")

    return presets_by_category


def get_preset_names() -> List[str]:
    """Get list of all preset names."""
    names = []
    for preset_dir in get_preset_directories():
        if not preset_dir.exists():
            continue
        for preset_file in preset_dir.glob("*.yaml"):
            if preset_file.stem not in names:
                names.append(preset_file.stem)
    return sorted(names)


def create_preset_from_dict(name: str, config: Dict, user_preset: bool = True) -> Path:
    """
    Create a new preset from a configuration dictionary.

    Parameters
    ----------
    name : str
        Preset name
    config : Dict
        Configuration dictionary
    user_preset : bool
        If True, save to user presets directory

    Returns
    -------
    Path
        Path to created preset file
    """
    if user_preset:
        preset_dir = Path.home() / ".viroforge" / "presets"
        preset_dir.mkdir(parents=True, exist_ok=True)
    else:
        preset_dir = Path(__file__).parent.parent / "presets"

    preset_file = preset_dir / f"{name}.yaml"

    with open(preset_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    return preset_file


def validate_preset(preset: Dict) -> tuple[bool, List[str]]:
    """
    Validate a preset configuration.

    Parameters
    ----------
    preset : Dict
        Preset configuration

    Returns
    -------
    tuple[bool, List[str]]
        (is_valid, error_messages)
    """
    errors = []

    # Required fields
    if 'name' not in preset:
        errors.append("Missing 'name' field")

    if 'parameters' not in preset:
        errors.append("Missing 'parameters' field")
    else:
        params = preset['parameters']

        # Check for required parameters
        if 'collection_id' not in params:
            errors.append("Missing 'collection_id' in parameters")

        if 'platform' not in params and 'short_platform' not in params:
            errors.append("Missing 'platform' or 'short_platform' in parameters")

        # Validate platform
        valid_platforms = ['novaseq', 'miseq', 'hiseq', 'pacbio-hifi', 'nanopore']
        platform = params.get('platform') or params.get('short_platform')
        if platform and platform not in valid_platforms:
            errors.append(f"Invalid platform: {platform}")

    return len(errors) == 0, errors
