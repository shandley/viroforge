"""
ViroForge: Synthetic virome data generator
"""

import re
from pathlib import Path

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


def _read_version() -> str:
    """Single source of truth: read __version__ from viroforge/__init__.py.

    Parsed via regex rather than imported, so building the package does not
    require importing viroforge (and its runtime dependencies).
    """
    init = Path(__file__).parent / "viroforge" / "__init__.py"
    match = re.search(r'^__version__\s*=\s*["\']([^"\']+)["\']', init.read_text(encoding="utf-8"), re.M)
    if not match:
        raise RuntimeError("Unable to find __version__ in viroforge/__init__.py")
    return match.group(1)


setup(
    name="viroforge",
    version=_read_version(),
    author="Scott Handley Lab",
    author_email="scott.handley@wustl.edu",
    description="A comprehensive mock metavirome data generator for benchmarking and validation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/shandley/viroforge",
    project_urls={
        "Bug Tracker": "https://github.com/shandley/viroforge/issues",
        "Documentation": "https://github.com/shandley/viroforge/docs",
        "Source Code": "https://github.com/shandley/viroforge",
    },
    packages=find_packages(),
    package_data={
        "viroforge": ["data/references/*.fasta"],
    },
    include_package_data=True,
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "biopython>=1.79",
        "pyyaml>=5.4.0",
        "scipy>=1.7.0",
        "rich>=13.0.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.2.0",
            "pytest-cov>=2.12.0",
            "black>=21.6b0",
            "flake8>=3.9.0",
            "mypy>=0.910",
        ],
        "docs": [
            "sphinx>=4.0.0",
            "sphinx-rtd-theme>=0.5.2",
        ],
        "web": [
            "flask>=2.0.0",
        ],
        "benchmark": [
            "mappy>=2.24",  # minimap2 bindings for assembly benchmarking
        ],
    },
    entry_points={
        "console_scripts": [
            "viroforge=viroforge.cli:main",
        ],
    },
)
