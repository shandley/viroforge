"""
ViroForge: Synthetic virome data generator
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="viroforge",
    version="0.1.0",
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
    },
    entry_points={
        "console_scripts": [
            "viroforge=viroforge.cli:main",
        ],
    },
)
