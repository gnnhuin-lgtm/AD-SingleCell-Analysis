"""
AD-SingleCell-Analysis: Single-cell RNA-seq analysis for Alzheimer's Disease research

This package provides tools and workflows for preprocessing, analyzing, and visualizing
single-cell transcriptomic data in the context of Alzheimer's Disease research.
"""

__version__ = "0.1.0"
__author__ = "gnnhuin-lgtm"
__email__ = "gnnhuin-lgtm@github.com"

from . import preprocessing
from . import analysis
from . import visualization
from . import utils

__all__ = [
    "preprocessing",
    "analysis",
    "visualization",
    "utils",
]
