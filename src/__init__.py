"""
AD-SingleCell-Analysis: Single-cell transcriptomics analysis for Alzheimer's Disease.

This package provides a complete workflow for analyzing single-cell RNA-seq data
from Alzheimer's Disease (AD) studies, including quality control, preprocessing,
clustering, cell type annotation, differential expression, and visualization.
"""

__version__ = "0.1.0"
__author__ = "gnnhuin-lgtm"

from src.preprocessing import QualityControl, Normalize, VariableGeneSelection
from src.analysis import Clustering, DifferentialExpression
from src.annotation import CellTypeAnnotation, ADGeneSetAnalyzer
from src.visualization import Visualizer
from src.utils import setup_logging, load_config, save_config, create_directories
from src.pipeline import run_pipeline, ADSingleCellPipeline

__all__ = [
    # Version
    "__version__",
    "__author__",
    # Classes
    "QualityControl",
    "Normalize",
    "VariableGeneSelection",
    "Clustering",
    "DifferentialExpression",
    "CellTypeAnnotation",
    "ADGeneSetAnalyzer",
    "Visualizer",
    "ADSingleCellPipeline",
    # Functions
    "setup_logging",
    "load_config",
    "save_config",
    "create_directories",
    "run_pipeline",
]
