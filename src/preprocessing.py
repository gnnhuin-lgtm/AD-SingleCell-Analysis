"""
Data preprocessing module for single-cell analysis.

This module provides functions for quality control, normalization, and preprocessing
of single-cell RNA-seq data.
"""

import logging
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from typing import Optional, Tuple

logger = logging.getLogger(__name__)


class QualityControl:
    """Quality control for single-cell data."""

    def __init__(self, adata: ad.AnnData, config: Optional[dict] = None):
        """
        Initialize QualityControl object.

        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        config : dict, optional
            Configuration dictionary with QC parameters
        """
        self.adata = adata
        self.config = config or {}

    def calculate_qc_metrics(self) -> ad.AnnData:
        """Calculate QC metrics for cells and genes."""
        logger.info("Calculating QC metrics...")

        # Add mitochondrial gene percentage
        self.adata.var["mt"] = self.adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            self.adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False
        )

        logger.info(f"QC metrics calculated. Data shape: {self.adata.shape}")
        return self.adata

    def filter_cells(
        self,
        min_genes: int = 200,
        max_genes: int = 10000,
        percent_mt: float = 20,
    ) -> ad.AnnData:
        """
        Filter cells based on QC metrics.

        Parameters
        ----------
        min_genes : int
            Minimum number of genes per cell
        max_genes : int
            Maximum number of genes per cell
        percent_mt : float
            Maximum percentage of mitochondrial genes

        Returns
        -------
        ad.AnnData
            Filtered data
        """
        logger.info(
            f"Filtering cells (min_genes={min_genes}, max_genes={max_genes}, percent_mt={percent_mt})"
        )

        initial_count = self.adata.n_obs

        sc.pp.filter_cells(self.adata, min_genes=min_genes)
        sc.pp.filter_cells(self.adata, max_genes=max_genes)

        self.adata = self.adata[self.adata.obs["pct_counts_mt"] < percent_mt]

        filtered_count = self.adata.n_obs
        logger.info(
            f"Cells filtered: {initial_count} -> {filtered_count} ({100*filtered_count/initial_count:.1f}%)"
        )

        return self.adata

    def filter_genes(self, min_cells: int = 3) -> ad.AnnData:
        """
        Filter genes based on expression.

        Parameters
        ----------
        min_cells : int
            Minimum number of cells expressing the gene

        Returns
        -------
        ad.AnnData
            Filtered data
        """
        logger.info(f"Filtering genes (min_cells={min_cells})")

        initial_count = self.adata.n_vars
        sc.pp.filter_genes(self.adata, min_cells=min_cells)
        filtered_count = self.adata.n_vars

        logger.info(
            f"Genes filtered: {initial_count} -> {filtered_count} ({100*filtered_count/initial_count:.1f}%)"
        )

        return self.adata


class Normalize:
    """Normalization for single-cell data."""

    def __init__(self, adata: ad.AnnData, config: Optional[dict] = None):
        """
        Initialize Normalize object.

        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        config : dict, optional
            Configuration dictionary with normalization parameters
        """
        self.adata = adata
        self.config = config or {}

    def normalize(self, target_sum: float = 1e4) -> ad.AnnData:
        """
        Normalize library sizes.

        Parameters
        ----------
        target_sum : float
            Target sum for normalization

        Returns
        -------
        ad.AnnData
            Normalized data
        """
        logger.info(f"Normalizing data (target_sum={target_sum})")
        sc.pp.normalize_total(self.adata, target_sum=target_sum)
        return self.adata

    def log_transform(self, base: int = 2) -> ad.AnnData:
        """
        Apply log transformation.

        Parameters
        ----------
        base : int
            Logarithm base

        Returns
        -------
        ad.AnnData
            Log-transformed data
        """
        logger.info(f"Applying log transformation (base={base})")
        sc.pp.log1p(self.adata)
        return self.adata

    def normalize_and_log(self, target_sum: float = 1e4) -> ad.AnnData:
        """Apply both normalization and log transformation."""
        self.normalize(target_sum)
        self.log_transform()
        return self.adata


class VariableGeneSelection:
    """Highly variable gene selection."""

    def __init__(self, adata: ad.AnnData, config: Optional[dict] = None):
        """
        Initialize VariableGeneSelection object.

        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        config : dict, optional
            Configuration dictionary
        """
        self.adata = adata
        self.config = config or {}

    def select_hvg(
        self, n_top_genes: int = 2000, flavor: str = "seurat_v3"
    ) -> ad.AnnData:
        """
        Select highly variable genes.

        Parameters
        ----------
        n_top_genes : int
            Number of highly variable genes to select
        flavor : str
            Method for HVG selection

        Returns
        -------
        ad.AnnData
            Data with HVG selection
        """
        logger.info(f"Selecting HVGs (n_top_genes={n_top_genes}, flavor={flavor})")
        sc.pp.highly_variable_genes(
            self.adata, n_top_genes=n_top_genes, subset=False, flavor=flavor
        )
        return self.adata
