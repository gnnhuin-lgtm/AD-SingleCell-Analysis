"""
Analysis module for single-cell data.

This module provides functions for clustering, annotation, and differential expression analysis.
"""

import logging
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from typing import Optional, List

logger = logging.getLogger(__name__)


class Clustering:
    """Clustering analysis for single-cell data."""

    def __init__(self, adata: ad.AnnData, config: Optional[dict] = None):
        """
        Initialize Clustering object.

        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        config : dict, optional
            Configuration dictionary
        """
        self.adata = adata
        self.config = config or {}

    def pca(self, n_pcs: int = 50, random_state: int = 42) -> ad.AnnData:
        """
        Perform PCA.

        Parameters
        ----------
        n_pcs : int
            Number of principal components
        random_state : int
            Random seed

        Returns
        -------
        ad.AnnData
            Data with PCA results
        """
        logger.info(f"Performing PCA (n_pcs={n_pcs})")
        sc.pp.pca(self.adata, n_comps=n_pcs, random_state=random_state)
        return self.adata

    def umap(
        self, n_neighbors: int = 15, min_dist: float = 0.1, metric: str = "euclidean"
    ) -> ad.AnnData:
        """
        Perform UMAP dimensionality reduction.

        Parameters
        ----------
        n_neighbors : int
            Number of neighbors for UMAP
        min_dist : float
            Minimum distance parameter
        metric : str
            Distance metric

        Returns
        -------
        ad.AnnData
            Data with UMAP coordinates
        """
        logger.info(f"Computing UMAP (n_neighbors={n_neighbors}, min_dist={min_dist})")

        if "X_pca" not in self.adata.obsm:
            logger.warning("PCA not found, computing PCA first")
            self.pca()

        sc.pp.neighbors(
            self.adata, n_neighbors=n_neighbors, use_rep="X_pca", metric=metric
        )
        sc.tl.umap(self.adata, min_dist=min_dist)

        return self.adata

    def leiden(
        self, resolution: float = 0.5, random_state: int = 42
    ) -> ad.AnnData:
        """
        Perform Leiden clustering.

        Parameters
        ----------
        resolution : float
            Resolution parameter for clustering
        random_state : int
            Random seed

        Returns
        -------
        ad.AnnData
            Data with clustering results
        """
        logger.info(f"Performing Leiden clustering (resolution={resolution})")

        if "neighbors" not in self.adata.obsp:
            logger.warning("Neighbors not computed, computing neighbors first")
            sc.pp.neighbors(self.adata)

        sc.tl.leiden(
            self.adata, resolution=resolution, random_state=random_state, key_added="leiden"
        )

        logger.info(f"Found {self.adata.obs['leiden'].nunique()} clusters")
        return self.adata

    def louvain(
        self, resolution: float = 1.0, random_state: int = 42
    ) -> ad.AnnData:
        """
        Perform Louvain clustering.

        Parameters
        ----------
        resolution : float
            Resolution parameter
        random_state : int
            Random seed

        Returns
        -------
        ad.AnnData
            Data with clustering results
        """
        logger.info(f"Performing Louvain clustering (resolution={resolution})")

        if "neighbors" not in self.adata.obsp:
            logger.warning("Neighbors not computed, computing neighbors first")
            sc.pp.neighbors(self.adata)

        sc.tl.louvain(
            self.adata,
            resolution=resolution,
            random_state=random_state,
            key_added="louvain",
        )

        logger.info(f"Found {self.adata.obs['louvain'].nunique()} clusters")
        return self.adata


class DifferentialExpression:
    """Differential expression analysis."""

    def __init__(self, adata: ad.AnnData, config: Optional[dict] = None):
        """
        Initialize DifferentialExpression object.

        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        config : dict, optional
            Configuration dictionary
        """
        self.adata = adata
        self.config = config or {}

    def run_de(
        self,
        groupby: str,
        group1: str,
        group2: str,
        method: str = "wilcoxon",
    ) -> pd.DataFrame:
        """
        Run differential expression analysis.

        Parameters
        ----------
        groupby : str
            Column name for grouping
        group1 : str
            First group
        group2 : str
            Second group
        method : str
            Method for DE analysis

        Returns
        -------
        pd.DataFrame
            Differential expression results
        """
        logger.info(f"Running DE analysis: {group1} vs {group2} (method={method})")

        # Create comparison groups
        comparison_groups = [group1, group2]

        # Subset data
        adata_sub = self.adata[self.adata.obs[groupby].isin(comparison_groups)]

        # Run DE analysis
        sc.tl.rank_genes_groups(
            adata_sub,
            groupby=groupby,
            groups=[group2],
            reference=group1,
            method=method,
        )

        # Extract results
        result = sc.get.rank_genes_groups_df(adata_sub, group=group2)
        logger.info(f"Found {len(result)} differentially expressed genes")

        return result
