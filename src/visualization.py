"""
Visualization module for single-cell analysis.

This module provides functions for creating various plots and visualizations.
"""

import logging
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import anndata as ad
import scanpy as sc
from typing import Optional, List, Tuple

logger = logging.getLogger(__name__)


class Visualizer:
    """Visualization utilities for single-cell analysis."""

    def __init__(self, adata: ad.AnnData, figsize: Tuple[int, int] = (10, 8)):
        """
        Initialize Visualizer object.

        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        figsize : tuple
            Figure size (width, height)
        """
        self.adata = adata
        self.figsize = figsize
        sns.set_style("whitegrid")

    def plot_umap(
        self,
        color: Optional[str] = None,
        save_path: Optional[str] = None,
        **kwargs,
    ) -> None:
        """
        Plot UMAP visualization.

        Parameters
        ----------
        color : str, optional
            Column name for coloring
        save_path : str, optional
            Path to save the figure
        **kwargs
            Additional arguments for sc.pl.umap
        """
        logger.info(f"Plotting UMAP (color={color})")

        plt.figure(figsize=self.figsize)
        sc.pl.umap(self.adata, color=color, save=save_path, **kwargs)
        plt.tight_layout()

        if save_path:
            logger.info(f"Saved figure to {save_path}")

    def plot_violin(
        self,
        keys: List[str],
        groupby: str,
        save_path: Optional[str] = None,
        **kwargs,
    ) -> None:
        """
        Plot violin plots.

        Parameters
        ----------
        keys : list
            Gene names to plot
        groupby : str
            Column name for grouping
        save_path : str, optional
            Path to save the figure
        **kwargs
            Additional arguments
        """
        logger.info(f"Plotting violin plot for {len(keys)} genes")

        plt.figure(figsize=self.figsize)
        sc.pl.violin(self.adata, keys=keys, groupby=groupby, save=save_path, **kwargs)
        plt.tight_layout()

        if save_path:
            logger.info(f"Saved figure to {save_path}")

    def plot_dot(
        self,
        var_names: List[str],
        groupby: str,
        save_path: Optional[str] = None,
        **kwargs,
    ) -> None:
        """
        Plot dotplot.

        Parameters
        ----------
        var_names : list
            Gene names to plot
        groupby : str
            Column name for grouping
        save_path : str, optional
            Path to save the figure
        **kwargs
            Additional arguments
        """
        logger.info(f"Plotting dotplot for {len(var_names)} genes")

        plt.figure(figsize=self.figsize)
        sc.pl.dotplot(
            self.adata, var_names=var_names, groupby=groupby, save=save_path, **kwargs
        )
        plt.tight_layout()

        if save_path:
            logger.info(f"Saved figure to {save_path}")

    def plot_heatmap(
        self,
        var_names: List[str],
        groupby: str,
        save_path: Optional[str] = None,
        **kwargs,
    ) -> None:
        """
        Plot heatmap.

        Parameters
        ----------
        var_names : list
            Gene names to plot
        groupby : str
            Column name for grouping
        save_path : str, optional
            Path to save the figure
        **kwargs
            Additional arguments
        """
        logger.info(f"Plotting heatmap for {len(var_names)} genes")

        plt.figure(figsize=self.figsize)
        sc.pl.heatmap(
            self.adata, var_names=var_names, groupby=groupby, save=save_path, **kwargs
        )
        plt.tight_layout()

        if save_path:
            logger.info(f"Saved figure to {save_path}")

    def plot_qc_metrics(
        self, save_path: Optional[str] = None
    ) -> None:
        """
        Plot QC metrics.

        Parameters
        ----------
        save_path : str, optional
            Path to save the figure
        """
        logger.info("Plotting QC metrics")

        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        # n_genes
        axes[0].hist(self.adata.obs["n_genes"], bins=50)
        axes[0].set_xlabel("Number of genes")
        axes[0].set_ylabel("Frequency")

        # n_counts
        axes[1].hist(np.log10(self.adata.obs["n_counts"] + 1), bins=50)
        axes[1].set_xlabel("Log10(Counts)")
        axes[1].set_ylabel("Frequency")

        # percent_mt
        if "pct_counts_mt" in self.adata.obs:
            axes[2].hist(self.adata.obs["pct_counts_mt"], bins=50)
            axes[2].set_xlabel("Percent MT")
            axes[2].set_ylabel("Frequency")

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            logger.info(f"Saved figure to {save_path}")
