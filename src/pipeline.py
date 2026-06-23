"""
End-to-end pipeline for AD single-cell analysis.

This module provides a complete analysis pipeline that orchestrates all steps
from data loading through preprocessing, clustering, annotation, and reporting.
"""

import logging
import os
from pathlib import Path
from typing import Optional, Dict, Any

import anndata as ad
import scanpy as sc

from src.preprocessing import QualityControl, Normalize, VariableGeneSelection
from src.analysis import Clustering, DifferentialExpression
from src.annotation import CellTypeAnnotation, ADGeneSetAnalyzer
from src.visualization import Visualizer
from src.utils import setup_logging, load_config, create_directories

logger = logging.getLogger(__name__)

# Suppress scanpy figure saving (we handle it ourselves)
sc.settings.autosave = False
sc.settings.autoshow = False


def run_pipeline(
    input_file: str,
    config_path: str = "config/config.yaml",
    output_dir: str = "results",
) -> ad.AnnData:
    """
    Run the complete AD single-cell analysis pipeline.

    Parameters
    ----------
    input_file : str
        Path to input h5ad file
    config_path : str
        Path to YAML configuration file
    output_dir : str
        Base output directory

    Returns
    -------
    ad.AnnData
        Fully processed AnnData object
    """
    pipeline = ADSingleCellPipeline(config_path, output_dir)
    result = pipeline.load_data(input_file).run()
    return result


class ADSingleCellPipeline:
    """Orchestrated pipeline for AD single-cell RNA-seq analysis."""

    def __init__(
        self,
        config_path: str = "config/config.yaml",
        output_dir: str = "results",
    ):
        """
        Initialize the pipeline.

        Parameters
        ----------
        config_path : str
            Path to YAML configuration
        output_dir : str
            Base output directory
        """
        self.config = load_config(config_path) if os.path.exists(config_path) else {}
        self.output_dir = Path(output_dir)
        self.adata: Optional[ad.AnnData] = None

        # Setup
        log_cfg = self.config.get("logging", {})
        setup_logging(
            log_file=os.path.join(output_dir, log_cfg.get("file", "logs/analysis.log")),
            level=log_cfg.get("level", "INFO"),
        )
        create_directories(output_dir)

    def load_data(self, file_path: str) -> "ADSingleCellPipeline":
        """
        Load data from a file.

        Parameters
        ----------
        file_path : str
            Path to input file (.h5ad, .h5, .csv, .loom)

        Returns
        -------
        ADSingleCellPipeline
            Self for chaining
        """
        logger.info(f"Loading data from {file_path}")
        ext = Path(file_path).suffix.lower()

        if ext == ".h5ad":
            self.adata = ad.read_h5ad(file_path)
        elif ext == ".h5":
            self.adata = ad.read_h5(file_path)
        elif ext in (".csv", ".tsv"):
            import pandas as pd
            df = pd.read_csv(file_path, index_col=0)
            self.adata = ad.AnnData(X=df.values, obs=pd.DataFrame(index=df.index))
            self.adata.var_names = df.columns
        elif ext == ".loom":
            self.adata = sc.read_loom(file_path)
        else:
            raise ValueError(f"Unsupported file format: {ext}")

        logger.info(f"Loaded data: {self.adata.shape[0]} cells x {self.adata.shape[1]} genes")
        return self

    def run_quality_control(self) -> "ADSingleCellPipeline":
        """Run quality control steps."""
        logger.info("=== Step 1: Quality Control ===")
        qc_config = self.config.get("quality_control", {})

        qc = QualityControl(self.adata)
        qc.calculate_qc_metrics()

        # Plot QC metrics
        vis = Visualizer(self.adata)
        qc_plot_path = str(self.output_dir / "figures" / "qc_metrics.png")
        vis.plot_qc_metrics(save_path=qc_plot_path)

        # Filter cells
        self.adata = qc.filter_cells(
            min_genes=qc_config.get("min_genes", 200),
            max_genes=qc_config.get("max_genes", 10000),
            percent_mt=qc_config.get("percent_mt", 20),
        )

        # Filter genes
        self.adata = qc.filter_genes(
            min_cells=qc_config.get("min_cells", 3),
        )

        logger.info(f"After QC: {self.adata.shape[0]} cells x {self.adata.shape[1]} genes")
        return self

    def run_normalization(self) -> "ADSingleCellPipeline":
        """Run normalization steps."""
        logger.info("=== Step 2: Normalization ===")
        norm_config = self.config.get("normalization", {})

        normalizer = Normalize(self.adata)
        self.adata = normalizer.normalize_and_log(
            target_sum=norm_config.get("target_sum", 1e4),
        )

        # HVG selection
        hvg_config = self.config.get("hvg", {})
        hvg = VariableGeneSelection(self.adata)
        self.adata = hvg.select_hvg(
            n_top_genes=hvg_config.get("n_top_genes", 2000),
            flavor=hvg_config.get("flavor", "seurat_v3"),
        )

        # Subset to HVGs for downstream analysis
        self.adata = self.adata[:, self.adata.var["highly_variable"]].copy()

        logger.info(f"After HVG selection: {self.adata.shape[0]} cells x {self.adata.shape[1]} genes")
        return self

    def run_dimensionality_reduction(self) -> "ADSingleCellPipeline":
        """Run dimensionality reduction and clustering."""
        logger.info("=== Step 3: Dimensionality Reduction & Clustering ===")
        dr_config = self.config.get("dim_reduction", {})
        cluster_config = self.config.get("clustering", {})
        umap_config = self.config.get("umap", {})

        clustering = Clustering(self.adata)

        # PCA
        self.adata = clustering.pca(
            n_pcs=dr_config.get("n_pcs", 50),
            random_state=dr_config.get("random_state", 42),
        )

        # UMAP
        self.adata = clustering.umap(
            n_neighbors=umap_config.get("n_neighbors", 15),
            min_dist=umap_config.get("min_dist", 0.1),
            metric=umap_config.get("metric", "euclidean"),
        )

        # Clustering
        method = cluster_config.get("method", "leiden")
        if method == "leiden":
            self.adata = clustering.leiden(
                resolution=cluster_config.get("resolution", 0.5),
                random_state=cluster_config.get("random_state", 42),
            )
        else:
            self.adata = clustering.louvain(
                resolution=cluster_config.get("resolution", 1.0),
                random_state=cluster_config.get("random_state", 42),
            )

        # Plot UMAP
        vis = Visualizer(self.adata)
        cluster_key = method  # "leiden" or "louvain"
        vis.plot_umap(
            color=cluster_key,
            save_path=str(self.output_dir / "figures" / f"umap_{method}.png"),
        )

        return self

    def run_annotation(self) -> "ADSingleCellPipeline":
        """Run cell type annotation."""
        logger.info("=== Step 4: Cell Type Annotation ===")

        annotator = CellTypeAnnotation(self.adata)
        _ = annotator.score_cell_types()

        # Plot annotated UMAP
        vis = Visualizer(self.adata)
        vis.plot_umap(
            color="cell_type",
            save_path=str(self.output_dir / "figures" / "umap_annotated.png"),
        )

        return self

    def run_ad_analysis(self) -> "ADSingleCellPipeline":
        """Run AD-specific gene set analysis."""
        logger.info("=== Step 5: AD-Specific Analysis ===")

        ad_analyzer = ADGeneSetAnalyzer(self.adata)
        _ = ad_analyzer.score_gene_sets()

        # Summary per cell type
        summary = ad_analyzer.get_gene_set_summary(groupby="cell_type")
        summary_path = str(self.output_dir / "tables" / "ad_gene_set_summary.csv")
        summary.to_csv(summary_path)
        logger.info(f"AD gene set summary saved to {summary_path}")

        return self

    def run_differential_expression(self) -> "ADSingleCellPipeline":
        """Run differential expression analysis if condition metadata exists."""
        logger.info("=== Step 6: Differential Expression ===")

        if "condition" not in self.adata.obs:
            logger.warning(
                "No 'condition' column found in adata.obs. "
                "Skipping differential expression. "
                "Add 'condition' with values like 'AD'/'Control' to enable this step."
            )
            return self

        de_cfg = self.config.get("differential_expression", {})
        de = DifferentialExpression(self.adata)

        groups = self.adata.obs["condition"].unique()
        if len(groups) < 2:
            logger.warning("Need at least 2 conditions for DE analysis")
            return self

        result = de.run_de(
            groupby="condition",
            group1=groups[0],
            group2=groups[1],
            method=de_cfg.get("method", "wilcoxon"),
        )

        # Filter significant genes
        significant = result[
            (result["pvals_adj"] < de_cfg.get("max_pval", 0.05))
            & (abs(result["logfoldchanges"]) > de_cfg.get("min_log2fc", 0.25))
        ]

        de_path = str(self.output_dir / "tables" / "differential_expression.csv")
        significant.to_csv(de_path, index=False)
        logger.info(f"DE results saved to {de_path}")

        return self

    def save_results(self, prefix: str = "final") -> "ADSingleCellPipeline":
        """Save the processed AnnData object."""
        output_path = str(self.output_dir / f"{prefix}_adata.h5ad")
        self.adata.write(output_path)
        logger.info(f"Final AnnData saved to {output_path}")
        return self

    def run(self) -> ad.AnnData:
        """
        Execute the complete pipeline.

        Returns
        -------
        ad.AnnData
            Fully processed AnnData object
        """
        if self.adata is None:
            raise ValueError("No data loaded. Call load_data() first.")

        self.run_quality_control()
        self.run_normalization()
        self.run_dimensionality_reduction()
        self.run_annotation()
        self.run_ad_analysis()
        self.run_differential_expression()
        self.save_results()

        logger.info("=" * 50)
        logger.info("Pipeline completed successfully!")
        logger.info("=" * 50)

        return self.adata
