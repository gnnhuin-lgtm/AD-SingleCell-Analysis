"""
Cell type annotation module for single-cell analysis.

This module provides automated cell type annotation and AD-specific gene set analysis.
"""

import logging
from typing import Optional, List, Dict

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

logger = logging.getLogger(__name__)

# Known cell type markers for brain tissues
BRAIN_CELL_MARKERS = {
    "Excitatory_Neuron": [
        "SLC17A7", "SATB2", "CUX2", "RORB", "TBR1", "FEZF2",
        "BCL11B", "CAMK2A", "NRGN", "SYT1",
    ],
    "Inhibitory_Neuron": [
        "GAD1", "GAD2", "SST", "PVALB", "VIP", "CALB1",
        "CALB2", "NPY", "CCK", "NOS1",
    ],
    "Microglia": [
        "TREM2", "TYROBP", "CSF1R", "CX3CR1", "CD33",
        "ITGAM", "FCGR3A", "P2RY12", "AIF1", "C1QA",
    ],
    "Astrocyte": [
        "GFAP", "S100B", "AQP4", "ALDH1L1", "GJA1",
        "SLC1A3", "SLC1A2", "GLUL", "EAAT1", "EAAT2",
    ],
    "Oligodendrocyte": [
        "MBP", "MOG", "PLP1", "MAG", "MOBP",
        "CNP", "OLIG1", "OLIG2", "SOX10", "CLDN11",
    ],
    "OPC": [
        "PDGFRA", "CSPG4", "SOX6", "PCDH15", "VCAN",
        "NG2", "PTPRZ1", "GPR17", "SOX5", "NKX2-2",
    ],
    "Endothelial": [
        "PECAM1", "VWF", "CDH5", "CLDN5", "FLT1",
        "KDR", "TIE1", "ICAM2", "EGFL7", "ESAM",
    ],
    "Pericyte": [
        "PDGFRB", "RGS5", "ANPEP", "ABCC9", "KCNJ8",
        "ACTA2", "CSPG4", "DES", "NOTCH3", "MYH11",
    ],
}


class CellTypeAnnotation:
    """Automated cell type annotation using marker gene expression."""

    def __init__(
        self,
        adata: ad.AnnData,
        marker_dict: Optional[Dict[str, List[str]]] = None,
        min_score: float = 0.1,
    ):
        """
        Initialize CellTypeAnnotation.

        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix (should have log-normalized counts)
        marker_dict : dict, optional
            Dictionary of cell types and their marker genes.
            Defaults to BRAIN_CELL_MARKERS.
        min_score : float
            Minimum score threshold for annotation
        """
        self.adata = adata
        self.marker_dict = marker_dict or BRAIN_CELL_MARKERS
        self.min_score = min_score

    def score_cell_types(self, key_added: str = "cell_type_score") -> pd.DataFrame:
        """
        Score each cell for each cell type using marker gene expression.

        Uses scanpy's sc.tl.score_genes for each cell type.

        Parameters
        ----------
        key_added : str
            Key to store scores in adata.obs

        Returns
        -------
        pd.DataFrame
            Cell type scores (cells x cell types)
        """
        logger.info(
            f"Scoring cell types using {len(self.marker_dict)} reference types"
        )

        score_matrix = pd.DataFrame(index=self.adata.obs_names)

        for cell_type, markers in self.marker_dict.items():
            # Find markers that exist in the data
            valid_markers = [g for g in markers if g in self.adata.var_names]

            if len(valid_markers) < 2:
                logger.warning(
                    f"Fewer than 2 markers found for {cell_type}, skipping"
                )
                continue

            try:
                sc.tl.score_genes(
                    self.adata, gene_list=valid_markers, score_name=f"score_{cell_type}"
                )
                score_matrix[cell_type] = self.adata.obs[f"score_{cell_type}"].values
            except Exception as e:
                logger.warning(f"Failed to score {cell_type}: {e}")
                score_matrix[cell_type] = 0.0

        self.adata.obsm[key_added] = score_matrix.values

        # Assign best cell type
        best_type = score_matrix.idxmax(axis=1)
        best_score = score_matrix.max(axis=1)
        self.adata.obs["cell_type"] = best_type
        self.adata.obs["cell_type_score"] = best_score

        # Mark low-confidence cells as "Unassigned"
        low_conf = best_score < self.min_score
        self.adata.obs.loc[low_conf, "cell_type"] = "Unassigned"
        self.adata.obs.loc[low_conf, "cell_type_score"] = 0.0

        n_assigned = (~low_conf).sum()
        n_total = len(self.adata.obs)
        logger.info(
            f"Assigned cell types: {n_assigned}/{n_total} "
            f"({100 * n_assigned / n_total:.1f}%)"
        )

        return score_matrix

    def plot_cell_type_umap(self, save_path: Optional[str] = None) -> None:
        """Plot UMAP colored by cell type annotation."""
        from src.visualization import Visualizer

        if "X_umap" not in self.adata.obsm:
            raise ValueError("UMAP not computed yet. Run Clustering.umap() first.")

        vis = Visualizer(self.adata)
        vis.plot_umap(color="cell_type", save_path=save_path)


class ADGeneSetAnalyzer:
    """Analysis of AD-specific gene sets across cell types."""

    # Core AD gene sets
    AD_GENE_SETS = {
        "Amyloid_Processing": [
            "APP", "PSEN1", "PSEN2", "BACE1", "BACE2",
            "ADAM10", "APH1A", "APH1B", "NCSTN", "PSENEN",
        ],
        "Tau_Pathology": [
            "MAPT", "GSK3B", "CDK5", "MARK1", "MARK2",
            "DYRK1A", "PP2A", "PIN1", "FKBP5", "Tau",
        ],
        "Neuroinflammation": [
            "IL1B", "IL6", "TNF", "TNFRSF1A", "NFKB1",
            "RELA", "TLR4", "MYD88", "NLRP3", "CCL2",
        ],
        "Oxidative_Stress": [
            "SOD1", "SOD2", "CAT", "GPX1", "GPX4",
            "PRDX1", "PRDX2", "TXNRD1", "GCLC", "GCLM",
        ],
        "Autophagy_Lysosome": [
            "SQSTM1", "LAMP1", "LAMP2", "CTSD", "CTSB",
            "BECN1", "ATG5", "ATG7", "LC3B", "GBA",
        ],
        "Synaptic_Plasticity": [
            "SYN1", "SYN2", "SYT1", "SYP", "DLG4",
            "GRIN1", "GRIN2A", "GRIN2B", "GABRA1", "GABRB2",
        ],
        "APOE_Signaling": [
            "APOE", "CLU", "LRP1", "LRP8", "VLDLR",
            "LDLR", "ABCA1", "ABCG1", "LPL", "TREM2",
        ],
    }

    def __init__(self, adata: ad.AnnData):
        """
        Initialize ADGeneSetAnalyzer.

        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        """
        self.adata = adata
        self.gene_sets = self.AD_GENE_SETS

    def add_custom_gene_sets(self, gene_sets: Dict[str, List[str]]) -> None:
        """Add custom gene sets for analysis."""
        self.gene_sets.update(gene_sets)
        logger.info(f"Added {len(gene_sets)} custom gene sets")

    def score_gene_sets(self, key_added: str = "AD_gene_set_scores") -> pd.DataFrame:
        """
        Score each cell for each AD-relevant gene set.

        Parameters
        ----------
        key_added : str
            Key to store scores

        Returns
        -------
        pd.DataFrame
            Gene set scores per cell
        """
        logger.info(f"Scoring {len(self.gene_sets)} AD gene sets")

        score_matrix = pd.DataFrame(index=self.adata.obs_names)

        for gs_name, genes in self.gene_sets.items():
            valid_genes = [g for g in genes if g in self.adata.var_names]

            if len(valid_genes) < 3:
                logger.debug(f"Fewer than 3 valid genes for {gs_name}, skipping")
                continue

            try:
                sc.tl.score_genes(
                    self.adata,
                    gene_list=valid_genes,
                    score_name=f"AD_{gs_name}",
                )
                score_matrix[gs_name] = self.adata.obs[f"AD_{gs_name}"].values
            except Exception as e:
                logger.warning(f"Failed to score {gs_name}: {e}")
                score_matrix[gs_name] = 0.0

        self.adata.obsm[key_added] = score_matrix.values
        return score_matrix

    def get_gene_set_summary(
        self, groupby: str = "cell_type"
    ) -> pd.DataFrame:
        """
        Get average gene set scores per group (e.g., per cell type).

        Parameters
        ----------
        groupby : str
            Column in adata.obs to group by

        Returns
        -------
        pd.DataFrame
            Mean gene set scores per group
        """
        score_cols = [
            col for col in self.adata.obs.columns if col.startswith("AD_")
        ]

        if not score_cols:
            logger.warning("No AD gene set scores found. Run score_gene_sets() first.")
            return pd.DataFrame()

        summary = self.adata.obs.groupby(groupby)[score_cols].mean()

        # Rename columns back for readability
        summary.columns = [
            col.replace("AD_", "") for col in summary.columns
        ]

        return summary

    def compare_ad_vs_control(
        self,
        ad_group: str = "AD",
        control_group: str = "Control",
        condition_key: str = "condition",
        groupby: str = "cell_type",
    ) -> pd.DataFrame:
        """
        Compare AD gene set scores between AD and control groups per cell type.

        Parameters
        ----------
        ad_group : str
            Label for AD group
        control_group : str
            Label for control group
        condition_key : str
            Column in adata.obs indicating condition
        groupby : str
            Column in adata.obs for cell type grouping

        Returns
        -------
        pd.DataFrame
            Differential gene set activity
        """
        logger.info(f"Comparing AD vs Control gene set activity")

        score_cols = [
            col for col in self.adata.obs.columns if col.startswith("AD_")
        ]

        results = []

        for cell_type in self.adata.obs[groupby].unique():
            mask_type = self.adata.obs[groupby] == cell_type
            mask_ad = self.adata.obs[condition_key] == ad_group
            mask_ctrl = self.adata.obs[condition_key] == control_group

            for score_col in score_cols:
                ad_scores = self.adata.obs.loc[mask_type & mask_ad, score_col]
                ctrl_scores = self.adata.obs.loc[mask_type & mask_ctrl, score_col]

                if len(ad_scores) < 5 or len(ctrl_scores) < 5:
                    continue

                from scipy.stats import mannwhitneyu

                stat, pval = mannwhitneyu(ad_scores, ctrl_scores, alternative="two-sided")

                results.append({
                    "cell_type": cell_type,
                    "gene_set": score_col.replace("AD_", ""),
                    "ad_mean": ad_scores.mean(),
                    "control_mean": ctrl_scores.mean(),
                    "log2_fc": np.log2(ad_scores.mean() / max(ctrl_scores.mean(), 1e-10)),
                    "p_value": pval,
                    "significant": pval < 0.05,
                })

        return pd.DataFrame(results).sort_values("p_value")
