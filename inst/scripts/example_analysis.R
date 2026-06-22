#!/usr/bin/env Rscript
#
# RNA-seq完整分析流程示例
# 使用生信技能包进行端到端的分析
#

library(bioinformaticsSkills)
library(tidyverse)

# ============================================================================
# 第1部分：数据准备
# ============================================================================

cat("\n========== 数据读取和预处理 ==========\n")

# 创建输出目录
dir.create("results", showWarnings = FALSE)
dir.create("results/qc", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)

# 这里应该替换为实际数据文件
# counts <- read_count_matrix("data/expression_matrix.csv")
# metadata <- read.csv("data/sample_metadata.csv", row.names = 1)

# 为演示使用模拟数据
set.seed(123)
n_genes <- 1000
n_samples <- 20

counts <- matrix(
  rpois(n_genes * n_samples, lambda = 10),
  nrow = n_genes,
  ncol = n_samples
)
rownames(counts) <- paste0("ENSG", sprintf("%010d", 1:n_genes))
colnames(counts) <- paste0("Sample_", 1:n_samples)

metadata <- data.frame(
  sample = colnames(counts),
  group = rep(c("control", "treatment"), each = 10),
  batch = rep(c("batch1", "batch2"), 10)
)
rownames(metadata) <- metadata$sample

cat("✓ 成功读取:", nrow(counts), "个基因,", ncol(counts), "个样本\n")

# ============================================================================
# 第2部分：质量控制
# ============================================================================

cat("\n========== 质量控制 ==========\n")

# 计算QC指标
qc_metrics <- calculate_qc_metrics(counts, metadata)
head(qc_metrics)

# 绘制QC报告
plot_qc_report(qc_metrics, output_dir = "results/qc/")

# 过滤低质量样本
filtered_counts <- filter_low_quality_samples(
  counts, qc_metrics,
  min_features = 500,
  min_counts = 5000
)

cat("✓ 过滤后保留:", ncol(filtered_counts), "个样本\n")

# ============================================================================
# 第3部分：数据标准化
# ============================================================================

cat("\n========== 数据标准化 ==========\n")

# CPM标准化
cpm_normalized <- normalize_data(filtered_counts, method = "cpm")

# 对数变换
log_counts <- log_transform(cpm_normalized, pseudocount = 1)

cat("✓ 数据标准化完成\n")

# ============================================================================
# 第4部分：批次效应校正（可选）
# ============================================================================

cat("\n========== 批次效应校正 ==========\n")

# 使用ComBat校正批次效应
batch_corrected <- batch_correction(
  log_counts,
  batch = metadata[colnames(log_counts), "batch"],
  method = "ComBat"
)

cat("✓ 批次效应校正完成\n")

# ============================================================================
# 第5部分：差异表达分析
# ============================================================================

cat("\n========== 差异表达分析 ==========\n")

# 提取对应样本的元数据
meta_filtered <- metadata[colnames(filtered_counts), ]

# 进行DE分析
de_results <- differential_expression(
  filtered_counts,
  metadata = meta_filtered,
  group_col = "group",
  contrast = c("treatment", "control"),
  method = "DESeq2"
)

# 统计显著基因
sig_genes <- de_results[de_results$padj < 0.05, ]
up_genes <- sum(sig_genes$log2FC > 1)
down_genes <- sum(sig_genes$log2FC < -1)

cat("✓ 差异表达分析完成\n")
cat("  显著基因总数:", nrow(sig_genes), "\n")
cat("  上调基因:", up_genes, "\n")
cat("  下调基因:", down_genes, "\n")

# ============================================================================
# 第6部分：可视化
# ============================================================================

cat("\n========== 结果可视化 ==========\n")

# 火山图
volcano_plot <- plot_volcano(
  de_results,
  fc_threshold = 1,
  pval_threshold = 0.05,
  title = "差异表达分析 - 火山图"
)
ggplot2::ggsave("results/plots/volcano_plot.pdf", volcano_plot, width = 8, height = 6)
cat("✓ 火山图已保存\n")

# MA图
ma_plot <- plot_ma(de_results, fc_threshold = 1, title = "MA图")
ggplot2::ggsave("results/plots/ma_plot.pdf", ma_plot, width = 8, height = 6)
cat("✓ MA图已保存\n")

# PCA图
pca_plot <- plot_pca(log_counts, meta_filtered, color_by = "group")
ggplot2::ggsave("results/plots/pca_plot.pdf", pca_plot, width = 8, height = 6)
cat("✓ PCA图已保存\n")

# 热图（选择top基因）
top_genes <- de_results[order(de_results$padj), ][1:50, "gene_id"]
heatmap_data <- log_counts[top_genes, ]
heatmap <- plot_heatmap(heatmap_data, title = "Top 50 差异表达基因热图")
grDevices::pdf("results/plots/heatmap.pdf", width = 10, height = 8)
print(heatmap)
grDevices::dev.off()
cat("✓ 热图已保存\n")

# ============================================================================
# 第7部分：富集分析
# ============================================================================

cat("\n========== 富集分析 ==========\n")

# 提取显著基因
sig_gene_ids <- de_results[de_results$padj < 0.05 & abs(de_results$log2FC) > 1, "gene_id"]

# GO富集分析
go_enrichment <- enrichment_analysis(
  sig_gene_ids,
  analysis_type = "GO",
  ont = "BP"
)

head(go_enrichment)
cat("✓ GO富集分析完成\n")

# ============================================================================
# 第8部分：结果导出
# ============================================================================

cat("\n========== 结果导出 ==========\n")

# 保存DE结果
save_result(de_results, "results/differential_expression_results.csv", format = "csv")
cat("✓ DE结果已保存\n")

# 保存富集分析结果
save_result(go_enrichment, "results/go_enrichment_results.csv", format = "csv")
cat("✓ GO富集结果已保存\n")

# 保存QC指标
save_result(qc_metrics, "results/qc_metrics.csv", format = "csv")
cat("✓ QC指标已保存\n")

# ============================================================================
# 完成
# ============================================================================

cat("\n========== 分析完成 ==========\n")
cat("所有结果已保存到 results/ 目录\n")
