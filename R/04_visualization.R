#' 数据可视化模块
#'
#' 提供高质量的生信分析图表绘制函数
#'
#' @keywords internal

#' 绘制火山图
#'
#' @param de_result 差异表达结果数据框
#' @param fc_threshold log2FC阈值
#' @param pval_threshold p值阈值
#' @param title 图表标题
#'
#' @return ggplot2对象
#'
#' @examples
#' \dontrun{
#' p <- plot_volcano(de_result, fc_threshold = 1, pval_threshold = 0.05)
#' print(p)
#' }
#'
#' @export
plot_volcano <- function(de_result, fc_threshold = 1, pval_threshold = 0.05, 
                         title = "Volcano Plot") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先安装 ggplot2 包")
  }
  
  de_result$Color <- ifelse(de_result$padj < pval_threshold & 
                             abs(de_result$log2FC) > fc_threshold,
                            ifelse(de_result$log2FC > 0, "Up-regulated", "Down-regulated"),
                            "Not significant")
  
  p <- ggplot2::ggplot(de_result, ggplot2::aes(x = log2FC, y = -log10(padj), color = Color)) +
    ggplot2::geom_point(alpha = 0.6, size = 2) +
    ggplot2::scale_color_manual(values = c("Up-regulated" = "red", 
                                           "Down-regulated" = "blue",
                                           "Not significant" = "gray")) +
    ggplot2::geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
                        linetype = "dashed", color = "black", alpha = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(pval_threshold), 
                        linetype = "dashed", color = "black", alpha = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "log2(Fold Change)", y = "-log10(adjusted p-value)",
                  color = "Regulation") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

#' 绘制MA图
#'
#' @param de_result 差异表达结果数据框
#' @param fc_threshold log2FC阈值
#' @param title 图表标题
#'
#' @return ggplot2对象
#'
#' @examples
#' \dontrun{
#' p <- plot_ma(de_result)
#' print(p)
#' }
#'
#' @export
plot_ma <- function(de_result, fc_threshold = 1, title = "MA Plot") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先安装 ggplot2 包")
  }
  
  de_result$A <- rowMeans(log2(de_result[, grep("count", names(de_result), ignore.case = TRUE)] + 1))
  
  de_result$Color <- ifelse(abs(de_result$log2FC) > fc_threshold, "Significant", "Not significant")
  
  p <- ggplot2::ggplot(de_result, ggplot2::aes(x = A, y = log2FC, color = Color)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "black", alpha = 0.3) +
    ggplot2::geom_hline(yintercept = c(-fc_threshold, fc_threshold), 
                        linetype = "dashed", color = "red", alpha = 0.5) +
    ggplot2::scale_color_manual(values = c("Significant" = "red", 
                                           "Not significant" = "gray")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "A (Average log2 intensity)", 
                  y = "M (log2 Fold Change)", color = "") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

#' 绘制热图
#'
#' @param data 表达矩阵
#' @param annotation 注释信息（可选）
#' @param title 图表标题
#'
#' @return ComplexHeatmap对象
#'
#' @examples
#' \dontrun{
#' p <- plot_heatmap(expr_matrix, title = "Gene Expression Heatmap")
#' print(p)
#' }
#'
#' @export
plot_heatmap <- function(data, annotation = NULL, title = "Heatmap") {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("请先安装 ComplexHeatmap 包: BiocManager::install('ComplexHeatmap')")
  }
  
  # 标准化数据
  data_scaled <- t(scale(t(data)))
  
  # 创建热图
  heatmap <- ComplexHeatmap::Heatmap(
    data_scaled,
    name = "Z-score",
    row_title = "Genes",
    column_title = "Samples",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = TRUE
  )
  
  return(heatmap)
}

#' 绘制PCA图
#'
#' @param data 表达矩阵
#' @param metadata 样本元数据
#' @param color_by 按哪个列进行着色
#'
#' @return ggplot2对象
#'
#' @examples
#' \dontrun{
#' p <- plot_pca(expr_matrix, metadata, color_by = "group")
#' print(p)
#' }
#'
#' @export
plot_pca <- function(data, metadata = NULL, color_by = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先安装 ggplot2 包")
  }
  
  # PCA分析
  pca_result <- stats::prcomp(t(data), scale. = TRUE)
  
  # 创建PCA结果数据框
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    sample = rownames(pca_result$x)
  )
  
  if (!is.null(metadata) && !is.null(color_by)) {
    pca_df <- merge(pca_df, metadata, by.x = "sample", by.y = "row.names", all.x = TRUE)
  }
  
  # 计算方差解释率
  var_exp <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
  
  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(size = 3)
  
  if (!is.null(color_by) && color_by %in% colnames(pca_df)) {
    p <- p + ggplot2::aes_string(color = color_by)
  }
  
  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "PCA Plot",
      x = paste0("PC1 (", round(var_exp[1], 2), "%)"),
      y = paste0("PC2 (", round(var_exp[2], 2), "%)"),
      color = color_by
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

#' 绘制箱线图
#'
#' @param data 表达矩阵或数据框
#' @param group_info 分组信息向量
#' @param title 图表标题
#'
#' @return ggplot2对象
#'
#' @examples
#' \dontrun{
#' p <- plot_boxplot(expr_matrix, group_info = group_vector)
#' print(p)
#' }
#'
#' @export
plot_boxplot <- function(data, group_info, title = "Boxplot") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先安装 ggplot2 包")
  }
  
  # 准备数据
  plot_data <- data.frame(
    value = as.numeric(data),
    group = rep(group_info, each = nrow(data)),
    gene = rep(rownames(data), ncol(data))
  )
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = value, fill = group)) +
    ggplot2::geom_boxplot(alpha = 0.6) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Group", y = "Expression Level", fill = "Group") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}