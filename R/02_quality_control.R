#' 质量控制模块
#'
#' 提供QC指标计算和质量评估函数
#'
#' @keywords internal

#' 计算QC指标
#'
#' @param counts 计数矩阵
#' @param metadata 样本元数据
#'
#' @return 包含QC指标的数据框
#'
#' @details
#' 计算的QC指标包括：
#' - nFeatures: 检测到的特征数
#' - nCounts: 总计数
#' - percentZero: 零值百分比
#' - librarySize: 文库大小
#'
#' @examples
#' \dontrun{
#' qc_metrics <- calculate_qc_metrics(counts, metadata)
#' }
#'
#' @export
calculate_qc_metrics <- function(counts, metadata = NULL) {
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop("counts 必须是矩阵或数据框")
  }
  
  # 转换为数值矩阵
  counts <- as.matrix(counts)
  
  qc_df <- data.frame(
    sample = colnames(counts),
    nFeatures = colSums(counts > 0),
    nCounts = colSums(counts),
    nZeros = colSums(counts == 0),
    percentZero = (colSums(counts == 0) / nrow(counts)) * 100,
    medianCount = apply(counts, 2, median),
    meanCount = colMeans(counts),
    row.names = NULL
  )
  
  if (!is.null(metadata)) {
    qc_df <- merge(qc_df, metadata, by = "row.names")
    rownames(qc_df) <- qc_df$row.names
    qc_df$row.names <- NULL
  }
  
  message("QC指标计算完成")
  return(qc_df)
}

#' 绘制QC报告
#'
#' @param qc_metrics QC指标数据框
#' @param output_dir 输出目录
#'
#' @return 隐形返回生成的图表文件列表
#'
#' @examples
#' \dontrun{
#' plot_qc_report(qc_metrics, output_dir = "./results/qc/")
#' }
#'
#' @export
plot_qc_report <- function(qc_metrics, output_dir = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先安装 ggplot2 包")
  }
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  plots <- list()
  
  # 特征数分布
  p1 <- ggplot2::ggplot(qc_metrics, ggplot2::aes(x = sample, y = nFeatures)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = "每个样本的特征数", x = "样本", y = "特征数")
  plots[["nFeatures"]] <- p1
  
  # 总计数分布
  p2 <- ggplot2::ggplot(qc_metrics, ggplot2::aes(x = sample, y = nCounts)) +
    ggplot2::geom_bar(stat = "identity", fill = "coral") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = "每个样本的总计数", x = "样本", y = "总计数")
  plots[["nCounts"]] <- p2
  
  # 零值百分比
  p3 <- ggplot2::ggplot(qc_metrics, ggplot2::aes(x = sample, y = percentZero)) +
    ggplot2::geom_bar(stat = "identity", fill = "mediumpurple") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = "每个样本的零值百分比", x = "样本", y = "零值百分比 (%)")
  plots[["percentZero"]] <- p3
  
  if (!is.null(output_dir)) {
    ggplot2::ggsave(file.path(output_dir, "01_nFeatures.pdf"), p1, width = 10, height = 6)
    ggplot2::ggsave(file.path(output_dir, "02_nCounts.pdf"), p2, width = 10, height = 6)
    ggplot2::ggsave(file.path(output_dir, "03_percentZero.pdf"), p3, width = 10, height = 6)
    message(paste("QC报告已保存到:", output_dir))
  }
  
  invisible(plots)
}

#' 过滤低质量样本
#'
#' @param counts 计数矩阵
#' @param qc_metrics QC指标数据框
#' @param min_features 最小特征数
#' @param min_counts 最小计数
#' @param max_percent_zero 最大零值百分比
#'
#' @return 过滤后的计数矩阵
#'
#' @examples
#' \dontrun{
#' filtered_counts <- filter_low_quality_samples(counts, qc_metrics, 
#'                                               min_features = 500)
#' }
#'
#' @export
filter_low_quality_samples <- function(counts, qc_metrics, 
                                        min_features = 500, 
                                        min_counts = 1000, 
                                        max_percent_zero = 90) {
  # 确定保留的样本
  keep <- (qc_metrics$nFeatures >= min_features) &
          (qc_metrics$nCounts >= min_counts) &
          (qc_metrics$percentZero <= max_percent_zero)
  
  n_removed <- sum(!keep)
  n_kept <- sum(keep)
  
  message(paste("移除", n_removed, "个样本，保留", n_kept, "个样本"))
  
  # 过滤计数矩阵
  filtered_counts <- counts[, keep]
  
  return(filtered_counts)
}

#' 离群值检测
#'
#' @param qc_metrics QC指标数据框
#' @param method 检测方法 ("mad", "iqr", "zscore")
#' @param threshold 阈值 (MAD/IQR方法用2.5, Z-score方法用3)
#'
#' @return 逻辑向量，TRUE表示离群值
#'
#' @examples
#' \dontrun{
#' outliers <- outlier_detection(qc_metrics, method = "mad")
#' }
#'
#' @export
outlier_detection <- function(qc_metrics, method = "mad", threshold = 2.5) {
  if (!(method %in% c("mad", "iqr", "zscore"))) {
    stop("method 必须为 'mad', 'iqr' 或 'zscore'")
  }
  
  nCounts <- qc_metrics$nCounts
  
  if (method == "mad") {
    median_val <- median(nCounts)
    mad_val <- stats::mad(nCounts)
    outliers <- abs(nCounts - median_val) > threshold * mad_val
  } else if (method == "iqr") {
    q1 <- stats::quantile(nCounts, 0.25)
    q3 <- stats::quantile(nCounts, 0.75)
    iqr <- q3 - q1
    outliers <- (nCounts < q1 - threshold * iqr) | (nCounts > q3 + threshold * iqr)
  } else if (method == "zscore") {
    z_scores <- abs(scale(nCounts))
    outliers <- z_scores > threshold
  }
  
  message(paste("检测到", sum(outliers), "个离群值"))
  return(outliers)
}