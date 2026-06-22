#' 工具函数集合
#'
#' 提供通用的数据处理和辅助函数
#'
#' @keywords internal

#' 检查输入数据
#'
#' @param data 输入数据框或矩阵
#' @param name 数据名称
#'
#' @return 检查通过返回TRUE，失败抛出错误
#'
#' @examples
#' \dontrun{
#' check_data(mtcars, "示例数据")
#' }
#'
check_data <- function(data, name = "data") {
  if (is.null(data)) {
    stop(paste(name, "不能为NULL"))
  }
  if (nrow(data) == 0) {
    stop(paste(name, "行数不能为0"))
  }
  if (ncol(data) == 0) {
    stop(paste(name, "列数不能为0"))
  }
  return(TRUE)
}

#' 获取数据摘要统计信息
#'
#' @param data 输入数据框或矩阵
#'
#' @return 包含摘要信息的列表
#'
#' @examples
#' \dontrun{
#' get_summary_stats(mtcars)
#' }
#'
get_summary_stats <- function(data) {
  check_data(data)
  
  list(
    n_samples = nrow(data),
    n_features = ncol(data),
    mean_value = mean(as.matrix(data), na.rm = TRUE),
    sd_value = sd(as.numeric(as.matrix(data)), na.rm = TRUE),
    min_value = min(as.matrix(data), na.rm = TRUE),
    max_value = max(as.matrix(data), na.rm = TRUE)
  )
}

#' 标准化矩阵
#'
#' @param matrix 输入矩阵
#' @param method 标准化方法 ("zscore", "minmax", "quantile")
#'
#' @return 标准化后的矩阵
#'
#' @examples
#' \dontrun{
#' normalized <- scale_matrix(expr_matrix, method = "zscore")
#' }
#'
scale_matrix <- function(matrix, method = "zscore") {
  if (!(method %in% c("zscore", "minmax", "quantile"))) {
    stop("method 必须为 'zscore', 'minmax' 或 'quantile'")
  }
  
  switch(method,
    "zscore" = scale(matrix),
    "minmax" = {
      minmax <- function(x) {
        (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
      }
      apply(matrix, 2, minmax)
    },
    "quantile" = {
      quantile_norm <- function(x) {
        rank(x) / length(x)
      }
      apply(matrix, 2, quantile_norm)
    }
  )
}

#' 保存分析结果
#'
#' @param result 结果对象
#' @param filepath 输出文件路径
#' @param format 输出格式 ("csv", "rds", "xlsx")
#'
#' @return 隐形返回保存的文件路径
#'
#' @examples
#' \dontrun{
#' save_result(analysis_result, "./results/output.csv", format = "csv")
#' }
#'
save_result <- function(result, filepath, format = "csv") {
  dir.create(dirname(filepath), showWarnings = FALSE, recursive = TRUE)
  
  switch(format,
    "csv" = write.csv(result, filepath, row.names = FALSE),
    "rds" = saveRDS(result, filepath),
    "xlsx" = {
      if (!requireNamespace("openxlsx", quietly = TRUE)) {
        install.packages("openxlsx")
      }
      openxlsx::write.xlsx(result, filepath)
    }
  )
  
  message(paste("结果已保存到:", filepath))
  invisible(filepath)
}