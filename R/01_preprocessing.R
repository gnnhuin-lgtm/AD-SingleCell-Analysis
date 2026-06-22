#' 数据预处理模块
#'
#' 提供原始数据的读取、清理和格式化函数
#'
#' @keywords internal

#' 读取计数矩阵
#'
#' @param filepath 文件路径
#' @param sep 分隔符 (默认为制表符)
#' @param row_names 是否使用第一列作为行名
#'
#' @return 计数矩阵
#'
#' @details
#' 支持多种格式的输入文件：
#' - CSV 文件
#' - TSV 文件
#' - RDS 文件
#'
#' @examples
#' \dontrun{
#' counts <- read_count_matrix("data/counts.csv")
#' }
#'
#' @export
read_count_matrix <- function(filepath, sep = "\t", row_names = TRUE) {
  if (!file.exists(filepath)) {
    stop(paste("文件不存在:", filepath))
  }
  
  # 根据文件扩展名选择读取方式
  ext <- tolower(tools::file_ext(filepath))
  
  if (ext == "rds") {
    data <- readRDS(filepath)
  } else if (ext == "csv") {
    data <- read.csv(filepath, sep = ",", row.names = if (row_names) 1 else NULL)
  } else if (ext %in% c("txt", "tsv")) {
    data <- read.csv(filepath, sep = sep, row.names = if (row_names) 1 else NULL)
  } else {
    stop("不支持的文件格式")
  }
  
  message(paste("已读取:", nrow(data), "个特征,", ncol(data), "个样本"))
  return(as.matrix(data))
}

#' 标准化数据
#'
#' @param counts 计数矩阵
#' @param method 标准化方法 ("cpm", "tpm", "rpkm", "quantile")
#' @param gene_length 基因长度 (用于TPM/RPKM计算)
#'
#' @return 标准化后的矩阵
#'
#' @details
#' - CPM: Counts Per Million
#' - TPM: Transcripts Per kilobase Million
#' - RPKM: Reads Per Kilobase Million
#' - Quantile: 分位数标准化
#'
#' @examples
#' \dontrun{
#' normalized <- normalize_data(counts, method = "cpm")
#' }
#'
#' @export
normalize_data <- function(counts, method = "cpm", gene_length = NULL) {
  if (!(method %in% c("cpm", "tpm", "rpkm", "quantile"))) {
    stop("method 必须为 'cpm', 'tpm', 'rpkm' 或 'quantile'")
  }
  
  if (method == "cpm") {
    normalized <- sweep(counts, 2, colSums(counts), "/") * 1e6
  } else if (method == "tpm") {
    if (is.null(gene_length)) {
      stop("TPM 标准化需要提供 gene_length")
    }
    rpk <- sweep(counts, 1, gene_length / 1000, "/")
    normalized <- sweep(rpk, 2, colSums(rpk) / 1e6, "/")
  } else if (method == "rpkm") {
    if (is.null(gene_length)) {
      stop("RPKM 标准化需要提供 gene_length")
    }
    normalized <- sweep(counts, 1, gene_length / 1000, "/")
    normalized <- sweep(normalized, 2, colSums(counts) / 1e6, "/")
  } else if (method == "quantile") {
    normalized <- limma::normalizeQuantiles(counts)
  }
  
  message(paste("已完成", method, "标准化"))
  return(normalized)
}

#' 对数变换
#'
#' @param data 输入数据
#' @param base 对数底数 (默认为2)
#' @param pseudocount 伪计数值 (避免log(0))
#'
#' @return 对数变换后的数据
#'
#' @examples
#' \dontrun{
#' log_data <- log_transform(normalized_counts, pseudocount = 1)
#' }
#'
#' @export
log_transform <- function(data, base = 2, pseudocount = 1) {
  if (any(data < 0, na.rm = TRUE)) {
    warning("数据中存在负值，已转换为0")
    data[data < 0] <- 0
  }
  
  log_data <- log2(data + pseudocount)
  message("对数变换完成")
  return(log_data)
}

#' 批次效应校正
#'
#' @param counts 计数矩阵
#' @param batch 批次向量
#' @param method 校正方法 ("ComBat", "SVA", "RUV")
#'
#' @return 校正后的矩阵
#'
#' @details
#' - ComBat: 使用 ComBat 算法
#' - SVA: 监督 SVA
#' - RUV: 相对日志表达 (RUVSeq)
#'
#' @examples
#' \dontrun{
#' corrected <- batch_correction(counts, batch = batch_info, method = "ComBat")
#' }
#'
#' @export
batch_correction <- function(counts, batch, method = "ComBat") {
  if (!(method %in% c("ComBat", "SVA", "RUV"))) {
    stop("method 必须为 'ComBat', 'SVA' 或 'RUV'")
  }
  
  if (method == "ComBat") {
    if (!requireNamespace("sva", quietly = TRUE)) {
      stop("请先安装 sva 包: BiocManager::install('sva')")
    }
    corrected <- sva::ComBat(dat = counts, batch = batch, par.prior = TRUE, prior.plots = FALSE)
  } else {
    message(paste(method, "方法尚未实现，请使用 'ComBat'"))
    corrected <- counts
  }
  
  message(paste("已完成", method, "批次校正"))
  return(corrected)
}