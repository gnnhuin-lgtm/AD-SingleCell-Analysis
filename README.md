# 生信技能 (Bioinformatics Skills)

一个专业的R语言生信分析框架，支持基因组数据处理（DNA-seq、RNA-seq）和蛋白质组学分析。

## 📋 项目概述

本项目提供了完整的生信分析工作流，包括：
- **数据预处理** - 原始数据的读取、清理和格式化
- **质量控制** - QC指标计算和质量评估
- **统计分析** - 差异表达分析、富集分析等
- **数据可视化** - 高质量图表生成

## 🎯 主要应用方向

### 1. 基因组数据处理
- **DNA-seq分析** - 变异检测、CNV分析
- **RNA-seq分析** - 转录组分析、差异表达、剪接体分析

### 2. 蛋白质组学
- 蛋白质定量分析
- 翻译后修饰（PTM）分析
- 蛋白质互作网络

## 📁 项目结构

```
├── README.md
├── .gitignore
├── DESCRIPTION              # R包元数据
├── NAMESPACE
├── R/                       # R函数源代码
│   ├── 01_preprocessing.R        # 数据预处理函数
│   ├── 02_quality_control.R      # QC函数
│   ├── 03_statistical_analysis.R # 统计分析函数
│   ├── 04_visualization.R        # 可视化函数
│   └── utils.R                   # 工具函数
├── man/                     # 函数文档
├── data/                    # 示例数据集
├── vignettes/               # 详细教程和案例
│   ├── 01_DNA-seq分析指南.Rmd
│   ├── 02_RNA-seq分析指南.Rmd
│   └── 03_蛋白质组学分析指南.Rmd
├── tests/                   # 单元测试
└── inst/
    └── scripts/             # 示例分析脚本
```

## 🚀 快速开始

### 安装依赖

```r
# 安装必要的R包
packages <- c(
  "tidyverse",      # 数据处理
  "ggplot2",        # 可视化
  "DESeq2",         # RNA-seq差异表达分析
  "limma",          # 蛋白质组和芯片数据分析
  "edgeR",          # 边际分析
  "GenomicRanges",  # 基因组数据处理
  "VariantAnnotation", # 变异注释
  "ComplexHeatmap"  # 复杂热图
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("DESeq2", "limma", "GenomicRanges", "VariantAnnotation", "ComplexHeatmap")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}
```

### 加载本包

```r
devtools::load_all()
# 或者
library(bioinformatics.skills)
```

## 📚 主要功能模块

### 模块一：数据预处理 (01_preprocessing.R)
- `read_count_matrix()` - 读取计数矩阵
- `normalize_data()` - 数据标准化
- `log_transform()` - 对数变换
- `batch_correction()` - 批次效应校正

### 模块二：质量控制 (02_quality_control.R)
- `calculate_qc_metrics()` - 计算QC指标
- `plot_qc_report()` - QC可视化报告
- `filter_low_quality_samples()` - 样本过滤
- `outlier_detection()` - 离群值检测

### 模块三：统计分析 (03_statistical_analysis.R)
- `differential_expression()` - 差异表达分析
- `enrichment_analysis()` - GO/KEGG富集分析
- `correlation_analysis()` - 相关性分析
- `survival_analysis()` - 生存分析

### 模块四：数据可视化 (04_visualization.R)
- `plot_volcano()` - 火山图
- `plot_ma()` - MA图
- `plot_heatmap()` - 热图
- `plot_pca()` - PCA图
- `plot_boxplot()` - 箱线图
- `plot_network()` - 网络图

## 📖 使用教程

详细的分析教程请参考 `vignettes/` 目录：

1. **DNA-seq分析指南** - 从原始数据到变异注释的完整流程
2. **RNA-seq分析指南** - 转录组分析的标准工作流
3. **蛋白质组学分析指南** - 蛋白质组数据的处理和分析

## 💾 示例数据

在 `data/` 目录中包含示例数据集，可用于测试和学习：
- `example_counts.rds` - 示例转录组计数矩阵
- `example_metadata.csv` - 示例样本信息
- `example_protein.csv` - 示例蛋白质组数据

## 🧪 单元测试

运行测试确保代码质量：

```r
devtools::test()
```

## 📝 代码风格指南

- 使用 `snake_case` 命名函数
- 添加详细的函数文档注释
- 遵循 Tidyverse 代码风格
- 包含错误处理和输入验证

## 🤝 贡献指南

欢迎提交问题报告和改进建议！

## 📄 许可证

MIT License

## 📧 联系方式

如有问题，欢迎提交 Issue 或 PR。

---

**最后更新**: 2026-06-22