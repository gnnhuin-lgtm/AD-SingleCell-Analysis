# AD-SingleCell-Analysis

## 项目概述

**AD-SingleCell-Analysis** 是一个针对阿尔茨海默病 (Alzheimer's Disease, AD) 的单细胞转录组分析项目。该项目提供了从原始数据处理到深度分析和可视化的完整工作流。

### 🎯 项目目标

- 对 AD 相关的单细胞测序数据进行质量控制和预处理
- 细胞聚类和注释分析
- 基因表达模式挖掘
- 病变和正常细胞对比分析
- 关键基因和通路鉴定

---

## 📋 目录结构

```
AD-SingleCell-Analysis/
├── .github/
│   └── workflows/                 # GitHub Actions 工作流
│       ├── tests.yml              # 自动化测试
│       └── lint.yml               # 代码质量检查
├── src/
│   ├── __init__.py
│   ├── preprocessing.py           # 数据预处理
│   ├── analysis.py                # 主要分析
│   ├── visualization.py           # 数据可视化
│   └── utils.py                   # 工具函数
├── notebooks/
│   ├── 01_data_loading.ipynb      # 数据加载和QC
│   ├── 02_preprocessing.ipynb     # 数据预处理
│   ├── 03_clustering.ipynb        # 聚类分析
│   ├── 04_annotation.ipynb        # 细胞注释
│   └── 05_comparative_analysis.ipynb  # 对比分析
├── tests/
│   ├── __init__.py
│   └── test_preprocessing.py      # 单元测试
├── data/
│   ├── raw/                       # 原始数据（gitignored）
│   └── processed/                 # 处理后的数据（gitignored）
├── results/
│   ├── figures/                   # 输出图表
│   ├── tables/                    # 输出表格
│   └── reports/                   # 分析报告
├── config/
│   └── config.yaml                # 项目配置文件
├── requirements.txt               # Python 依赖
├── setup.py                       # 包安装配置
├── .gitignore                     # Git 忽略文件
├── README.md                      # 项目说明（本文件）
└── LICENSE                        # 开源许可证
```

---

## 🚀 快速开始

### 前置要求

- Python 3.8+
- pip 或 conda

### 安装

1. **克隆仓库**
```bash
git clone https://github.com/gnnhuin-lgtm/AD-SingleCell-Analysis.git
cd AD-SingleCell-Analysis
```

2. **创建虚拟环境（推荐）**
```bash
# 使用 venv
python -m venv venv
source venv/bin/activate  # macOS/Linux
# 或
venv\Scripts\activate     # Windows

# 或使用 conda
conda create -n ad-analysis python=3.9
conda activate ad-analysis
```

3. **安装依赖**
```bash
pip install -r requirements.txt
```

4. **安装项目包**
```bash
pip install -e .
```

---

## 📊 使用方法

### 数据准备

将您的单细胞转录组数据放在 `data/raw/` 目录下。支持的格式：
- `.h5ad` (AnnData 格式)
- `.h5` (HDF5 格式)
- `.csv` / `.tsv`

### 运行分析

#### 方法 1：使用 Jupyter Notebooks（推荐新手）

```bash
jupyter lab notebooks/
```

按顺序运行以下笔记本：
1. `01_data_loading.ipynb` - 加载和初步质量控制
2. `02_preprocessing.ipynb` - 数据预处理和标准化
3. `03_clustering.ipynb` - 无监督聚类
4. `04_annotation.ipynb` - 细胞类型注释
5. `05_comparative_analysis.ipynb` - AD vs 正常对比

#### 方法 2：使用 Python 脚本

```python
from src.preprocessing import QualityControl, Normalize
from src.analysis import Clustering, Annotation
import anndata as ad

# 加载数据
adata = ad.read_h5ad('data/raw/your_data.h5ad')

# 质量控制
qc = QualityControl(adata)
adata = qc.filter_cells_genes()

# 标准化
normalizer = Normalize(adata)
adata = normalizer.normalize_and_log_transform()

# 聚类
clustering = Clustering(adata)
adata = clustering.run_clustering()

# 保存结果
adata.write('data/processed/analysis_results.h5ad')
```

---

## 📁 核心模块说明

### `src/preprocessing.py`
- **QualityControl**: 细胞和基因过滤
- **Normalize**: 数据标准化和对数转换
- **VariableGeneSelection**: 高变基因筛选

### `src/analysis.py`
- **Clustering**: PCA、UMAP、聚类分析
- **DifferentialExpression**: 差异表达分析

### `src/visualization.py`
- **Visualizer**: 包含 UMAP、小提琴图、热力图等可视化函数

---

## 🧪 测试

运行单元测试：

```bash
pytest tests/ -v
```

生成覆盖率报告：

```bash
pytest tests/ --cov=src --cov-report=html
```

---

## 📈 分析流程概览

```
原始数据
    ↓
质量控制 (QC)
    ↓
数据预处理 (归一化、对数变换)
    ↓
高变基因筛选
    ↓
降维分析 (PCA, UMAP)
    ↓
聚类分析
    ↓
细胞类型注释
    ↓
差异表达分析 (AD vs Control)
    ↓
功能富集分析
    ↓
可视化和报告
```

---

## 📚 关键参数

在 `config/config.yaml` 中可以配置：

```yaml
# 质量控制参数
quality_control:
  min_genes: 200
  max_genes: 10000
  min_cells: 3
  percent_mt: 20

# 聚类参数
clustering:
  method: "leiden"
  resolution: 0.5
  n_neighbors: 15

# 差异表达分析
differential_expression:
  method: "wilcoxon"
  min_log2fc: 0.25
  max_pval: 0.05
```

---

## 🔍 示例数据

您可以使用以下公开数据集进行测试：
- [Allen Brain Atlas](https://portal.brain-map.org/)
- [GEO Omnibus](https://www.ncbi.nlm.nih.gov/geo/)
- [10x Genomics](https://www.10xgenomics.com/datasets)

---

## 📝 引用

如果您在研究中使用了本项目，请引用：

```bibtex
@software{ad_singlecell_2024,
  title={AD-SingleCell-Analysis},
  author={gnnhuin-lgtm},
  year={2024},
  url={https://github.com/gnnhuin-lgtm/AD-SingleCell-Analysis}
}
```

---

## 🤝 贡献指南

欢迎提交 Issue 和 Pull Request！

### 提交步骤：
1. Fork 此仓库
2. 创建功能分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 打开 Pull Request

### 代码规范：
- 遵循 PEP 8 标准
- 添加相应的单元测试
- 更新文档

---

## 📄 许可证

本项目采用 MIT 许可证，详见 [LICENSE](LICENSE) 文件。

---

## 📧 联系方式

- GitHub: [@gnnhuin-lgtm](https://github.com/gnnhuin-lgtm)
- 有问题？请提交 [Issue](https://github.com/gnnhuin-lgtm/AD-SingleCell-Analysis/issues)

---

## 📖 参考资源

- [Scanpy 文档](https://scanpy.readthedocs.io/)
- [AnnData 文档](https://anndata.readthedocs.io/)
- [Seurat 单细胞分析指南](https://satijalab.org/seurat/)
- [Single-cell best practices](https://www.sc-best-practices.org/)

---

**最后更新**: 2024年6月  
**维护者**: gnnhuin-lgtm
