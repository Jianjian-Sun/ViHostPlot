# ViHostPlot 0.0.4

- **部分重构**
  - 重构了 circos 图中宿主染色体信息的读取逻辑，现在支持通过 `host` 参数指定常见宿主，也可以通过 `chrom_file` 传入自定义染色体长度文件。
  - 将宿主染色体相关函数独立到 `host_chrom_utils.R` 中，使 circos 绘图逻辑和宿主基因组解析逻辑分离。
  - 更改 circos 图中病毒序列的处理逻辑，现在要求用户显式传入 `virus_name` 和 `virus_length`，而不是在染色体文件中手动添加病毒序列。

- **新功能与优化**
  - 新增通过 UCSC 在线获取宿主染色体长度的功能，例如 `human` 可自动解析为 `hg38`。

- **移除与简化**
  - 移除了 `color_file` 和 `palettes` 参数。
  - `chrom_file` 现在只需要包含宿主染色体信息，不再需要包含病毒序列。
  - 染色体名称会统一去除 `chr` 前缀，例如 `chr1` 会被标准化为 `1`。

# ViHostPlot 0.0.3

* **部分重构**
    * 修复了 oncoplot 图在绘制图时有时候标签的标题不会显示的问题
    * 更改 circlize 中寻找病毒序列的逻辑，现在要求用户显式传入病毒名，而不是函数内部自己去寻找最短的序列作为病毒名
    * 对于 circlize 添加了逻辑`histgram`的逻辑，在散点层可续选择绘制直方图或者散点图来代表测序深度。


# ViHostPlot 0.0.2

* **部分重构**
    * 跳过生成中间 VCF 文件的步骤。现在 `visualize_viral_integration()` 通过新增的 `create_gi_from_table()` 函数，直接将输入的表格数据解析为 `GInteractions` 对象。
    * 将核心绘图引擎迁移至 `circlize` 的基因组专用函数族（包括 `circos.genomicInitialize`、`circos.genomicTrackPlotRegion`、`circos.genomicPoints`、`circos.genomicLink` 等），以实现更稳健、精确的基因组坐标处理。

* **精简依赖**
    * 移除了重量级的 Bioconductor 依赖包：`VariantAnnotation`、`SummarizedExperiment`。
    * 移除了 Tidyverse 依赖包：`dplyr`、`stringr`。

* **新功能与优化**
    * 新增 `create_host_virus_layout()` 函数，帮助用户快速生成标准的轨道布局列表。
    * 新增 `match_chr_style()` 内部辅助函数，可智能处理并统一宿主与病毒染色体命名风格的不一致（例如自动匹配 "chr1" 和 "1"）。
    * 新增 `make_axis_labels()` 内部辅助函数，可自动为 Ideogram 轨道格式化坐标轴单位（宿主染色体使用 Mb，病毒序列使用 kb）。
    * 新增输入验证函数 `validate_layout_list()`，用于防止错误或畸形的轨道布局定义导致绘图失败。


# ViHostPlot 0.0.1

* 重构 `oncoplot()` 布局逻辑并增强临床轨道功能
* 移除 `patchwork` / `ggplotify` 依赖，采用纯 `aplot` 原生布局
* 修改函数标签分类机制，移除用户指定 `discrete` 和 `continuous`，使用 `binning_numeric()` 函数，支持连续变量的自动分箱与标签美化
* 暴露 `n_breaks` 和 `colors` 参数以提升可视化自定义程度
