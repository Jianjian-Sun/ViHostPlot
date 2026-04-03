# ViHostPlot 0.0.1

* 重构 `oncoplot()` 布局逻辑并增强临床轨道功能
* 移除 `patchwork` / `ggplotify` 依赖，采用纯 `aplot` 原生布局
* 修改函数标签分类机制，移除用户指定 `discrete` 和 `continuous`，使用 `binning_numeric()` 函数，支持连续变量的自动分箱与标签美化
* 暴露 `n_breaks` 和 `colors` 参数以提升可视化自定义程度
