library(ggplot2)
library(patchwork)
library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)

# 确保我们在 SCT 之前，先处理 combined_raw
# 1. 拆分图层，模拟每个 Plate 是一个独立 Sample 的逻辑
# 这样你的循环就能跑起来，且每个 Plate 都有自己的图表
plate_list <- SplitObject(combined_raw, split.by = "Plate_ID")

filtered_seurat_objects <- list()

# ===== 🪭 2. QC 可视化循环 =====
for (plate_id in names(plate_list)) {
  obj <- plate_list[[plate_id]]
  
  # A. 计算基础指标 (V5 兼容)
  obj$nCount_RNA <- colSums(GetAssayData(obj, layer = "counts"))
  obj$nFeature_RNA <- colSums(GetAssayData(obj, layer = "counts") > 0)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # B. 预检查过滤阈值
  low_feature_cutoff <- 100
  high_feature_cutoff <- 7000
  mt_cutoff <- 25
  low_umi_cutoff <- 500
  
  pass_qc <- which(obj$nFeature_RNA > low_feature_cutoff & 
                     obj$nFeature_RNA < high_feature_cutoff & 
                     obj$nCount_RNA > low_umi_cutoff & 
                     obj$percent.mt < mt_cutoff)
  
  print(paste0("Plate: ", plate_id, " | Initial: ", ncol(obj), " | Pass QC: ", length(pass_qc)))
  
  # --- 安全检查点 1：细胞太少则跳过 ---
  if (length(pass_qc) < 20) {
    message(paste0("⚠️ Skipping ", plate_id, ": Too few cells pass QC."))
    next 
  }
  
  # 执行初步过滤
  obj_filtered <- obj[, pass_qc]
  
  # C. 可视化 (完整绘图逻辑)
  # 1. Feature vs. UMI 散点图
  p_scatter_feature_umi <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                                          pt.size = 0.8, cols = "#6F00FF") +
    geom_hline(yintercept = low_feature_cutoff, linetype = "dashed", color = "#6B3F69", size = 0.8) +
    geom_hline(yintercept = high_feature_cutoff, linetype = "dashed", color = "red", size = 0.8) +
    geom_vline(xintercept = low_umi_cutoff, linetype = "dotted", color = "#E45A92", size = 0.8) +
    ggtitle("Feature vs. UMI") + theme_minimal()
  
  # 2. MT% vs. UMI 散点图
  p_scatter_mt_umi <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt", 
                                     pt.size = 0.8, cols = "#6F00FF") +
    geom_hline(yintercept = mt_cutoff, linetype = "dashed", color = "#3B0270", size = 0.8) +
    ggtitle("MT% vs. UMI") + theme_minimal()
  
  # 3. nFeature 小提琴图
  p_vln_feature <- VlnPlot(obj, features = "nFeature_RNA", pt.size = 0.1, cols = "#6F00FF") +
    geom_hline(yintercept = low_feature_cutoff, linetype = "dashed", color = "#6B3F69") +
    ggtitle("nFeature Dist.") + NoLegend()
  
  # 4. percent.mt 小提琴图
  p_vln_mt <- VlnPlot(obj, features = "percent.mt", pt.size = 0.1, cols = "#6F00FF") +
    geom_hline(yintercept = mt_cutoff, linetype = "dashed", color = "#3B0270") +
    ggtitle("MT% Dist.") + NoLegend()
  
  # 组合图表并保存
  combined_qc_plots_final <- (p_scatter_feature_umi | p_vln_feature) / (p_scatter_mt_umi | p_vln_mt) +
    plot_annotation(title = paste0("QC for Plate: ", plate_id))
  
  ggsave(paste0("QC/QC_Plate_", plate_id, ".png"), combined_qc_plots_final, width = 12, height = 10)
  
  # D. Doublet 检测
  if (ncol(obj_filtered) > 50) {
    # 临时标准化以通过 scDblFinder 需要的聚类步骤
    obj_tmp <- NormalizeData(obj_filtered, verbose = FALSE)
    sce <- as.SingleCellExperiment(obj_tmp)
    
    # 运行 scDblFinder
    sce <- scDblFinder(sce, verbose = FALSE)
    
    # 回传结果到 Seurat 对象
    obj_final <- as.Seurat(sce)
    # 只保留单细胞 (Singlets)
    obj_final <- subset(obj_final, subset = scDblFinder.class == "singlet")
  } else {
    obj_final <- obj_filtered
  }
  
  # 存储处理好的板子
  filtered_seurat_objects[[plate_id]] <- obj_final
  
  # 清理内存
  rm(obj_tmp, sce); gc()
}

# ===== 🏁 3. 最终合并（V5 稳健版） =====
if (length(filtered_seurat_objects) > 0) {
  
  combined_qc <- merge(
    x = filtered_seurat_objects[[1]], 
    y = filtered_seurat_objects[-1], 
    add.cell.ids = names(filtered_seurat_objects)
  )
  
  DefaultAssay(combined_qc) <- "RNA"
  tryCatch({
    combined_qc <- JoinLayers(combined_qc)
    message("Layers joined successfully.")
  }, error = function(e) {
    message("Layers were already joined or don't need joining.")
  })
  
  print(paste0("✅ 全部处理完成！最终剩余细胞总数: ", ncol(combined_qc)))
  
} else {
  stop("❌ 严重错误：没有细胞通过 QC 流程！")
}


# ==== Sanity check ====
as.data.frame(table(combined_qc$Plate_ID))


# ===== CONSTRUCTION ZONE =====
# ==== 1. SCTransform ====
# 之前的操作只完成了物理合并，现在我们要进行跨板子的标准化
# vars.to.regress 排除掉线粒体干扰，vst.flavor = "v2" 提高灵敏度
combined_qc <- SCTransform(combined_qc, 
                           vars.to.regress = "percent.mt", 
                           vst.flavor = "v2", 
                           verbose = TRUE)

# ==== 2. Reduc & batch corr ====
DefaultAssay(combined_qc) <- "SCT"

# 跑 PCA (主成分分析)
combined_qc <- RunPCA(combined_qc, npcs = 50, verbose = FALSE)

# 【关键】Seurat V5 独有的 split 逻辑：按 Plate 分层以准备集成
# 这一步不会增加内存负担，只是在 SCT Assay 内部标记好批次来源
combined_qc[["SCT"]] <- split(combined_qc[["SCT"]], f = combined_qc$Plate_ID)

# ==============================================================================
# ==== 3. 执行集成 (IntegrateLayers) ====
# 目标：抹平 43 个板子间的批次差异，让 T1 和 T2 的生物学信号凸显出来
combined_qc <- IntegrateLayers(
  object = combined_qc, 
  method = RPCAIntegration,  
  orig.reduction = "pca", 
  new.reduction = "rpca",
  assay = "SCT",
  verbose = TRUE
)

# ==== 4. 聚类与生成 UMAP 坐标 ====

# 基于校正后的空间 (integrated.rpca) 找邻居
combined_qc <- FindNeighbors(combined_qc, reduction = "integrated.rpca", dims = 1:30)

# 聚类：寻找细胞群
combined_qc <- FindClusters(combined_qc, resolution = 0.4)

# 生成可视化坐标
combined_qc <- RunUMAP(combined_qc, reduction = "integrated.rpca", dims = 1:30, 
                       reduction.name = "umap.rpca")

# ==== 5. 最终可视化：T1 (Monocytes) vs T2 (Macrophages) ====
# 看看你的细胞在 7 天里经历了什么
p_timepoint <- DimPlot(combined_qc, 
                       reduction = "umap.rpca", 
                       group.by = "Timepoint", 
                       cols = c("T1" = "#00BFC4", "T2" = "#F8766D"), # 经典配色
                       pt.size = 0.6) + 
  ggtitle("Single-cell Map: T1 (D0) vs T2 (D7)") +
  theme_minimal()

print(p_timepoint)