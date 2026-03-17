# ====== 🫐 INFO ======
  #
  # This script is (attemping) analyzing the scRNA-seq data of trained monocytes.
  # Based on previous script, this one now further sub-divides the sample and only
  # Uses monocytes.
  #
  # Original paper: https://www.jci.org/articles/view/147719
  #
  # Original database: https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-9702/sdrf 
  # 
  # Original code: https://github.com/CiiM-Bioinformatics-group/trained_immunity 
  #
  
  
# ==== 🛟 Questions to be solved ====
"Myeloid and PBMC separation from metadata - remove all PBMCs"
"Play with clustering parameters to see if we can get better result"
"Implement Helen He QC"

# =====🍇 Library ====
library(Seurat)
library(dplyr)
library(harmony)
library(readr)
library(ggplot2)



# === Process data ===
setwd("C:/Users/denxiack/Documents/Bowen_zhang/Bowen_Zhang")
files <- list.files(pattern = "\\.ReadCounts\\.tsv$")
files

csv_test <- read_tsv("RMC-SM-003_HM3JFBGX9_S5_R2.ReadCounts.tsv")
tail(csv_test$GENEID)

metadata_bowen <- read_csv("Metadata_TrainedImmunity.csv")
summary(metadata_bowen)
View(metadata_bowen)
metadata_bowen$`Cell type`  # Get all yourmonocytes

# Select only monocytes
process_plate <- function(i) {
  row <- metadata_bowen[i, ]
  if (row$`Cell type` != "Monocytes") {
    message(paste("Skipping PBMC plate:", row$Plate_ID))
    return(NULL)
  }
  
  plate_id <- row$Plate_ID
  files <- list.files(pattern = "\\.ReadCounts\\.tsv$")
  file_path <- files[grep(plate_id, files)]
  
  if(length(file_path) == 0) return(NULL)
  
  # Read and trim like before
  counts <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1)
  rownames(counts) <- gsub("__chr.*", "", rownames(counts))
  colnames(counts) <- paste0(plate_id, "_", colnames(counts))
  
  obj <- CreateSeuratObject(counts = counts, project = "TrainedMonocytes")
  obj$stimulus <- row$`Training stimulus`
  obj$donor    <- row$Donor
  obj$timepoint <- row$Timepoint
  obj$origin_cell_type <- row$`Cell type`
  
  return(obj)
}

obj_list_mono <- lapply(1:nrow(metadata_bowen), process_plate)
obj_list_mono <- Filter(Negate(is.null), obj_list_mono) # PBMCs = NULL here

# Make monocyte seurat obj
combined_seu_mono <- merge(obj_list_mono[[1]], 
                           y = obj_list_mono[2:length(obj_list_mono)], 
                           add.cell.ids = NULL)
# QC
table(combined_seu_mono$origin_cell_type) # only monocytes now


# ==== Clustering ====
combined_seu_mono <- SCTransform(combined_seu_mono, verbose = FALSE)
combined_seu_mono <- RunPCA(combined_seu_mono, npc = 30, verbose = FALSE)

# correct for donor by harmony
combined_seu_mono <- RunHarmony(combined_seu_mono, group.by.vars = "donor")

combined_seu_mono <- RunUMAP(combined_seu_mono, reduction = "harmony", dims = 1:20)

# Increase distance...
#combined_seu_mono <- RunUMAP(combined_seu_filtered, 
#        reduction = "harmony", 
#        dims = 1:20, 
#        min.dist = 0.5,  
#        spread = 1.5)

names(combined_seu_mono@graphs) # "SCT_snn" = old harmony, no spread; "SCT_nn" =  spread

combined_seu_mono <- FindNeighbors(combined_seu_mono, reduction = "harmony", dims = 1:20)
combined_seu_mono <- FindClusters(combined_seu_mono, graph.name = "SCT_snn", resolution = 0.8) # higher since its only monocyte now

# Fact check: T / B cells?
#VlnPlot(combined_seu_mono, features = c("CD14", "CD68", "CD3E", "MS4A1"), stack = TRUE)

# ==== QC Filtering ====
# Finding out mitochondrial genes
grep("^MT", rownames(combined_seu_mono), value = TRUE)
combined_seu_mono[["percent.mt"]] <- PercentageFeatureSet(combined_seu_mono, pattern = "MT-")
combined_seu_mono <- JoinLayers(combined_seu_mono)

# For plate mode specifically: Merge plates 
combined_seu_mono <- JoinLayers(combined_seu_mono)
Layers(combined_seu_mono)

# Manual cehck ---
qc_data <- combined_seu_mono@meta.data
qc_data[is.na(qc_data)] <- 0

p1 <- ggplot(qc_data, aes(x = "Monocytes", y = nFeature_RNA)) + 
  geom_violin(fill = "#A6CEE3") + 
  geom_jitter(width = 0.2, size = 0.1, alpha = 0.1) +
  theme_classic() + labs(x = NULL, title = "nFeature_RNA")

p2 <- ggplot(qc_data, aes(x = "Monocytes", y = nCount_RNA)) + 
  geom_violin(fill = "#B2DF8A") + 
  geom_jitter(width = 0.2, size = 0.1, alpha = 0.1) +
  theme_classic() + labs(x = NULL, title = "nCount_RNA")

p3 <- ggplot(qc_data, aes(x = "Monocytes", y = percent.mt)) + 
  geom_violin(fill = "#FB9A99") + 
  geom_jitter(width = 0.2, size = 0.1, alpha = 0.1) +
  theme_classic() + labs(x = NULL, title = "Percent MT")

p1 + p2 + p3 # No ideal, go implement your Helen QC

# Decide on cut offs:
combined_seu_filtered <- subset(combined_seu_mono, 
                                subset = nFeature_RNA > 500 & 
                                  nCount_RNA > 1000 & 
                                  percent.mt < 15)

summary(combined_seu_mono$percent.mt)


# ---- Monocyte Fix: Join layers  ----
combined_seu_mono <- JoinLayers(combined_seu_mono)

combined_seu_mono[["percent.mt"]] <- PercentageFeatureSet(
  combined_seu_mono, 
  pattern = "^MT-", 
  assay = "RNA"
)

combined_seu_mono$nFeature_RNA <- colSums(GetAssayData(combined_seu_mono, layer = "counts") > 0)
combined_seu_mono$nCount_RNA   <- colSums(GetAssayData(combined_seu_mono, layer = "counts"))

summary(combined_seu_mono$percent.mt)
summary(combined_seu_mono$nFeature_RNA)

combined_seu_filtered <- subset(
  combined_seu_mono, 
  subset = nFeature_RNA > 2000 &   # from 500 -> 2000, much better output 
    nCount_RNA > 1000 & 
    percent.mt < 15
)

print(paste("过滤后存活细胞数:", ncol(combined_seu_filtered)))
# Debug:
print(paste("总细胞数:", ncol(combined_seu_mono)))
print(paste("nFeature > 500 剩余:", sum(combined_seu_mono$nFeature_RNA > 2500, na.rm = TRUE)))
print(paste("nCount > 1000 剩余:", sum(combined_seu_mono$nCount_RNA > 1000, na.rm = TRUE)))
print(paste("percent.mt < 15 剩余:", sum(combined_seu_mono$percent.mt < 15, na.rm = TRUE)))
# % MT is issue
# This is a strict parameter. Only ~ 3000 cells left....
# ==== Clustering visualization ====
combined_seu_filtered$percent.mt
p1 <- DimPlot(combined_seu_filtered, reduction = "umap", group.by = "stimulus",
              pt.size = 0.5)
p2 <- DimPlot(combined_seu_filtered, reduction = "umap", group.by = "donor",
              pt.size = 0.5)
p3 <- DimPlot(combined_seu_filtered, reduction = "umap", group.by = "seurat_clusters",
              pt.size = 0.5)
p4 <- DimPlot(combined_seu_filtered, reduction = "umap", group.by = "timepoint",
              pt.size = 0.5)
p5 <- DimPlot(combined_seu_filtered, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.5)

FeaturePlot(combined_seu_filtered, features = "nFeature_RNA") # scattering shows the outsiders are low features

p_sct_sum <- p1 + p2 + p3 + p4 + p5
print(p_sct_sum)
ggsave("bowen_zhang_myeloid_sct_res08.png", p_sct_sum, height = 8, width = 16, dpi = 300)
getwd()

# saveRDS(combined_seu_filtered, "monocyte_df_harm_08_clustering.rds")

#### work with this version that has 10 clusters.
#==== Azimuth  ====
combined_seu_filtered <- RunAzimuth(
  query = combined_seu_filtered, 
  reference = "pbmcref")

head(combined_seu_filtered@meta.data) #Works

dp <- DimPlot(combined_seu_filtered, 
              reduction = "umap",  
              group.by = "predicted.celltype.l3", 
              label = TRUE, 
              repel = TRUE, 
              pt.size = 1.5,
              raster=TRUE) + 
  ggtitle("Annotated Azimuth Bowen")

dp <- DimPlot_scCustom(
  seurat_object = combined_seu_filtered, 
  reduction = "umap", 
  group.by = "predicted.celltype.l3", 
  label = TRUE,          
  repel = TRUE,
  pt.size = 1,        
  label.size = 3,
  colors_use = DiscretePalette_scCustomize(num_colors = 13,
                                           palette = "varibow")
) + 
  ggtitle("Annotated Azimuth Bowen (L3)") +
  theme(plot.title = element_text(hjust = 0.5))

print(dp)
ggsave("Bowen_myeloid_10Cluster_Azimuth.png", height = 8, width = 9, dpi = 300)


# ==== Check the monocytes ====
mono_genes <- c(
  "CD14", "FCGR3A", "IL1B",
  "PTGS2", "CXCL10",
  "CXCL11", "CD68",
  "CD83"
)
summary(rownames(combined_seu_filtered) == "FCGR3A")  # CD16 does not exist in the panel...?

mono_genes_present <- intersect(mono_genes, rownames(combined_seu_filtered))

# 绘图检查分布
mono_check <- FeaturePlot_scCustom(
  seurat_object = combined_seu_filtered, 
  features = mono_genes_present, 
  reduction = "umap", 
  pt.size = 0.8,
  order = FALSE  
)

mono_check
ggsave("mono_genes_on_10clusters_myeloid.png", mono_check, height = 15, width = 15, dpi = 300)

# ==== Monocytes by timepoingt ====
mono_genes <- c("CD14", "FCGR3A", "IL1B", "PTGS2", "CXCL10", "CXCL11", "CD68", "CD83")
mono_genes_present <- intersect(mono_genes, rownames(combined_seu_filtered))

mono_split_check <- FeaturePlot_scCustom(
  seurat_object = combined_seu_filtered, 
  features = mono_genes_present, 
  reduction = "umap", 
  pt.size = 0.8,
  order = TRUE,                
  split.by = "timepoint",      
  num_columns = 2              
)


print(mono_split_check)

ggsave("mono_genes_split_by_timepoint.png", 
       plot = mono_split_check, 
       height = length(mono_genes_present) * 6, # 每行基因给 4 单位高度
       width = 14,                             # 宽度给 15 确保横向三个点能排开
       limitsize = FALSE,
       dpi = 300)

# ==== So many erythocytes...====
ery_genes <- c(
  "HBA1", "HBA2", "HBB",  # 血红蛋白系列 (Hemoglobins)
  "SLC4A1",                              # 带3蛋白 (Anion exchanger 1)
  "CA1"                                  # 碳酸酐酶 (Carbonic anhydrase I)
)

ery_genes_present <- intersect(ery_genes, rownames(combined_seu_filtered))

# 绘图检查分布
ery_check <- FeaturePlot_scCustom(
  seurat_object = combined_seu_filtered, 
  features = ery_genes_present, 
  reduction = "umap", 
  pt.size = 0.8,
  order = TRUE  
)

ggsave("erythocytes_genes_on_10clusters_myeloid.png", ery_check, height = 21, width = 15, dpi = 300)


FeaturePlot_scCustom(
  seurat_object = combined_seu_filtered, 
  features = "HSPA1A", 
  reduction = "umap", 
  pt.size = 0.8,
  order = TRUE  
)

# ===== Check things by timepoint =====
times <- unique(combined_seu_filtered$timepoint)

for (t in times) {
  sub_obj <- subset(combined_seu_filtered, subset = timepoint == t)
  
  # 2. 核心：检测哪些基因在当前子集里有表达
  # 计算每个基因在当前子集中的最大表达量，只有 > 0 的才保留
  exp_data <- GetAssayData(sub_obj, layer = "data") # 确保已经跑过 NormalizeData
  valid_genes <- mono_genes_present[rowSums(exp_data[mono_genes_present, , drop = FALSE]) > 0]
  
  if (length(valid_genes) == 0) {
    message("Timepoint ", t, " has no expression for the selected genes. Skipping...")
    next
  }
  
  # 3. 绘制密度图 (只用有表达的基因，防止 KDE 报错)
  p <- Plot_Density_Custom(
    seurat_object = sub_obj,
    features = valid_genes,
    reduction = "umap",
    pt.size = 1
  )
  
  # 4. 美化与保存
  p_final <- p + plot_annotation(
    title = paste("Density Plot - Timepoint:", t),
    subtitle = paste("Showing", length(valid_genes), "active markers"),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
  )
  
  file_name <- paste0("density_plot_monocytes_", t, ".png")
  ggsave(file_name, p_final, width = 12, height = 10, dpi = 300)
  
  message("Successfully saved: ", file_name, " (using genes: ", paste(valid_genes, collapse = ", "), ")")
}


# ===== Check things by stimulus =====
stims <- unique(combined_seu_filtered$stimulus)

for (s in stims) {
  sub_obj <- subset(combined_seu_filtered, subset = stimulus == s)
  
  exp_data <- GetAssayData(sub_obj, layer = "data") # 确保已经跑过 NormalizeData
  valid_genes <- mono_genes_present[rowSums(exp_data[mono_genes_present, , drop = FALSE]) > 0]
  
  if (length(valid_genes) == 0) {
    message("Stims ", s, " has no expression for the selected genes. Skipping...")
    next
  }
  
  # 3. 绘制密度图 (只用有表达的基因，防止 KDE 报错)
  p <- Plot_Density_Custom(
    seurat_object = sub_obj,
    features = valid_genes,
    reduction = "umap",
    pt.size = 2
  )
  
  # 4. 美化与保存
  p_final <- p + plot_annotation(
    title = paste("Density Plot - Stimulus:", s),
    subtitle = paste("Showing", length(valid_genes), "active markers"),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
  )
  
  print(p_final)
  
  file_name <- paste0("density_plot_monocytes_stimulus", s, ".png")
  ggsave(file_name, p_final, width = 12, height = 10, dpi = 300)
  
  message("Successfully saved: ", file_name, " (using genes: ", paste(valid_genes, collapse = ", "), ")")
}

# ==== violin plot to visualize changes in corresponding cluster ====
DefaultAssay(combined_seu_filtered) <- "SCT" # 或者 "RNA"

table(combined_seu_filtered$stimulus)

for (gene in mono_genes_present) {
  
  # 使用 VlnPlot_scCustom 绘制
  p_vln <- VlnPlot_scCustom(
    seurat_object = combined_seu_filtered,
    features = gene,
    group.by = "predicted.celltype.l3",   # 横轴：你的 Cluster
    split.by = "stimulus",          # 核心：在 Cluster 内按刺激组切分
    pt.size = 0,                 # 散点大小，900个细胞建议小一点
    colors_use = c("#F075AE", "#F7DB91", "#FD8A6B", "#9BC264", "#79C9C5") # 经典蓝红配色（可根据 stimulus 数量调整）
  ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Gene:", gene, "- Split by Stimulus"))
  
  # 打印预览
  print(p_vln)
  
  # 自动保存
  file_name_vln <- paste0("bowen_myeloid_10cluster_", gene, "_.png")
  ggsave(file_name_vln, p_vln, width = 10, height = 6, dpi = 300)
  
  message("Saved Violin Plot for: ", gene)
}
