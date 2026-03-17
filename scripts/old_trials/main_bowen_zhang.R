# ========================== 🫐 INFO =============================
# Project: Trained Immunity in Monocytes (Bowen Zhang Dataset)
# Author: Xiaoyue Deng
# Date: 17-03-2026
# 
# Script Objective:
#   1. Process Plate-based scRNA-seq data from E-MTAB-9702.
#   2. Isolate purified Monocyte populations from heterogeneous PBMC plates.
#   3. Characterize Trained Immunity signatures (HIF-1/Glycolysis, Mac1, Mac2).
#   4. Validate findings in cross-disease myeloid contexts (AD / Psoriasis).
#
# Key functions of this script:
#   
#   - Custom Plate Loader: Iteratively process .tsv plates with metadata-based filtering.
#   
#   - Dual QC: Used stricter filtering standards and added plate-specific controls 
#     (ERCC Spike-in removal).
#  
#   - State-specific pathway analysis (Classical, Non-classical, HIF-Glycolysis, 
#     Inflammatory/Repair Macrophage).
# 
#   - SCTransform-based integration split by stimulus.
#
# Goal Identifying if Trained Immunity signatures (e.g., IL1B, HIF1A) 
# are more pronounced in a specific sub-section of myeloid populations in AD.
#
# References:
#   - Publication: https://www.jci.org/articles/view/147719
#   - Data: https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-9702
#   - Raw Logic: https://github.com/CiiM-Bioinformatics-group/trained_immunity
# ============ 🛟 Questions to be solved ===============
#
#   Check dataset from other AD samples
#
# =====🍇 Library ====
library(Seurat)
library(dplyr)
library(harmony)
library(readr)
library(ggplot2)
library(AUCell)
library(scCustomize)

# === Process data ===
set.seed(666)
setwd("C:/Users/denxiack/Documents/Bowen_zhang/Bowen_Zhang")
files <- list.files(pattern = "\\.ReadCounts\\.tsv$")


csv_test <- read_tsv("RMC-SM-003_HM3JFBGX9_S5_R2.ReadCounts.tsv")
tail(csv_test$GENEID)

metadata_bowen <- read_csv("Metadata_TrainedImmunity.csv")
summary(metadata_bowen)
metadata_bowen$`Cell type`  # Get all monocytes


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
# Peepping 
table(combined_seu_mono$origin_cell_type)
head(rownames(combined_seu_mono))
head(combined_seu_mono@meta.data)
grep("^ERCC-", rownames(combined_seu_mono), value = TRUE)
table(combined_seu_mono$stimulus)
head(combined_seu_mono@meta.data) 

# ====== 👑 calculate QC thresholds ======
combined_seu_mono[["percent.mt"]] <- PercentageFeatureSet(combined_seu_mono, pattern = "^MT-")
combined_seu_mono[["percent.ercc"]] <- PercentageFeatureSet(combined_seu_mono, pattern = "^ERCC-")

# ====== 👑 filtering ======
# remove underpresented genes 
combined_seu_mono <- JoinLayers(combined_seu_mono)
counts <- GetAssayData(combined_seu_mono, layer = "counts")

keep_genes <- rowSums(counts > 0) >= 3
combined_seu_mono <- subset(combined_seu_mono, features = names(keep_genes[keep_genes]))

# filter the way paper did
combined_seu_filtered <- subset(
  combined_seu_mono, 
  subset = nFeature_RNA >= 100 & 
    nFeature_RNA <= 7000 & 
    percent.mt <= 25 & 
    percent.ercc <= 90
)

# Rm erythocytes such as ERCC
non_ercc_genes <- grep("^ERCC-", rownames(combined_seu_filtered), value = TRUE, invert = TRUE)
combined_seu_filtered <- subset(combined_seu_filtered, features = non_ercc_genes)

# ====== 👑 SCTransform by condition ======
obj.list <- SplitObject(combined_seu_filtered, split.by = "stimulus")

obj.list <- lapply(X = obj.list, FUN = SCTransform, vars.to.regress = "percent.mt", verbose = FALSE)

features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 2000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", anchor.features = features)
combined_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

Assays(combined_integrated)

# ====== 👑 deconclustering ======
DefaultAssay(combined_integrated) <- "integrated"

combined_integrated <- RunPCA(combined_integrated, npcs = 30, verbose = FALSE)
ElbowPlot(combined_integrated, ndims = 30)

combined_integrated <- RunUMAP(combined_integrated, reduction = "pca", dims = 1:20)
combined_integrated <- FindNeighbors(combined_integrated, reduction = "pca", dims = 1:20)
combined_integrated <- FindClusters(combined_integrated, resolution = 0.5)

p1 <- DimPlot(combined_integrated, reduction = "umap", group.by = "timepoint") + ggtitle("T1 vs T2")
p2 <- DimPlot(combined_integrated, reduction = "umap", group.by = "stimulus") + ggtitle("Stimuli")
p3 <- DimPlot(combined_integrated, reduction = "umap", group.by = "seurat_clusters") + ggtitle("Clusters")
p4 <- DimPlot(combined_integrated, reduction = "umap", group.by = "donor") + ggtitle("Donor")
umap_combined <- p1/p4 | p2 / p3
umap_combined
ggsave("umap_bowen_unannotated_merge.png", umap_combined, height = 12, width = 15, dpi = 300)

# ===== 👑 Annotation ==== 
# Switch back to SCT assay
DefaultAssay(combined_integrated) <- "SCT"

# --- T1 Monocytes ---
# Classical: CD14, S100A9
# Non-classical: FCGR3A (CD16), VMO1
# Intermediate: CD14, FCGR3A, HLA-DRA
t1_markers <- c("CD14", "S100A9", "FCGR3A", "VMO1", "HLA-DRA")

# --- T2 Macrophages
# Mac-1: CD36, CD83, IL1B
# Mac-2: IL1RN, MSR1, SPP1
t2_markers <- c("CD36", "CD83", "IL1B", "IL1RN", "MSR1", "SPP1")

t1_feature <- FeaturePlot(combined_integrated, features = t1_markers, ncol = 3)
t2_feature <- FeaturePlot(combined_integrated, features = t2_markers, ncol = 3)

ggsave("t1_marekers_on_feature_plot.png", t1_feature, height = 12, width = 18, dpi = 300)
ggsave("t2_marekers_on_feature_plot.png", t2_feature, height = 12, width = 18, dpi = 300)

# ==== 👑 Find All Markers ====
combined_integrated <- PrepSCTFindMarkers(combined_integrated)
all_markers <- FindAllMarkers(combined_integrated, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)


top5 <- all_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
print(top5)
View(top5)

# ====== 👑 Top 10 Markers ======
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DefaultAssay(combined_integrated) <- "SCT"
combined_sub <- subset(combined_integrated, downsample = 100)
combined_sub <- ScaleData(combined_sub, features = top10_markers$gene)

bowen_heatmap <- DoHeatmap(combined_sub, 
          features = top10_markers$gene, 
          size = 3,           
          angle = 0,          
          draw.lines = TRUE) + 
  theme(axis.text.y = element_text(size = 6)) 

ggsave("heatmap.png", bowen_heatmap, height = 8, width = 12, dpi = 300)


# ===== 👑 Annotattion based on Featureplot gene expression =====
DefaultAssay(combined_integrated) <- "SCT"

p_mono <- FeaturePlot_scCustom(combined_integrated, 
                               features = c("CD14", "S100A9", "VCAN", "S100A8", "FCGR3A",
                                            "VMO1", "LST1", "MS4A7"),
                               num_columns = 3, order = FALSE)
p_mono  # Classical : high S100A9, S100A8, VCAN | low GCFR3A
        # Non-classical : VMO1, LST1, MS4A7 | low CD14
ggsave("p_mono_featureplot.png", p_mono, height = 9, width = 15, dpi = 300)

p_hif <- FeaturePlot_scCustom(combined_integrated,  
                              features = c("HIF1A", "LDHA", "ENO1", "GAPDH"), 
                              colors_use = viridis_plasma_dark_high, 
                              num_columns = 2, order = FALSE)
p_hif
ggsave("p_hif_featureplot.png", p_mono, height = 9, width = 12, dpi = 300)

# C. Pro-inflammatory Mac1
p_mac1 <- FeaturePlot_scCustom(combined_integrated, 
                               features = c("IL1B", "CD83", "CCL3", "IL23A"), 
                               colors_use = viridis_plasma_dark_high, 
                               num_columns = 2, order = FALSE)

p_mac1
ggsave("p_mac1_featureplot.png", p_mac1, height = 9, width = 12, dpi = 300)

# D. T2 Mac2
p_mac2 <- FeaturePlot_scCustom(combined_integrated, 
                               features = c("IL1RN", "MSR1", "SPP1", "LPL"), 
                               colors_use = viridis_plasma_dark_high, 
                               num_columns = 2, order = FALSE)
p_mac2
ggsave("p_mac2_featureplot.png", p_mac2, height = 9, width = 12, dpi = 300)

# UGDH-AS1 cells
p_lnc <- FeaturePlot_scCustom(combined_integrated, 
                              features = c("UGDH-AS1", "KCNQ1OT1"), 
                              colors_use = viridis_plasma_dark_high, 
                              num_columns = 2, order = FALSE)
p_lnc
ggsave("p_lnc_featureplot.png", p_lnc, height = 6, width = 12, dpi = 300)


mo_dc <- FeaturePlot_scCustom(combined_integrated, 
                     features = c("CD86", "HLA-DRA", "CD1C", "FCER1A"), 
                     colors_use = viridis_plasma_dark_high,
                     num_columns = 2, order = FALSE)
mo_dc
ggsave("p_mo_dc_featureplot.png", mo_dc, height = 9, width = 12, dpi = 300)


# ==== 👑 save the feature plots ====

saveRDS(combined_integrated, file = "combined_seu_mono_manual_annotation.rds")
# combined_integrated <- readRDS("combined_seu_mono_final_annotated.rds")

# ===== 👑 Manual Annotation - WIP =====
new_names <- c(
  "0" = "intermediate monocytes",
  "1" = "Macrophages_1",
  "2" = "Macrophages_1/2 mix",
  "3" = "classical monocytes",
  "4" = "unpolarized",
  "5" = "macrophages_2",
  "6" = "Macrophages_1",
  "7" = "macrophages_2",
  "8" = "UGDH-AS1+ cells",
  "9" = "HIF1 signaling cells",
  "10" = "non-classical monocytes",
  "11" = "monocytic DCs"
  
)

Idents(combined_integrated) <- "seurat_clusters"
combined_integrated <- RenameIdents(combined_integrated, new_names)
combined_integrated$manual_annotation <- Idents(combined_integrated)

head(combined_integrated@meta.data)
table(combined_integrated$manual_annotation) 

p5 <- DimPlot(combined_integrated,
                                  reduction = "umap",
                                  group.by = "manual_annotation",
                                  label = TRUE,
                                  repel = TRUE,
                                  pt.size = 0.5) +
  scale_color_manual(values = c("#F0E491","#BBC863","#F08787","#FFDBB6","#C2E2FA",
                                "#F5D3C4","#9B5DE0","#FCB53B","#3396D3","#67C090"))+
  ggtitle("Annotation Manual")
p5
ggsave("annotated UMAP.png", p5, height = 9, width = 12, dpi = 300)


# ==== focus: unpolarized population ====

# If T1 has high classical mono score but low HIF_score and moDC -> resting monocytes
# If T2 has low Mac1 and Mac2 -> unpolarized macrophage

# ===== DEBUG: Using modulescore ====
module_list <- list(
  Classical_Score = c("CD14", "S100A9", "VCAN", "S100A8"),
  NonClassical_Score = c("FCGR3A", "VMO1", "LST1", "MS4A7"),
  HIF_Glycolysis_Score = c("HIF1A", "LDHA", "ENO1", "GAPDH"),
  Mac1_Inflam_Score = c("IL1B", "CD83", "CCL3", "IL23A"),
  Mac2_Repair_Score = c("IL1RN", "MSR1", "SPP1", "LPL"),
  moDC_Score = c("CD86", "HLA-DRA", "CD1C", "FCER1A")
)

combined_integrated <- AddModuleScore(
  object = combined_integrated,
  features = module_list,
  name = "Module_Score_"
)

score_names <- names(module_list)
for(i in 1:length(score_names)){
  colnames(combined_integrated@meta.data)[grep(paste0("Module_Score_", i, "$"), colnames(combined_integrated@meta.data))] <- score_names[i]
}


p_all_scores <- FeaturePlot_scCustom(
  combined_integrated, 
  features = score_names, 
  colors_use = viridis_plasma_dark_high,
  num_columns = 3, 
  order = TRUE)     

p_all_scores
ggsave("module_score_analysis.png", p_all_scores, height = 6, width = 12, dpi = 300 )
# ==== DEBUG : Using AUCell  -  WIP ====


# ==== 🥭 PBMC whole: WIP =====
# Reasons for halt: No need for PBMC genes. 
# Will probably see a lot of erythoblasts / mRNA containmination
set.seed(666)
setwd("C:/Users/denxiack/Documents/Bowen_zhang/Bowen_Zhang")
files <- list.files(pattern = "\\.ReadCounts\\.tsv$")


csv_test <- read_tsv("RMC-SM-003_HM3JFBGX9_S5_R2.ReadCounts.tsv")
tail(csv_test$GENEID)

metadata_bowen <- read_csv("Metadata_TrainedImmunity.csv")
summary(metadata_bowen)
metadata_bowen$`Cell type`  # Get all yourmonocytes

# NOT selecting only monocytes, but whole 


# ==== 🥭 To Do List ====
# Original paper finding : HIF-1 pathway particularly triggered by Beta-glucan.
#
# To Do 1:
# Check HIF-1 Glycolysis score in AD datasets 
#
# To Do 2:
# Check Macrophage 1 / Macrophage 2 ratio in AD
#
#
# To Do 3: 
# Cell to cell communication check 

# ==== 🥭 Dataset Index ====
# Choice of Datasets 
# Zhang_blood 


# ===== 😸 Zhang ====== 
getwd()
zhang_pbmc <- readRDS("transformed_blood_renamed.rds")
zhang_pbmc$named_clusters
table(zhang_pbmc$named_clusters)

zhang_p <- DimPlot(zhang_pbmc, 
                   reduction = "umap",  
                   group.by = "named_clusters", 
                   label = TRUE, 
                   repel = TRUE, 
                   pt.size = 0.5,
                   raster=TRUE) + 
  ggtitle("Zhang Blood")
print(zhang_p)

myeloid_zhang <- readRDS("myeloid_zhang_Annotated.rds")
myeloid_zhang$llm_prediction
myeloid_z <- DimPlot(myeloid_zhang, 
                   reduction = "umap",  
                   group.by = "llm_prediction", 
                   label = TRUE, 
                   repel = TRUE, 
                   pt.size = 0.5,
                   raster=TRUE) + 
  ggtitle("Myeloid Zhang")
print(myeloid_z)

# ===== 😸  Marker in AD vs HC ======
summary(myeloid_zhang$disease) # BUT no HC.... Liu and Rojahn are both skin.
# Research question: Are there trained immunity signatures that differentiates AD and Psoriasis

pz_mono <- FeaturePlot_scCustom(myeloid_zhang, 
                               features = c("CD14", "S100A9", "VCAN", "S100A8", "FCGR3A",
                                            "VMO1", "LST1", "MS4A7"),
                               num_columns = 3, order = FALSE)
pz_mono  # Classical : high S100A9, S100A8, VCAN | low GCFR3A
# Non-classical : VMO1, LST1, MS4A7 | low CD14
ggsave("zhang/p_zhang_mono_featureplot.png", p_mono, height = 9, width = 15, dpi = 300)

pz_hif <- FeaturePlot_scCustom(myeloid_zhang,  
                              features = c("HIF1A", "LDHA", "ENO1", "GAPDH"), 
                              colors_use = viridis_plasma_dark_high, 
                              num_columns = 2, order = FALSE)
pz_hif # high in all....
ggsave("zhang/p_zhang_hif_featureplot.png", p_mono, height = 9, width = 12, dpi = 300)

# C. Pro-inflammatory macrophage
pz_mac1 <- FeaturePlot_scCustom(myeloid_zhang, 
                               features = c("IL1B", "CD83", "CCL3", "IL23A"), 
                               colors_use = viridis_plasma_dark_high, 
                               num_columns = 2, order = FALSE)

pz_mac1
ggsave("zhang/p_zhang_mac1_featureplot.png", p_mac1, height = 9, width = 12, dpi = 300)

# D. T2 Mac2
pz_mac2 <- FeaturePlot_scCustom(myeloid_zhang, 
                               features = c("IL1RN", "MSR1", "SPP1", "LPL","IL1B"), # should have no IL1B 
                               colors_use = viridis_plasma_dark_high, 
                               num_columns = 2, order = TRUE)
pz_mac2
ggsave("zhang/p_mac2_featureplot.png", p_mac2, height = 9, width = 12, dpi = 300)

# UGDH-AS1 cells
pz_lnc <- FeaturePlot_scCustom(myeloid_zhang, 
                              features = c("UGDH-AS1", "KCNQ1OT1"), 
                              colors_use = viridis_plasma_dark_high, 
                              num_columns = 2, order = FALSE)
pz_lnc
ggsave("zhang/p_zhang_lnc_featureplot.png", p_lnc, height = 6, width = 12, dpi = 300)


moz_dc <- FeaturePlot_scCustom(myeloid_zhang, 
                              features = c("CD86", "HLA-DRA", "CD1C", "FCER1A"), 
                              colors_use = viridis_plasma_dark_high,
                              num_columns = 2, order = FALSE)
moz_dc
ggsave("zhang/p_zhang_mo_dc_featureplot.png", mo_dc, height = 9, width = 12, dpi = 300)


# ==== Myeloid Bowen Annotation ====
myeloid_zhang$seurat_clusters
DimPlot(myeloid_zhang, 
        reduction = "umap",  
        group.by = "seurat_clusters", 
        label = TRUE, 
        repel = TRUE, 
        pt.size = 1,
        raster=TRUE) 

# Double check what's in 8
Idents(myeloid_zhang) <- "seurat_clusters" 

small_cluster_markers <- FindMarkers(myeloid_zhang, 
                                     ident.1 = "8", 
                                     only.pos = TRUE, 
                                     min.pct = 0.25, 
                                     logfc.threshold = 0.25)

head(small_cluster_markers, n = 15)

manual_names <- c(
  "0" = "classical monocytes",
  "1" = "inflammatory monocytes",
  "2" = "classical monocytes",
  "3" = "CD16+ monocytes", # edge breaching into intermediate monocytes
  "5" = "moDCs",
  "8" = "Activated monocytes" # IL1B, CCL3...
)

Idents(myeloid_zhang) <- "seurat_clusters"
myeloid_zhang <- RenameIdents(myeloid_zhang, manual_names)
myeloid_zhang$manual_names <- Idents(myeloid_zhang)

head(myeloid_zhang@meta.data)
table(myeloid_zhang$manual_names) 

pz_mo_combined <- DimPlot(myeloid_zhang,
              reduction = "umap",
              group.by = "manual_names",
              label = TRUE,
              repel = TRUE,
              pt.size = 0.5) +
  scale_color_manual(values = c("#F0E491","#BBC863","#F08787","#C2E2FA",
                                "#9B5DE0","#FCB53B","#3396D3","#67C090"))+
  ggtitle("Annotation Manual")

pz_mo_combined
ggsave("zhang/manual_annotated UMAP.png", pz_mo_combined, height = 9, width = 12, dpi = 300)

# ===== Check Signature Genes on AD vs Pso ====
table(myeloid_zhang$disease)
features_to_plot <- c("IL1B", "HIF1A","PDE4B", "TNFAIP3", "CCL3")
sv <- Stacked_VlnPlot(seurat_object = myeloid_zhang, 
                features = features_to_plot, 
                group.by = "manual_names", 
                split.by = "disease", 
                plot_legend = TRUE,
                colors_use = c("#537D96", "#EB4C4C"))
sv
ggsave("zhang/stacked_violin.png", sv, height = 8, width = 12, dpi = 300)

# ==== IL1B featureplot ====
p <- FeaturePlot_scCustom(myeloid_zhang, 
                 features = "IL1B", 
                 split.by = "disease", 
                 order = TRUE, 
                 pt.size = 0.5) & 
  theme(legend.position = "right")

p
ggsave("zhang/IL1B.png", p, height = 6, width = 12, dpi = 300)


# ==== HIF ====
"HIF1A"
p <- FeaturePlot_scCustom(myeloid_zhang, 
                          features = "HIF1A", 
                          split.by = "disease", 
                          order = TRUE, 
                          pt.size = 0.5) & 
  theme(legend.position = "right")
p

ggsave("zhang/HIF1A.png", p, height = 6, width = 12, dpi = 300)
