# ====== 🫐 INFO ======
#
# This script is (attemping) analyzing the scRNA-seq data of trained monocytes.
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

# =====🍇 Library ====
library(Seurat)
library(dplyr)
library(harmony)
library(readr)

# Notes
#In a 10x experiment, you get one matrix for 5,000 cells. Here, you have 43 matrices for ~384 cells each. 
#To treat this like 10x data, you must:
#Iteratively Load each plate.
#Attach Metadata (Stimulus, Donor, Timepoint) to each cell immediately.
#Merge everything into one "Big Seurat Object."
# Standardize Gene Names (remove those __chr suffixes so they match common gene lists).


setwd("C:/Users/denxiack/Documents/Bowen_zhang/Bowen_Zhang")
files <- list.files(pattern = "\\.ReadCounts\\.tsv$")
files

csv_test <- read_tsv("RMC-SM-003_HM3JFBGX9_S5_R2.ReadCounts.tsv")
tail(csv_test$GENEID)

metadata_bowen <- read_csv("Metadata_TrainedImmunity.csv")
summary(metadata_bowen)
View(metadata_bowen)

  process_plate <- function(i) {
    row <- metadata_bowen[i, ]
    plate_id <- row$Plate_ID
    
    #Match File with Plate ID...
    files <- list.files(pattern = "\\.ReadCounts\\.tsv$")
    file_path <- files[grep(plate_id, files)]
    
    if(length(file_path) == 0) return(NULL)
    
    counts <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1)
    rownames(counts) <- gsub("__chr.*", "", rownames(counts))
    
    # Prefix cell names to remove duplicates
    colnames(counts) <- paste0(plate_id, "_", colnames(counts))
    
    # Create Seurat Object + add meta data
    obj <- CreateSeuratObject(counts = counts, project = "TrainedImmunity")
    obj$stimulus <- row$`Training stimulus`
    obj$donor    <- row$Donor
    obj$timepoint <- row$Timepoint
    
    return(obj)
  }

# Loop
obj_list <- lapply(1:nrow(metadata_bowen), process_plate)
obj_list <- Filter(Negate(is.null), obj_list) # Remove empty entries

# Merging them into seurat obj? - here no longer a marker for whether PBMC or monocytes.
combined_seu <- merge(obj_list[[1]], y = obj_list[2:length(obj_list)], add.cell.ids = NULL)
saveRDS(combined_seu, "results/bowen_merged_raw.rds")

# ====🪀 Manuall filtering =====
combined_seu_filtered <- subset(combined_seu, subset = nFeature_RNA > 500 & nCount_RNA > 1000)
ncol(combined_seu_filtered)

# ====💜Seurat procedure ====
combined_seu_filtered <- SCTransform(combined_seu_filtered, verbose = FALSE)
combined_seu_filtered <- RunPCA(combined_seu_filtered, npc = 20, verbose = FALSE)

ElbowPlot(combined_seu_filtered)

combined_seu_filtered <- RunUMAP(combined_seu_filtered, dims = 1:20)
combined_seu_filtered <- FindNeighbors(combined_seu_filtered, dims = 1:20)
combined_seu_filtered <- FindClusters(combined_seu_filtered, resolution = 0.4)

p1 <- DimPlot(combined_seu_filtered, reduction = "umap", group.by = "stimulus")
p2 <- DimPlot(combined_seu_filtered, reduction = "umap", group.by = "donor")
p3 <- DimPlot(combined_seu_filtered, reduction = "umap", group.by = "seurat_clusters")

p_sct_sum <- p1 + p2 + p3
p_sct_sum  #!!!! really does not look integrated! 
#ggsave("sct_transformed_clustering.png", p_sct_sum, height = 7, width = 21, dpi = 300)
getwd()

# =====🪻 LogNorm, no Seurat ====
combined_seu_filtered <- NormalizeData(combined_seu_filtered)
combined_seu_filtered <- FindVariableFeatures(combined_seu_filtered, selection.method = "vst", nfeatures = 2000)
combined_seu_filtered <- ScaleData(combined_seu_filtered)

combined_seu_filtered <- RunPCA(combined_seu_filtered, npcs = 30, verbose = FALSE)

# This forces Donor A and Donor B to occupy the same space !!!
combined_seu_filtered <- RunHarmony(combined_seu_filtered, group.by.vars = "donor")

combined_seu_filtered <- RunUMAP(combined_seu_filtered, reduction = "harmony", dims = 1:15) 
combined_seu_filtered <- FindNeighbors(combined_seu_filtered, reduction = "harmony", dims = 1:15)
combined_seu_filtered <- FindClusters(combined_seu_filtered, resolution = 0.5) 

p1 <- DimPlot(combined_seu_filtered, reduction = "umap", group.by = "stimulus")
p2 <- DimPlot(combined_seu_filtered, reduction = "umap", group.by = "donor")
p3 <- DimPlot(combined_seu_filtered, reduction = "umap", group.by = "seurat_clusters")

p1 + p2 + p3
p_log_sum <- p1 + p2 + p3
ggsave("log_transformed_clustering.png", p_log_sum, height = 7, width = 21, dpi = 300)

getwd()

# ==== 🧿 No original anno.. ====
#metadata_bowen$`Cell type`  #only has monocytes or PBMCs, does not help. Where is their annotated clusters?

# --- check publication gene list ---
#gene_list <- c("IL1B", "PTGS2", "CXCL10","CXCL11","CD14","CD83", "CD68","C1QA","C1QB","C1QC","MERTK",
#               "SLC16A1")
#FeaturePlot_scCustom(seurat_object = combined_seu_filtered, 
#                     features = gene_list, 
#                     num_columns = 3, 
#                     pt.size = 0.5,
#                     order = FALSE)
# -- check clusters --
#combined_seu_filtered$seurat_clusters
# ==== 🧤 Find top markers ====
#combined_seu_filtered <- PrepSCTFindMarkers(combined_seu_filtered)
#all_markers <- FindAllMarkers(combined_seu_filtered, 
#                              only.pos = TRUE, 
#                              min.pct = 0.25, 
#                              logfc.threshold = 0.25,
#                              verbose = FALSE)

#write.csv(all_markers, "bowen_combined_seu_all_markers.csv", row.names = FALSE)

#top10_markers <- all_markers %>%
#  group_by(cluster) %>%
#  slice_max(n = 10, order_by = avg_log2FC)

#top20_markers <- all_markers %>%
#  group_by(cluster) %>%
#  slice_max(n = 20, order_by = avg_log2FC)


#top10_genes <- unique(top10_markers$gene)

#seurat_integrated <- ScaleData(combined_seu_filtered, 
#                               features = top10_genes, 
#                               verbose = FALSE)

#p_dotplot <- DotPlot(
#  object = combined_seu_filtered, 
#  features = top10_genes, 
#  group.by = "seurat_clusters", 
#  cols = c("#0F2854", "#F375C2"), 
#  dot.scale = 6                
#) + 
#  RotatedAxis() +              
#  theme(
#    axis.text.x = element_text(size = 8, face = "italic"),
#    legend.position = "right"
#  ) +
#  labs(
#    title = "Top 10 Marker Genes Dot Plot per Cluster",
#    subtitle = "Size: % Expressed | Color: Avg Expression"
#  )

#print(p_dotplot)
#ggsave("bowen_myeloid_subcluster_topmarkers.png",p_dotplot, height = 8, width = 14, dpi = 300) 
# Weird lookking, do it the old fashion again:

#heatmap_plot <- DoHeatmap(combined_seu_filtered, features = top10_markers$gene) +
#  ggtitle("Top 10 Marker Genes per Cluster")
#print(heatmap_plot) # Probably failed
#View(top10_markers)
#ggsave("bowen_myeloid_subcluster_topmarkers_heatmaps.png",p_dotplot, height = 8, width = 14, dpi = 300) 

#saveRDS(combined_seu_filtered, "results/bowen_integrated_processed.rds")

# ==== 🎮 Azimuth ====
#combined_seu_filtered <- RunAzimuth(
#  query = combined_seu_filtered, 
#  reference = "pbmcref")#

#head(combined_seu_filtered@meta.data) #Works

#dp <- DimPlot(combined_seu_filtered, 
#              reduction = "umap",  
#              group.by = "predicted.celltype.l3", 
#              label = TRUE, 
#              repel = TRUE, 
#              pt.size = 1.5,
#              raster=TRUE) + 
#  ggtitle("Annotated Azimuth Bowen")

#dp <- DimPlot_scCustom(
#  seurat_object = combined_seu_filtered, 
#  reduction = "umap", 
#  group.by = "predicted.celltype.l3", 
#  label = TRUE,          
#  repel = TRUE,
#  pt.size = 1.5,        
#  label.size = 3,
#  colors_use = DiscretePalette_scCustomize(num_colors = 13,
#                                           palette = "varibow")
#) + 
#  ggtitle("Annotated Azimuth Bowen (L3)") +
#  theme(plot.title = element_text(hjust = 0.5))

#print(dp)
#saveRDS(combined_seu_filtered, "bowen_final_azimuth_nnotated.rds")

# ==== 🟣 LLM Top Markers ====
#new_cluster_ids <- c(
#  "0"  = "TREM2+ Lipid-Mac",
#  "1"  = "IFN-gamma-Resp Mac",
#  "2"  = "Inflammatory-Chemokine Mac",
#  "3"  = "Technical Artifact (ERCC)",
#  "4"  = "Non-classical Monocyte",
#  "5"  = "ISG-high Mac",
#  "6"  = "Glycolytic/Hypoxic Mac",
#  "7"  = "Epithelial-like Contam",
#  "8"  = "Rare Myeloid",
#  "9"  = "moDC",
# "10" = "Stress-Response Mac",
#  "11" = "VSTM4+ Myeloid"
#)
#combined_seu_filtered$llm_annotation <- unname(new_cluster_ids[as.character(combined_seu_filtered$seurat_clusters)])
#head(combined_seu_filtered$llm_annotation)

#p_llm <- DimPlot_scCustom(
#  seurat_object = combined_seu_filtered, 
#  reduction = "umap", 
#  group.by = "llm_annotation", 
#  label = TRUE,          
#  repel = TRUE,
#  pt.size = 1.5,        
#  label.size = 3,
#  colors_use = DiscretePalette_scCustomize(num_colors = 12,
#                                           palette = "alphabet")
#) + 
#  ggtitle("Annotated mLLm Bowen") +
#  theme(plot.title = element_text(hjust = 0.5))

#p_llm | dp
#combi <- p_llm + dp

#ggsave("combined_mLLm_vs_azimuth_depth.png", combi, height = 7, width = 18, dpi = 300)

# TWO ANNOTATIONS SAVED HERE
#saveRDS(combined_seu_filtered, "bowen_final_annotated_v2.rds")


# ==== ++++ Refine with Publication ++++ ====