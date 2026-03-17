library(Seurat)
library(dplyr)
library(harmony)
library(readr)
setwd("C:/Users/denxiack/Documents/Bowen_zhang/Bowen_Zhang")

# 1. 数据加载与元数据注入 (核心：处理 43 个 plates)
metadata_bowen <- read_csv("Metadata_TrainedImmunity.csv")
head(metadata_bowen)
files <- list.files(pattern = "\\.ReadCounts\\.tsv$")

process_plate_refined <- function(i) {
  row_data <- metadata_bowen[i, ]
  plate_id <- as.character(row_data$Plate_ID)
  
  file_path <- files[grep(plate_id, files)]
  if(length(file_path) == 0) return(NULL)
  
  counts <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1)
  rownames(counts) <- gsub("__chr.*", "", rownames(counts))
  colnames(counts) <- paste0(plate_id, "_", colnames(counts))
  
  obj <- CreateSeuratObject(counts = counts, project = "TrainedImmunity")
  for(col in colnames(metadata_bowen)) {
    obj[[col]] <- as.character(row_data[[col]])
  }
  
  return(obj)
}

# Merge into a Seurat.
obj_list <- lapply(1:nrow(metadata_bowen), process_plate_refined)
obj_list <- Filter(Negate(is.null), obj_list)
combined_raw <- merge(obj_list[[1]], y = obj_list[2:length(obj_list)])

print(colnames(combined_raw@meta.data))
print(unique(combined_raw$cell_group))
combined_raw@meta.data

# === Add back nCount and nFeature...? ===

combined_raw <- JoinLayers(combined_raw) # PLATE BASED WILL BE TREATED AS DIFFERENT LAYERS FOR SEURAT

combined_raw$nCount_RNA <- colSums(GetAssayData(combined_raw, assay = "RNA", layer = "counts"))
combined_raw$nFeature_RNA <- colSums(GetAssayData(combined_raw, assay = "RNA", layer = "counts") > 0)

combined_raw[["percent.mt"]] <- PercentageFeatureSet(combined_raw, pattern = "^MT-")

# Sanity Check 
head(combined_raw@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt")])


combined_qc <- subset(combined_raw, 
                      subset = nFeature_RNA > 100 & 
                        nFeature_RNA < 7000 & 
                        percent.mt < 25)

# Check how many cells, and how many genes left


# Make layers for plates
combined_qc[["RNA"]] <- split(combined_qc[["RNA"]], f = combined_qc$Plate_ID)

# paper sct
combined_qc <- SCTransform(combined_qc, 
                           vars.to.regress = "percent.mt", 
                           vst.flavor = "v2", 
                           verbose = TRUE)
# Add sct back to meta
DefaultAssay(combined_qc) <- "SCT"

combined_qc$nCount_SCT <- colSums(GetAssayData(combined_qc, assay = "SCT", layer = "counts"))
combined_qc$nFeature_SCT <- colSums(GetAssayData(combined_qc, assay = "SCT", layer = "counts") > 0)

VlnPlot(combined_qc, features = c("nCount_SCT", "nFeature_SCT"), group.by = "Timepoint")


# ==== Helen way of QC ====
