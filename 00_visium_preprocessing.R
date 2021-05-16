### Meylan et al 2021
### pre-processing of visium data from Space Ranger V1.1.0
### max.meylan@gmail.com

library(Seurat)
library(dplyr)
library(MCPcounter)
source("~/script/maxime.utils.r")

res_folder <- "~/projects/visium/results/"
#mcp var names an signatures
cell_types <- make.names(rownames(MCPcounter.estimate(MCPcounterExampleData)))
signatures <- get_signature(get_list = T)
################## 
### code for systematic analysis 
slide_list <- c("a_15v",
                "a_4",
                "b_1",
                "b_17",
                "b_12",
                "c_10",
                "c_3",
                "c_26",
                "c_2",
                "a_1",
                "a_18",
                "b_6",)
spatial_list <- sapply(slide_list,function(slide){
  print(slide)
  raw_data_directory <- paste0("/data/visium_ccRCC/processed_data/",slide,"/outs/")
  spatial_object <- Seurat::Load10X_Spatial(raw_data_directory)
  
  # Collect all genes coded on the mitochondrial genome
  mt.genes <- grep(pattern = "^MT-", x = rownames(spatial_object), value = TRUE)
  spatial_object$percent.mito <- (Matrix::colSums(spatial_object@assays$Spatial@counts[mt.genes, ])/Matrix::colSums(spatial_object@assays$Spatial@counts))*100
  
  #remove mt genes
  genes_to_keep <- setdiff(names(which(Matrix::rowSums(spatial_object@assays$Spatial@counts )>5)),mt.genes)
  
  spatial_object_subset <- subset(spatial_object,features =genes_to_keep, subset = nFeature_Spatial > 300 & percent.mito < 30)
  cat("Spots removed: ", ncol(spatial_object) - ncol(spatial_object_subset), "\n")
  cat("Genes kept: ", length(genes_to_keep),"from",nrow(spatial_object), "\n") 
  
  spatial_object_subset <- SCTransform(spatial_object_subset, assay = "Spatial", verbose = T)

  #compute mcp scores
  mcp_scores <- MCPcounter.estimate(spatial_object_subset[["SCT"]]@data,featuresType = "HUGO_symbols")
  rownames(mcp_scores) <- make.names(rownames(mcp_scores))
  cell_types <- rownames(mcp_scores)
  for(cell_type in cell_types){
    spatial_object_subset <- AddMetaData(object= spatial_object_subset,metadata = mcp_scores[cell_type,],col.name = cell_type)
  }
  # Unsupervised analysis
  spatial_object_subset <- RunICA(spatial_object_subset, assay = "SCT", verbose = FALSE)
  spatial_object_subset <- FindNeighbors(spatial_object_subset, reduction = "ica")
  spatial_object_subset <- FindClusters(spatial_object_subset, verbose = FALSE,resolution = 0.4)
  spatial_object_subset <- RunUMAP(spatial_object_subset, reduction = "ica", dims = 1:30)
  return(spatial_object_subset)
})
save.image(paste0(res_folder,Sys.Date(),"_","seurat_processed.RData"))
