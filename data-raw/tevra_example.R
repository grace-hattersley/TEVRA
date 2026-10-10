
library(SingleCellExperiment)

sce_obj <- readRDS('/Users/gracehattersley/Library/Mobile Documents/com~apple~CloudDocs/Nhood_VDJ/R_analysis/Dataset_1/RDS/SCE.rds')

colnames(sce_obj@colData)
sce_obj@colData <- sce_obj@colData[,!(colnames(sce_obj@colData) %in% c('RNA_snn_res.0.5', 'seurat_clusters', 'old.ident', 'RNA_snn_res.0.8'))]
colnames(sce_obj@colData)

names(sce_obj@metadata)
sce_obj@metadata[c('nhood_umap', 'nhood_diversity', 'annotations','merged_nhood_analysis',
                   'nhood_cdr3_distances_levenshtein', 'cdr3_merged_levensthein',
                   'nhood_cdr3_distances_blosum', 'cdr3_merged_blosum')] <- NULL
names(sce_obj@metadata)

assayNames(sce_obj)
assay(sce_obj, "scaledata") <- NULL
assay(sce_obj, "counts") <- NULL

dec  <- scran::modelGeneVar(sce_obj)
hvgs <- scran::getTopHVGs(dec, n = 2000)
sce_obj <- sce_obj[rownames(sce_obj) %in% hvgs, ]
dim(sce_obj)

format(object.size(sce_obj), units = "MB")

tevra_example <- sce_obj

usethis::use_data(tevra_example, compress = "xz", overwrite = TRUE)

file.size("data/tevra_example.rda") / 1e6


