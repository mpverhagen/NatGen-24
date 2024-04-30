library(Seurat)
library(ggplot2)
library(dplyr)
library(Signac)
path <- ("/mnt/data/User/Mathijs/Becker_polyps_CRC_singlecell")

### scATACseq Becker et al.: GSE201336 ####
####Import Fragment files ####
epithelial_celltypes_atac <- read.delim(paste(path, "/scATAC_GEO_Normal/epithelial_celltypes_atac.tsv", sep = ""))
epithelial_celltypes_atac$Cell_short <- sapply(strsplit( epithelial_celltypes_atac$Cell, "#"), "[", 2 )

# Make fragment objects for all conditions
# 1/8
B001_A_301_cells <- epithelial_celltypes_atac %>% dplyr::filter(Sample  == "B001-A-301-D")
B001_A_301 <- CreateFragmentObject(path = "/mnt/data/User/Mathijs/Becker_polyps_CRC_singlecell/scATAC_GEO_Normal/B001-A-301/GSM6058850_B001-A-301-D_20200804_fragments.tsv.gz", 
                                   cells = sapply(strsplit(B001_A_301_cells$Cell, "#"), "[", 2 ),
                                   verbose = TRUE,
                                   validate.fragments = TRUE)
B001_A_301_cells$ID <- paste("B001_A_301", sapply(strsplit(B001_A_301_cells$Cell, "#"), "[", 2 ), sep = "_")
B001_A_301_cells$short <- sapply(strsplit(B001_A_301_cells$Cell, "#"), "[", 2 )

# 2/8
B001_A_302_cells <- epithelial_celltypes_atac %>% dplyr::filter(Sample  == "B001-A-302-D")
B001_A_302 <- CreateFragmentObject(path = "/mnt/data/User/Mathijs/Becker_polyps_CRC_singlecell/scATAC_GEO_Normal/B001-A-302/GSM6058851_B001-A-302-D_fragments.tsv.gz", 
                                   cells = sapply(strsplit(B001_A_302_cells$Cell, "#"), "[", 2 ),
                                   verbose = TRUE,
                                   validate.fragments = TRUE)
B001_A_302_cells$ID <- paste("B001_A_302", sapply(strsplit(B001_A_302_cells$Cell, "#"), "[", 2 ), sep = "_")

# 3/8
B001_A_401_cells <- epithelial_celltypes_atac %>% dplyr::filter(Sample  == "B001-A-401-D")
B001_A_401 <- CreateFragmentObject(path = "/mnt/data/User/Mathijs/Becker_polyps_CRC_singlecell/scATAC_GEO_Normal/B001-A-401/GSM6058852_B001-A-401-D_fragments.tsv.gz", 
                                   cells = sapply(strsplit(B001_A_401_cells$Cell, "#"), "[", 2 ),
                                   verbose = TRUE,
                                   validate.fragments = TRUE)
# 4/8
B001_A_406_cells <- epithelial_celltypes_atac %>% dplyr::filter(Sample  == "B001-A-406-D")
B001_A_406 <- CreateFragmentObject(path = "/mnt/data/User/Mathijs/Becker_polyps_CRC_singlecell/scATAC_GEO_Normal/B001-A-406/GSM6058853_B001-A-406-D_fragments.tsv.gz", 
                                   cells = sapply(strsplit(B001_A_406_cells$Cell, "#"), "[", 2 ),
                                   verbose = TRUE,
                                   validate.fragments = TRUE)

# 5/8
B001_A_501_cells <- epithelial_celltypes_atac %>% dplyr::filter(Sample  == "B001-A-501-D")
B001_A_501 <- CreateFragmentObject(path = "/mnt/data/User/Mathijs/Becker_polyps_CRC_singlecell/scATAC_GEO_Normal/B001-A-501/GSM6058854_B001-A-501-D_fragments.tsv.gz", 
                                   cells = sapply(strsplit(B001_A_501_cells$Cell, "#"), "[", 2 ),
                                   verbose = TRUE,
                                   validate.fragments = TRUE)

# 6/8: R1
B004_A_004_R1_cells <- epithelial_celltypes_atac %>% dplyr::filter(Sample  == "B004-A-004-D")
B004_A_004_R1 <- CreateFragmentObject(path = "/mnt/data/User/Mathijs/Becker_polyps_CRC_singlecell/scATAC_GEO_Normal/B004-A-004_R1/GSM6058855_B004-A-004-D_20200715_fragments.tsv.gz", 
                                      cells = sapply(strsplit(B004_A_004_R1_cells$Cell, "#"), "[", 2 ),
                                      verbose = TRUE,
                                      validate.fragments = TRUE)

# 6/8: R2
B004_A_004_R2_cells <- epithelial_celltypes_atac %>% dplyr::filter(Sample  == "B004-A-004-D-R2")
B004_A_004_R2 <- CreateFragmentObject(path = "/mnt/data/User/Mathijs/Becker_polyps_CRC_singlecell/scATAC_GEO_Normal/B004-A-004_R2/GSM6058856_B004-A-004-D_20200817_fragments.tsv.gz", 
                                      cells = sapply(strsplit(B004_A_004_R2_cells$Cell, "#"), "[", 2 ),
                                      verbose = TRUE,
                                      validate.fragments = TRUE)

# 7/8
B004_A_008_cells <- epithelial_celltypes_atac %>% dplyr::filter(Sample  == "B004-A-008-D")
B004_A_008 <- CreateFragmentObject(path = "/mnt/data/User/Mathijs/Becker_polyps_CRC_singlecell/scATAC_GEO_Normal/B004-A-008/GSM6058857_B004-A-008-D_20200817_fragments.tsv.gz", 
                                   cells = sapply(strsplit(B004_A_008_cells$Cell, "#"), "[", 2 ),
                                   verbose = TRUE,
                                   validate.fragments = TRUE)

# 8/8
B004_A_204_cells <- epithelial_celltypes_atac %>% dplyr::filter(Sample  == "B004-A-204-D")
B004_A_204 <- CreateFragmentObject(path = "/mnt/data/User/Mathijs/Becker_polyps_CRC_singlecell/scATAC_GEO_Normal/B004-A-204/GSM6058858_B004-A-204-D_20200702_fragments.tsv.gz", 
                                   cells = sapply(strsplit(B004_A_204_cells$Cell, "#"), "[", 2 ),
                                   verbose = TRUE,
                                   validate.fragments = TRUE)

# Peaks based on on Sample or in /Combined folder, add the different fragment files
fragpath = dir('/mnt/data/User/Mathijs/Becker_polyps_CRC_singlecell/scATAC_GEO_Normal/Combined/',"*gz$",full.names=TRUE)

#### Start from 1Mb bed file ####
# Convert Windows.bed to peaks all format
windows.bed <- as.data.frame(read.table("/mnt/data/User/Mathijs/Becker_polyps_CRC_singlecell/Signac/GRCh38_no_alt.1Mbwindows.min0.92umapk36.noTelomere.noCentromere.autosome.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
colnames(windows.bed) <- c("Chr", "Start", "End")
windows.bed <- GenomicRanges::makeGRangesFromDataFrame(windows.bed)

# Create peak count matrix: All samples
bin.mat <- FeatureMatrix(fragments = list(B001_A_301,B001_A_302, B001_A_401, B001_A_406, B001_A_501, B004_A_004_R1, B004_A_004_R2, B004_A_008, B004_A_204),features = windows.bed)
#### Create Chromatin Assay  ####
chrom_assay <- CreateChromatinAssay(
  counts = bin.mat,
  sep = c(":", "-"),
  fragments = list(B001_A_301,B001_A_302, B001_A_401, B001_A_406, B001_A_501, B004_A_004_R1, B004_A_004_R2, B004_A_008, B004_A_204),
  min.cells = 10,
  min.features = 200
)

# Create Signac Object
Colon.atac <- CreateSeuratObject(counts = chrom_assay, assay = "bins")

#### Dim reduction ####
# Normalization and linear dimensional reduction
Colon.atac <- RunTFIDF(Colon.atac)
Colon.atac <- FindTopFeatures(Colon.atac, min.cutoff = 'q0')
Colon.atac <- RunSVD(Colon.atac)
DepthCor(Colon.atac)

# Non-linear dimension reduction and clustering
Colon.atac <- RunUMAP(object = Colon.atac, reduction = 'lsi', dims = 2:30)
epithelial_celltypes_atac_take <- epithelial_celltypes_atac %>% dplyr::filter(Sample %in% c("B001-A-301-D", "B001-A-302-D", "B001-A-401-D", "B001-A-406-D", "B001-A-501-D", "B004-A-004-D", "B004-A-004-D-R2","B004-A-008-D", "B004-A-204-D" ))

# Add meta 
Colon.atac$Sample <- NA
Colon.atac$Sample <-  epithelial_celltypes_atac_take$Sample[match(colnames(Colon.atac), epithelial_celltypes_atac_take$Cell_short )]

# Add Cell type information
Colon.atac$Celltype <- "Other"
Colon.atac$Celltype <-  epithelial_celltypes_atac$CellType[match(colnames(Colon.atac), epithelial_celltypes_atac$Cell_short )]

# Merge two replicate samples
Colon.atac$Batch <- Colon.atac$Sample
Colon.atac$Batch[Colon.atac$Batch == "B004-A-004-D-R2"] <- "B004-A-004-D"
table(Colon.atac$Batch)

# Plot UMAP
DimPlot(object = Colon.atac, group.by = "Celltype", label = TRUE) + NoLegend()
DimPlot(object = Colon.atac, group.by = "Sample", label = TRUE) + NoLegend()
DimPlot(object = Colon.atac, group.by = "Batch", label = TRUE) + NoLegend()

saveRDS(Colon.atac, paste(path, "/Signac/Colon_Healthy_Signac.rds", sep = ""))

#### Group major cell types ####
cols.celltype <- c("#B995E6", "#56B6A3", "#4AACF0", "#D188AF", "#F3A0DE", "#AFD655", "#8E1B7A", "#A458B0", "#EA6363", "#FFBC4F", "#FFBC4F")
cols.cellmajor <- c("#B995E6", "#4AACF0", "#D188AF", "#F3A0DE", "#EA6363", "#AFD655")
cols.lineage <- c( "#4AACF0", "#F3A0DE", "#AFD655")

# With Niko bins
Colon.atac <- readRDS(paste(path, "/Signac/Colon_Healthy_Signac.rds", sep = ""))
Colon.atac$Celltype_major <- Colon.atac$Celltype
Colon.atac$Celltype_major[Colon.atac$Celltype_major == "Immature Goblet"] <- "Goblet"
Colon.atac$Celltype_major[Colon.atac$Celltype_major %in% c("TA1", "TA2", "Secretory TA", "Enterocyte Progenitors")] <- "TA"
Colon.atac$Celltype_major[Colon.atac$Celltype_major %in% c("Immature Enterocytes")] <- "Enterocytes"

# Main celltypes
Colon.atac$Celltype_lineage <- Colon.atac$Celltype_major
Colon.atac$Celltype_lineage[Colon.atac$Celltype_lineage %in% c("Stem", "TA") ] <- "Undifferentiated"
Colon.atac$Celltype_lineage[Colon.atac$Celltype_lineage %in% c("Goblet", "Enteroendocrine") ] <- "Secretory"
Colon.atac$Celltype_lineage[Colon.atac$Celltype_lineage %in% c("Enterocytes", "Best4+ Enterocytes") ] <- "Absorptive"

DimPlot(object = Colon.atac,  group.by = "Celltype", cols = cols.celltype) +
  DimPlot(object = Colon.atac,  group.by = "Celltype_major", cols = cols.cellmajor) +
  DimPlot(object = Colon.atac,  group.by = "Celltype_lineage", cols = cols.lineage) 

#### Batch correction ####
# Integrate the gene activity matrix
DefaultAssay(Colon.atac)

# Merge two replicate samples, one is too tiny for separate integration
Colon.atac$Batch <- Colon.atac$Sample
Colon.atac$Batch[Colon.atac$Batch == "B004-A-004-D-R2"] <- "B004-A-004-D"
table(Colon.atac$Batch)
Colon.list <- SplitObject(Colon.atac, split.by = "Batch")

# normalize and identify variable features for each dataset independently
Colon.list <- lapply(X = Colon.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = Colon.list)
colon.anchors <- FindIntegrationAnchors(object.list = Colon.list, anchor.features = features)
colon.combined <- IntegrateData(anchorset = colon.anchors, verbose = TRUE)
DefaultAssay(colon.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
colon.combined <- ScaleData(colon.combined, verbose = TRUE)
colon.combined <- RunPCA(colon.combined, npcs = 30, verbose = TRUE)
colon.combined <- RunUMAP(colon.combined, reduction = "pca", dims = 1:30)

DimPlot(object = colon.combined,  group.by = "Celltype", cols = cols.celltype) +
  DimPlot(object = colon.combined,  group.by = "Celltype_major", cols = cols.cellmajor) +
  DimPlot(object = colon.combined,  group.by = "Celltype_lineage", cols = cols.lineage) 
