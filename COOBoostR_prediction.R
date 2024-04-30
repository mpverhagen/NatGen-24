library(ggplot2)
library(dplyr)
library(reshape2)
library(Seurat)
library(Signac)

### Installation ####
rm(list=ls())
# Set the working directory as COOBoostR-main/COOBoostR/ for package installation.
packagePath <- "/mnt/data/Projects/IBD-CRC_WGS_Finland/COOBoostR-main/COOBoostR-main/COOBoostR"
setwd(packagePath)
set.seed(1708)

# Installing packages using devtools and load COOBoostR.
devtools::install()
library(COOBoostR)

### Prepare data ####
# Import scATAC data set
Colon.atac <- readRDS("/mnt/data/User/Mathijs/Becker_polyps_CRC_singlecell/Signac/Colon_Healthy_1Mb_bins_Niko_Signac.rds")
Colon.sum <- AggregateExpression(Colon.atac, group.by = "Celltype", slot = "counts")
Colon.sum <- Colon.sum$bins

# Normalize matrix to equalize at 10 Million counts per cell type
Colon.sum <- sweep(Colon.sum, MARGIN=2, c(1e+07/colSums(Colon.sum)), `*`)
Colon.sum <- data.frame(round(Colon.sum))
colnames(Colon.sum)[1] <- "Best4.Enterocytes"

# Import Patient data
Patient.dat <- read.csv("/mnt/data/Projects/IBD-CRC_WGS_Finland/Data/MSSCRC_IBDCRC_varcounts.csv.gz", row.names=1)
Patient.dat <- t(Patient.dat)

new.labels <- paste(rownames(Patient.dat), "000000", sep = "")
new.labels <- paste(new.labels, as.character(as.integer(as.numeric(sapply(strsplit(as.character(new.labels), "_"), "[", 2 )) + 1000000)), sep = "_")
new.labels <- gsub("_", "-", new.labels)

rownames(Patient.dat) <- new.labels
table(rownames(Patient.dat) %in% rownames(Colon.sum))

IBDCRC_meta <- read.csv("/mnt/data/Projects/IBD-CRC_WGS_Finland/Data/MSSCRC_IBDCRC_meta.csv")
table(colnames(Patient.dat) == IBDCRC_meta$sn)
IBDCRC_meta$Label <- NA
IBDCRC_meta$Label[IBDCRC_meta$IBDstatus == "yes"] <- paste("IBD-CRC", 1:25, sep = "_")
IBDCRC_meta$Label[IBDCRC_meta$IBDstatus == "no"] <- paste("sCRC", 1:257, sep = "_")

# Replace colnames
table(IBDCRC_meta$sn == colnames(Patient.dat))
colnames(Patient.dat) <- IBDCRC_meta$Label
Patient.dat <- Patient.dat[match(rownames(Colon.sum), rownames(Patient.dat)),]
table(rownames(Patient.dat) == rownames(Colon.sum))

write.csv(Colon.sum, "/mnt/data/Projects/IBD-CRC_WGS_Finland/COOBoostR-results/2_Patientdata_default_10Aug/Data/COO_epimarker.csv", row.names = F)
write.csv(Patient.dat, "/mnt/data/Projects/IBD-CRC_WGS_Finland/COOBoostR-results/2_Patientdata_default_10Aug/Data/COO_mutation.csv", row.names = F)

### Cell-of-Origin predictions ####
#### COOBoost: Celltype major, mEta = 0.3, mdepth = 6 ####
library(Seurat)
library(Signac)

setwd("/mnt/data/Projects/IBD-CRC_WGS_Finland/COOBoostR-results/8_Patientdata_Celltype_major_learning0.3_depth6")
# sourcePath : Directory with data needed for analysis(epimarker, mutation). Must end with a delimiter (/).
sourcePath <- "/mnt/data/Projects/IBD-CRC_WGS_Finland/COOBoostR-results/6_Patientdata_Celltype_major_Lancet_parameters/Data/"
# resultsPath : The directory where the results will be saved after analysis. Must end with a delimiter (/).
resultsPath <- "/mnt/data/Projects/IBD-CRC_WGS_Finland/COOBoostR-results/8_Patientdata_Celltype_major_learning0.3_depth6/Results/" 

# Data generated through 1mb preprocessing. Requires csv format.(You need a csv file name, not a data frame or variable.)
epimarker <- "COO_epimarker.csv"
mutation <- "COO_mutation.csv"

# Analysis with COOBoostR
COOBoostR(sourcePath = sourcePath, resultPath = resultsPath, epimarker_rawdata = epimarker, mutation_rawdata = mutation, mEta = 0.3, mdepth = 6)

#### Inspect Results ####
# Colors
cols.celltype <- c("#B995E6", "#56B6A3", "#4AACF0", "#D188AF", "#F3A0DE", "#AFD655", "#8E1B7A", "#A458B0", "#EA6363", "#FF8427", "#FFBC4F")
cols.cellmajor <- c("#B995E6", "#4AACF0", "#D188AF", "#F3A0DE", "#EA6363", "#AFD655")
cols.lineage <- c( "#4AACF0", "#F3A0DE", "#AFD655")
cols.ibs <- c("#D76A03", "#32908F")

resultsPath <- "/mnt/data/Projects/IBD-CRC_WGS_Finland/COOBoostR-results/8_Patientdata_Celltype_major_learning0.3_depth6/Results/COO_mutation/" 
files <- list.files(resultsPath, pattern = "output_")
files.take <- files[-grep("extract",files)]

# Make loop to summarize data
data.out <- data.frame()

for (i in 1:length(files.take)){
  print(i)
  temp <- read.csv(paste(resultsPath, files.take[i], sep = ""), sep="")
  data.temp <- data.frame(Sample = gsub(".txt", "",paste(sapply(strsplit(files.take[i], "_",), "[", 2), sapply(strsplit(files.take[i], "_",), "[", 3), sep = "_")),Classification = temp$name)
  data.out <- rbind(data.out, data.temp)
  rm(temp, data.temp)}

data.out$Sample <- gsub("IBD.CRC", "IBD-CRC", data.out$Sample)
data.out$Group <- sapply(strsplit(as.character(data.out$Sample), "_"), "[", 1 )

table(data.out$Classification, data.out$Group)
test <- fisher.test(table(data.out$Classification, data.out$Group))
test

data.out %>% dplyr::group_by(Group, Classification) %>% summarize(N = n()) %>% group_by(Group) %>% mutate(freq = N/sum(N)) %>%
  ggplot(aes(x = Classification, y = N, fill = Classification)) + geom_bar(stat = "identity", width = 0.7) + 
  facet_wrap(~ Group, scales = "free_y") + 
  xlab("") + scale_fill_manual(values = cols.cellmajor[c(1:2,4:6)]) + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  theme(strip.background = element_blank()) + Seurat::NoLegend() + Seurat::RotatedAxis() +
  ylab("COOBoost predicted Cell-of-Origin")+ 
  geom_text(aes(label=paste0(round(freq*100,0),"%"),y=ifelse(N > 10, N+5, N+0.5)), size=4) +
  ggtitle(paste("Fisher's exact test, P =", round(test$p.value, 2), sep = " ")) + 
  theme(plot.title = element_text(hjust = 0.5)) 

# write data
write.csv(data.out, "/mnt/data/Projects/IBD-CRC_WGS_Finland/COOBoostR-results/8_Patientdata_Celltype_major_learning0.3_depth6/Results/Patient_COO_prediction.csv")