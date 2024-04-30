library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(factoextra)

path <- "/mnt/data/Projects/IBD-CRC_WGS_Finland"

# Colors
cols.ibs <- c("#D76A03", "#32908F")

### Import data ####
# Import Patient data
Patient.dat <- read.csv(paste(path, "/Data/MSSCRC_IBDCRC_varcounts.csv.gz", sep = ""), row.names=1)
Patient.dat <- t(Patient.dat)

new.labels <- paste(rownames(Patient.dat), "000000", sep = "")
new.labels <- paste(new.labels, as.character(as.integer(as.numeric(sapply(strsplit(as.character(new.labels), "_"), "[", 2 )) + 1000000)), sep = "_")
new.labels <- gsub("_", "-", new.labels)

rownames(Patient.dat) <- new.labels
table(rownames(Patient.dat) %in% rownames(Average_dat))

IBDCRC_meta <- read.csv(paste(path, "/Data/MSSCRC_IBDCRC_meta.csv", sep = ""))
table(colnames(Patient.dat) == IBDCRC_meta$sn)
IBDCRC_meta$Label <- NA
IBDCRC_meta$Label[IBDCRC_meta$IBDstatus == "yes"] <- paste("IBD-CRC", 1:25, sep = "_")
IBDCRC_meta$Label[IBDCRC_meta$IBDstatus == "no"] <- paste("sCRC", 1:257, sep = "_")

# Replace colnames
table(IBDCRC_meta$sn == colnames(Patient.dat))
colnames(Patient.dat) <- IBDCRC_meta$Label

### Differential bins: IBD vs sCRC ####
# Normalize per patient
Patient.dat.norm <- t(apply(Patient.dat, 1, function(x) x/colMeans(Patient.dat)))

Differential_bins <- data.frame(GSALightning::wilcoxTest(eset = Patient.dat.norm, fac = sapply(strsplit( colnames(Patient.dat.norm), "_"), "[", 1 ), tests = "unpaired"))

cut.off <- 0.01
Differential_bins_sig <- Differential_bins %>% dplyr::filter(`p.value.up.regulated.in.IBD.CRC` <= cut.off | `p.value.up.regulated.in.sCRC` <= cut.off)

Differential_bins_sig$pct.IBD <- rowSums(Patient.dat.norm[rownames(Differential_bins_sig), 1:25] > 0)/25
Differential_bins_sig$pct.sCRC <- rowSums(Patient.dat.norm[rownames(Differential_bins_sig), 26:ncol(Patient.dat.norm)] > 0)/257
Differential_bins_sig$pct.delta <- Differential_bins_sig$pct.sCRC -Differential_bins_sig$pct.IBD

Differential_bins_sig$Direction_up <- "sCRC"
Differential_bins_sig$Direction_up[Differential_bins_sig$p.value.up.regulated.in.IBD.CRC <= cut.off] <- "IBD-CRC"
table(Differential_bins_sig$Direction_up)

Differential_bins_sig$bin <- rownames(Differential_bins_sig)
Differential_bins_sig$Chr <- sapply(strsplit( Differential_bins_sig$bin, "-"), "[", 1 )

Differential_bins_sig %>%
  group_by(Direction_up, Chr) %>%
  summarise(N = n()) %>%
  mutate(N = ifelse(Direction_up == "IBD-CRC", N*-1, N)) %>%
  mutate(Chr = forcats::fct_reorder(Chr,Direction_up, .desc = T)) %>%
  mutate(Chr = forcats::fct_reorder(Chr,N, .desc = F)) %>%
  ggplot(aes(x = Chr, y = N, fill = Direction_up)) + geom_bar(stat = "identity", position = "stack", width = 0.7) +
  theme_classic()  + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_manual("Group",values = cols.ibs) + xlab("") + ylab("Differential Genomic Bins") + coord_flip() +
  geom_hline(yintercept = 0, color = "black") +
  ggtitle("Wilcox Test, P < 0.01") +  theme(plot.title = element_text(hjust = 0.5))

Differential_bins_sig_qval <- Differential_bins_sig %>% dplyr::filter(`q.value.up.regulated.in.IBD.CRC` <= cut.off | `q.value.up.regulated.in.sCRC` <= cut.off)

### Average Mut /Mb in regions ####
Mut_dat <- data.frame(Group = sapply(strsplit( IBDCRC_meta$Label, "_"), "[", 1 ),
                      Whole_genome = colMeans(Patient.dat),
                      IBD_specific = colMeans(Patient.dat.norm[rownames(Differential_bins_sig[Differential_bins_sig$Direction_up == "IBD-CRC",]),]),
                      sCRC_specific = colMeans(Patient.dat.norm[rownames(Differential_bins_sig[Differential_bins_sig$Direction_up == "sCRC",]),]))

Mut_dat %>% reshape2::melt(id.vars = "Group") %>%
  ggplot(aes(x = Group, y = value, color = Group)) + 
  geom_violin(trim = FALSE, aes(color = Group), position = position_dodge(width = 0.8), width = 0.75) + 
  ggbeeswarm::geom_quasirandom(method = "pseudorandom", size = 0.2, dodge.width=.2, varwidth = T) +
  facet_wrap(~variable, scales = "free") + xlab("") + ylab("Average Mutations per Mb") + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(strip.background = element_blank()) + 
  scale_color_manual(values = cols.ibs) + RotatedAxis() + NoLegend() +
  ggpubr::stat_compare_means(method = "t.test", label = "p.format", label.x = 1.25) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 0.2, size = 0.2, color = "black")

