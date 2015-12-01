library(readxl)
library(dplyr)
library(edgeR)
library(ggplot2)

raw_data <- read_excel("data/nature13990-s4.xlsx", sheet = "Raw shRNA counts")

# experiment design info
exp.info <- matrix(unlist(strsplit(colnames(raw_data[,-(1:3)]), "_")),
                   ncol=4, byrow = TRUE)
cellstage <- factor(exp.info[,1])
treatment <- factor(exp.info[, 3])
treatment <- relevel(treatment, ref = "24h")
timepoints <- factor(exp.info[,4])

groups <- factor(paste(cellstage, treatment, sep="."))
data <- DGEList(counts = raw_data[,-(1:3)], genes = raw_data[,1:3],
                group = groups)

# Select only hairpins with CPM > 0.5 in more than two samples
shrnas.to.keep <- rowSums(edgeR::cpm(data) > 0.5) >= 2
data <- data[shrnas.to.keep, , keep.lib.sizes=FALSE]

data <- edgeR::calcNormFactors(data)

# attempt to reproduce extended data figure 3f
col.order <- c("NE_sh_minus_4", "NE_sh_minus_3", "NE_sh_minus_2",
               "ERG_sh_minus_3", "ERG_sh_minus_1", "ERG_sh_minus_2",
               "NE_sh_minus_1", "ERG_sh_plus_3", "ERG_sh_plus_2",
               "ERG_sh_plus_1", "NE_sh_plus_2", "NE_sh_plus_1", "NE_sh_plus_4",
               "NE_sh_plus_3", "MRG_sh_minus_2", "MRG_sh_minus_1",
               "MRG_sh_minus_3", "MRG_sh_minus_4", "MRG_sh_plus_2",
               "MRG_sh_plus_1", "MRG_sh_plus_3", "MRG_sh_plus_4",
               "ERG_sh_24h_2", "NE_sh_24h_2", "ERG_sh_24h_3", "NE_sh_24h_3",
               "NE_sh_24h_4", "MRG_sh_24h_3", "MRG_sh_24h_1", "MRG_sh_24h_4",
               "MRG_sh_24h_2", "ERG_sh_24h_1", "NE_sh_24h_1")
log.cpms <- edgeR::cpm(data, log=TRUE, prior.count = 0.25)
log.cpms <- log.cpms[ ,rev(col.order)]
heatmap(cor(log.cpms), Rowv = NA, Colv = NA, scale = "none")

# construct design matrix
design <- model.matrix(~treatment + treatment:cellstage)
rownames(design) <- colnames(data)

# estimate common, tagwise, and trend dispersion and plot the biological
# coefficient of variation (BCV)
data <- edgeR::estimateDisp(data, design, robust=TRUE)
plotBCV(data)

