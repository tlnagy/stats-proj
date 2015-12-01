library(readxl)
library(dplyr)
library(edgeR)
library(ggplot2)

raw_data <- read_excel("data/nature13990-s4.xlsx", sheet = "Raw shRNA counts")

# experiment design info
exp.info <- matrix(unlist(strsplit(colnames(raw_data[,-(1:3)]), "_")), ncol=4, byrow = TRUE)
cellstage <- factor(exp.info[,1])
treatment <- factor(exp.info[, 3])
treatment <- relevel(treatment, ref = "24h")
timepoints <- factor(exp.info[,4])

groups <- factor(paste(cellstage, treatment, sep="."))
data <- DGEList(counts = raw_data[,-(1:3)], genes = raw_data[,1:3], group = groups)

# Select only hairpins with CPM > 0.5 in more than two samples
shrnas.to.keep <- rowSums(edgeR::cpm(data) > 0.5) >= 2
data <- data[shrnas.to.keep, , keep.lib.sizes=FALSE]

data <- edgeR::calcNormFactors(data)

# construct design matrix
design <- model.matrix(~treatment + treatment:cellstage)
rownames(design) <- colnames(data)

# estimate common, tagwise, and trend dispersion
data <- edgeR::estimateDisp(data, design, robust=TRUE)
plotBCV(data)

