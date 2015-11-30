library(readxl)
library(dplyr)
library(edgeR)
library(ggplot2)

raw_data <- read_excel("data/nature13990-s4.xlsx", sheet = "Raw shRNA counts")

cpms <- edgeR::cpm(raw_data[,-(1:3)])
log_cpms <- edgeR::cpm(raw_data[,-(1:3)], log = TRUE, prior.count = 0.5)
rownames(cpms) <- raw_data$`Target Gene Symbol`
rownames(log_cpms) <- raw_data$`Target Gene Symbol`

# Select only hairpins with CPM > 0.5 in more than two samples
cpm_mask <- rowSums(cpms > 0.5) >= 2
cpms <- cpms[cpm_mask, ]

log_cors <- cor(log_cpms)
heatmap(log_cors, Rowv = NA, Colv = NA, scale = "none")
