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
groups <- factor(paste(cellstage, treatment, sep="."))
design <- model.matrix(~0+groups)
colnames(design) <- levels(groups)

# estimate common, tagwise, and trend dispersion and plot the biological
# coefficient of variation (BCV)
data <- edgeR::estimateDisp(data, design, robust=TRUE)
plotBCV(data)

# construct contrasts (aka comparison vectors)
contrasts <- makeContrasts(ERGvsControl = ERG.24h-ERG.plus,
                           ERGvsMinus = ERG.minus - ERG.plus,
                           MRGvsControl = MRG.24h-MRG.plus,
                           MRGvsMinus = MRG.minus - MRG.plus,
                           NEvsControl = NE.24h-NE.plus,
                           NEvsMinus = NE.minus - NE.plus, levels=design2)

fit <- edgeR::glmFit(data, design = design)
result <- list()

comparisons <- colnames(contrasts)
for (val in comparisons){
  lrt <- edgeR::glmLRT(fit, contrast = contrasts[, val])
  sig.genes <- data.frame(topTags(lrt, n=100000, p.value = 0.05))

  by.gene <- group_by(sig.genes, Target.Gene.Symbol)
  same.as.majority.sign <- function(x) as.numeric(sign(x) == sign(mean(sign(x))))
  gene.info <- summarize(by.gene, count = n(), means = mean(sign(logFC)),
                         wmeanfc = sum(logFC*logCPM*same.as.majority.sign(logFC)/sum(logCPM*same.as.majority.sign(logFC))))
  gene.info <- filter(gene.info, count >= 2 & means != 0)
  gene.info <- gene.info[c("Target.Gene.Symbol", "wmeanfc")]
  result[[val]] <- gene.info
}

tophits <- list()

# The PlusvsControl and PlusvsMinus samples for each stage by taking the
# maximum abs value in either as the score on a gene by gene basis
for (i in seq(1, length(comparisons), 2)){
  combined <- rbind(result[[comparisons[i]]], result[[comparisons[i+1]]])
  abs.ordered <- combined[order(combined$Target.Gene.Symbol,
                                -abs(combined$wmeanfc)), ]
  take.max.abs <- abs.ordered[!duplicated(abs.ordered$Target.Gene.Symbol), ]
  tophits[[strsplit(comparisons[i], "vs")[[1]][1]]] <- take.max.abs
}

final.data <- merge(tophits$NE, tophits$ERG, by = "Target.Gene.Symbol", all = TRUE)
final.data <- merge(final.data, tophits$MRG, by = "Target.Gene.Symbol", all = TRUE)
colnames(final.data) <- c("Target.Gene.Symbol", "NE", "ERG", "MRG")
final.data[is.na(final.data)] <- 0