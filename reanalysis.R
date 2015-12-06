library(readxl)
library(dplyr)
library(edgeR)
library(gplots)

raw_data <- read_excel("data/nature13990-s4.xlsx", sheet = "Raw shRNA counts")
paper.data <- read_excel("data/nature13990-s4.xlsx", sheet = "All_stages_KD_0.05_results")

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

# data <- edgeR::calcNormFactors(data)

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
                           NEvsMinus = NE.minus - NE.plus, levels=design)

fit <- edgeR::glmFit(data, design = design)
result <- list()

comparisons <- colnames(contrasts)
for (val in comparisons){
  lrt <- edgeR::glmLRT(fit, contrast = contrasts[, val])
  sig.genes <- data.frame(topTags(lrt, n=100000, p.value = 0.05))

  by.gene <- group_by(sig.genes, Target.Gene.Symbol)
  same.as.majority.sign <- function(x) as.numeric(sign(x) == sign(mean(sign(x))))
  gene.info <- summarize(by.gene, count = n(), means = mean(sign(logFC)),
                         # wmeanfc = sum(logFC*logCPM*same.as.majority.sign(logFC)/sum(logCPM*same.as.majority.sign(logFC))))
                         wmeanfc = sum(logFC*logCPM/sum(logCPM)))
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
colnames(final.data) <- c("Gene", "NE", "ERG", "MRG")
final.data[is.na(final.data)] <- 0

final.mat <- as.matrix(final.data[, 2:4])
rownames(final.mat) <- final.data[, 1]
# drop OLIG1 because the value is funky, look into this TODO
final.mat <- final.mat[rownames(final.mat) != "OLIG1", ]
final.mat <- final.mat[order(rownames(final.mat)), ]

heatmap.2(final.mat[rowSums(final.mat > 0) >= 1, ], col = greenred(75), Colv = F,
          dendrogram = "none", density.info = "none", key = TRUE, trace = "none",
          main="Reproduced Figure", lhei = c(1, 5), margins = c(4, 5), Rowv = F,
          key.xlab = "Depletion Score", key.title = NA, cexRow=0.5, cexCol = 1)

paper.mat <- as.matrix(paper.data[, 2:4])
rownames(paper.mat) <- paper.data$Gene
paper.mat <- paper.mat[order(rownames(paper.mat)), ]

heatmap.2(paper.mat[rowSums(paper.mat > 0) >= 1, ], col = greenred(75), Colv = F,
          dendrogram = "none", density.info = "none", key = TRUE, trace = "none",
          main="Data from paper", lhei = c(1, 5), margins = c(4, 5), Rowv = F,
          key.xlab = "Depletion Score", key.title = NA, cexRow=0.5, cexCol = 1)

shared <- intersect(rownames(paper.mat), rownames(final.mat))