library(org.Hs.eg.db)
library(zoo)
library(readxl)
library(dplyr)
library(amap)

proc_data <- read_excel("data/nature13990-s4.xlsx", sheet = "All_stages_KD_0.05_results")
tf_data <- read.table("trans_fact.txt", sep="\t", stringsAsFactors = F,
                      fill = TRUE, quote="", header = TRUE, fileEncoding = "ISO-8859-2")
tf_data[tf_data == ""] <- NA
tf_data$Gene.ID <- na.locf(tf_data$Gene.ID)
tf_data$Description <- na.locf(tf_data$Description)
colnames(tf_data)[3] <- "Classification"

# GFP is a negative control
proc_data <- proc_data[proc_data$Gene != "GFP", ]

# Map Entrez ids to data
proc_data$entrez <- mget(proc_data$Gene,  org.Hs.egALIAS2EG )
# Pick only one id for each gene
proc_data$entrez <- lapply(proc_data$entrez, function (x) x[1])
proc_data <- transform(proc_data, entrez = as.numeric(entrez))

mydata <- proc_data[, 2:4]

wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(Kmeans(data, i, iter.max = 25)$withinss)}
  plot(1:nc, log10(wss), type="b", xlab="Number of Clusters",
       ylab="Log10 Within groups sum of squares")}

wssplot(mydata)


# K-Means Cluster Analysis
fit <- Kmeans(mydata, 8) # 8 cluster solution
proc_data$cluster <- fit$cluster
write.csv(proc_data, file="8clusters.csv")

clusters <- group_by(proc_data, cluster)

