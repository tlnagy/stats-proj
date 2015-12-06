library(org.Hs.eg.db)
library(zoo)


#tf_data <- read.table("trans_fact.txt",sep="\t")
library(readxl)
proc_data <- read_excel("data/nature13990-s4.xlsx", sheet = "All_stages_KD_0.05_results")
tf_data <- read.table("trans_fact.txt", sep="\t", stringsAsFactors = F,
                      fill = TRUE, quote="", header = TRUE, fileEncoding = "ISO-8859-2")
tf_data[tf_data == ""] <- NA
tf_data$Gene.ID <- na.locf(tf_data$Gene.ID)
tf_data$Description <- na.locf(tf_data$Description)
View(tf_data)




 # Convert the object to a list
xx <- as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]
if(length(xx) > 0){
  # The entrez gene identifiers for the first two elements of XX
  xx[1:2]
  # Get the first one
  xx[[1]]
}

proc_data$entrez <- lapply(proc_data$Gene, function(x) xx[[x]])


mydata <- proc_data

# Determine number of clusters
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
fit <- kmeans(mydata, 5) # 5 cluster solution
# get cluster means 
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster)

# Ward Hierarchical Clustering
d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward") 
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")

# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit <- pvclust(mydata, method.hclust="ward",
               method.dist="euclidean")
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)

# Model Based Clustering
library(mclust)
fit <- Mclust(mydata)
plot(fit) # plot results 
summary(fit) # display the best model

# K-Means Clustering with 5 clusters
fit <- kmeans(mydata, 5)

# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
library(cluster) 
clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)

# Centroid Plot against 1st 2 discriminant functions
library(fpc)
plotcluster(mydata, fit$cluster)
