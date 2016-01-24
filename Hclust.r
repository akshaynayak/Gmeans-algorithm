
##################
library(mclust)
# Run the function to see how many clusters
# it finds to be optimal, set it to search for
# at least 1 model and up 20.
d_clust <- Mclust(as.matrix(X), G=1:20)
m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")
# 4 clusters
plot(d_clust)

###############
X<-hw5.3d.data
kosDist = dist(X, method="euclidean")

kosHierClust = hclust(kosDist, method="ward.D")
plot(kosHierClust)

rect.hclust(kosHierClust, k = 12, border = "red")
kosClusters = cutree(kosHierClust, k = 15)
str(kosClusters)

plot(kosClusters)

