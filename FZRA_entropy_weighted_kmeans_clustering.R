# Entropy-weighted k-means clustering
# Matt Irish
# Feb 2020
#
# Takes name of .MAT file containing k-means input matrix as command line argument

library(wskm) #installed from source with tarball from CRAN
library(R.matlab)

# Test ewkm:
myewkm <- ewkm(iris[1:4], 3, lambda=0.5, maxiter=100)
plot(iris[1:4], col=myewkm$cluster)
# For comparative testing
mykm <- kmeans(iris[1:4], 3)
plot(iris[1:4], col=mykm$cluster)

########################################################################################

setwd('/Users/mattirish/Documents/MATLAB/freezing_rain_trends_project')

# Import data from MATLAB ##############################################################
X_mat <- readMat("X_for_kmeans_612.mat")
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0){
  X_mat <- readMat(sprintf("%s.mat",args[1]))
  numclusters <- args[2]
}

# Find ideal number of clusters ########################################################
# Note: adapted from https://medium.com/codesmart/r-series-k-means-clustering-silhouette-794774b46586
silhouette_score <- function(k){
  km <- ewkm( t(X_mat[["X.scaled"]]), k, lambda=3, maxiter=100, delta=0.01, maxrestart=10)
  ss <- silhouette(km$cluster, dist(t(X_mat[["X.scaled"]])))
  mean(ss[, 3])
}
k <- 2:10
avg_sil <- sapply(k, silhouette_score)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)

numclusters = 6
#Run ewkm:
ewkm_result <- ewkm( t(X_mat[["X.scaled"]]), numclusters, lambda=3, maxiter=100, delta=0.01, maxrestart=10)


# Export data back to MATLAB ###########################################################
if(length(args)>0){
  writeMat(sprintf("%s_ewkm_results.mat",args[1]),
           IDX=ewkm_result$cluster,
           centroids=ewkm_result$centers,
           fixNames=T)
}else{
  writeMat("X_for_kmeans_612_ewkm_results.mat",
           IDX=ewkm_result$cluster,
           centroids=ewkm_result$centers,
           fixNames=T)
}

