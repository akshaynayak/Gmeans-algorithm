gmeans<-function(X=hw5.3d.data,significance_level=5){
library(MASS)
library(stats)
library(nortest)
library(fpc)
library(cluster)
library(jpeg)


if(significance_level==10)
critical_value<-0.631
if(significance_level==5)
  critical_value<-0.752
if(significance_level==2.5)
  critical_value<-0.873
if(significance_level==1)
  critical_value<-1.035
if(significance_level==0.5)
  critical_value<-1.159
if(significance_level==0.1)
  critical_value<-1.4607
if(significance_level==0.05)
  critical_value<-1.5903
if(significance_level==0.01)
  critical_value<-1.869

cat("\n critical value for this significance_level(in %) is:", critical_value)

#X<-hw5.3d.data
X<-as.matrix(X)
center<-colMeans(X)
plot(X)

center<-colMeans(X)
center<-matrix(center,ncol = ncol(X))

#***************  NOTE     **********
# ----- change path below to SAVE PLOTS generated during iterations 
#mypath <- file.path("C:","Users","AkshayN","Desktop","fds","r_codes","hw5_papers","cluster_plots",paste("cluster_generated_during_iteration_0.jpg", sep = ""))
#jpeg(file=mypath)
#mytitle = paste("clusters detected during iteration:0")
#plot(X,xlab = "dim1",ylab="dim2",main=mytitle)
#points(center,col="orange",pch=11,lwd=3)
#dev.off()
#**********************************
  
  
  
points(center,pch=11,col="red",lwd=3)
#for the first iteration we assume just one cluster with center as the mean of all the attributes and test if it follows gaussian distribution
prin_comp<-eigen(cov(X))
s=prin_comp$vectors[,1]
m=s*sqrt(2*prin_comp$values[1]/pi)

c1<-center+m
c2<-center-m
cat("\n Proposed centers(c1,c2) for cluster\n",c1,"\n",c2)

KMC = kmeans(X, centers=matrix(rbind(c1,c2), ncol = ncol(X)))
cat("\n centers found by using kmeans \n",KMC$centers[1,],"\n",KMC$centers[2,])
X<-as.matrix(X)
#str(KMC)
v<-KMC$centers[1,]-KMC$centers[2,]
normv<-sqrt(sum(v^2))
xdash<-X%*%cbind(v)/(normv^2)
scaled_xdash<-scale(xdash)
zxdash<-ecdf(scaled_xdash)
library(nortest)
test_result<-ad.test(zxdash(scaled_xdash))
cat("\n the test statistic is",test_result$statistic)
nsq<-nrow(X)^2
modified_statistic<-test_result$statistic*(1+(4/nrow(X))-(25/nsq))
cat("\nA*(Z)",modified_statistic)

if(test_result$statistic>critical_value)
{
  newcenters=rbind(KMC$centers[1,],KMC$centers[2,])
  cat("\n test statistic is greater than critical value we accept the split")
}else 
{
  print("\n the dataset just has one cluster")
  return
} 
iteration<-1
#while new centers are detected keep executing the while loop
while(nrow(newcenters)>nrow(center)){
  cat("\n-----------Next Iteration-------------",iteration)
  center<-newcenters
  newcenters<-{}
  KMC<-kmeans(X,centers = center)
  
  
  
  #*******************************
  #UNCOMMENT AND change the file path below to store the plots formed during the iterations
  #    WILL GIVE AN ERROR CHANGE PATH BELOW TO SAVE PLOTS
# mypath <- file.path("C:","Users","AkshayN","Desktop","fds","r_codes","hw5_papers","cluster_plots",paste("cluster_generated_during_iteration", iteration, ".jpg", sep = ""))
#  jpeg(file=mypath)
#  mytitle = paste("clusters detected during iteration:",iteration)
#  plot(X,col=KMC$cluster,xlab = "dim1",ylab="dim2",main=mytitle)
#  points(KMC$centers,col="orange",pch=11,lwd=3)
#  dev.off()
  #********************************
  
  
  
  for(i in 1:nrow(KMC$centers))
  {
    cat("\n\n--------considering cluster------",i)  
    ci<-KMC$centers[i,]
    cat("\n Current center for this cluster",ci)
    xi<-subset(X,KMC$cluster==i)
    prin_comp<-eigen(cov(xi))
    s=prin_comp$vectors[,1]
    m=s*sqrt(2*prin_comp$values[1]/pi)
    c1<-ci+m
    c2<-ci-m
    cat("\n proposed centers(c1,c2) for this cluster\n",c1,"\n",c2)
    kmcxi<-kmeans(xi,centers=matrix(rbind(c1,c2),ncol=ncol(xi)))
    cat("\n centers found by using kmeans \n",kmcxi$centers[1,],"\n",kmcxi$centers[2,])
    if(kmcxi$size[1]!=0&&kmcxi$size[2]!=0)
  {  v=kmcxi$centers[1,]-kmcxi$centers[2,]
    
    normv<-sqrt(sum(v^2))
    xdash<-xi%*%cbind(v)/(normv^2)
    scaled_xdash<-scale(xdash)
    zxdash<-ecdf(scaled_xdash)
    library(nortest)
    test_result<-ad.test(zxdash(scaled_xdash))
    cat("\n A(Z)",test_result$statistic)  
    nsq<-nrow(X)^2
    modified_statistic<-test_result$statistic*(1+(4/nrow(X))-(25/nsq))
    cat("\nA*(Z)",modified_statistic)
    
    #critical value for alpha=0.0001 is 1.8692
    if(modified_statistic>critical_value)
    {
      newcenters=rbind(newcenters,kmcxi$centers[1,],kmcxi$centers[2,])
      cat("\n since test statistic is greater than the critical value we reject H0 and keep{c1,c2}")
    }else 
    {
      newcenters=rbind(newcenters,ci)
      cat("\n since test statistic is smaller than cv we do not reject H0 and discard{c1,c2}")
    } 
  }else{
      if(kmcxi$size[1]!=0)
        modified_center<-kmcxi$centers[1,]
      else
        modified_center<-kmcxi$centers[2,]
    newcenters=rbind(newcenters,modified_center)
  }
      
  }
  iteration<-iteration+1
} 
cat("\n ------Final centers for the clusters-----\n")
for(i in 1:nrow(newcenters))
  cat("\n",newcenters[i,])

plot(X,col=KMC$cluster,xlab = "dim1",ylab="dim2")
points(KMC$centers,col="orange",pch=11)
clusplot(X,KMC$cluster,color = TRUE,shade = TRUE,labels=2,lines = 0)
}
