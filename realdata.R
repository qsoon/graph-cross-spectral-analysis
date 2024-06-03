library(gasper)
library(GNAR)
library(ggplot2)
library(cccd)
library(igraph)
library(igraphdata)
library(MASS)
library(pracma)
library(R.matlab)
library(geosphere)
library(grid)
library(gridBase)
library(gridExtra)

source("utils.R")


#################################
#### French meteorology data ####
#################################
###########################
data_frtemp <- readMat("Data/stationary/meteo_molene_t.mat")
fr.x.temp <- unlist(data_frtemp$info[[4]])
fr.y.temp <- unlist(data_frtemp$info[[3]])
fr.z.temp <- unlist(data_frtemp$info[[5]])

data_frhumid <- readMat("Data/stationary/meteo_molene_u.mat")
fr.x <- unlist(data_frhumid$info[[4]])
fr.y <- unlist(data_frhumid$info[[3]])
fr.z <- unlist(data_frhumid$info[[5]])

loc.ind.x <- c()
loc.ind.y <- c()
loc.ind.z <- c()
for(i in 1:22){
  loc.ind.x <- c(loc.ind.x, which(fr.x[i]==fr.x.temp))
  loc.ind.y <- c(loc.ind.y, which(fr.y[i]==fr.y.temp))
  loc.ind.z <- c(loc.ind.z, which(fr.z[i]==fr.z.temp))
}

data_frtemp$value <- data_frtemp$value[loc.ind.x,]
  

frmeteo <- list()
frmeteo$xyz <- cbind(fr.x, fr.y, 5*fr.z)
rownames(frmeteo$xyz) <- 1:nrow(frmeteo$xyz) 
frmeteo$temp <- data_frtemp$value
frmeteo$humid <- data_frhumid$value

frmeteo$temp.centered <- t(t(data_frtemp$value)- colMeans(data_frtemp$value))
frmeteo$humid.centered <- t(t(data_frhumid$value)- colMeans(data_frhumid$value))

plot(frmeteo$humid.centered[,1], type="l")
lines(frmeteo$temp.centered[,1], col="red")


distmat.fr <- as.matrix(dist(frmeteo$xyz, method="euclidean"))
A.fr <- c()
for(i in 1:(nrow(distmat.fr)-1)){
  for(j in (i+1):ncol(distmat.fr)){
    val <- distmat.fr[i,j]
    A.fr <- rbind(A.fr, c(i,j,val))
  }
}

# G.knn <- nng(dx=distmat.htemp, k=5, mutual=TRUE)
G.knn <- as.undirected(nng(dx=distmat.fr, k=6), mode="collapse")
edge.wt <- igraph::as_data_frame(G.knn, what="edges")
edge.wt <- sapply(edge.wt, as.numeric)
edge.wt <- cbind(edge.wt, 0)

for(i in 1:nrow(edge.wt)){
  edge.wt[i,3] <- distmat.fr[edge.wt[i,1], edge.wt[i,2]]
}  

wmat <- matrix(0, nrow=nrow(frmeteo$xyz), ncol=nrow(frmeteo$xyz))

colnames(wmat) <- 1:nrow(frmeteo$xyz) 
rownames(wmat) <- 1:nrow(frmeteo$xyz)

for(i in 1:nrow(edge.wt)){
  wmat[edge.wt[i,1], edge.wt[i,2]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
  wmat[edge.wt[i,2], edge.wt[i,1]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
}  

sp.wmat <- c()
for(i in 1:nrow(edge.wt)){
  sp.wmat <- rbind(sp.wmat, c(edge.wt[i,1], edge.wt[i,2], 
                              wmat[edge.wt[i,1], edge.wt[i,2]]))
}

# weight matrix
frmeteo$A <- wmat

# sparse weight matrix
frmeteo$sA <- sp.wmat

frmeteo$dist <- distmat.fr
frmeteo$sdist <- A.fr
frmeteo$xy <- frmeteo$xyz[,1:2]

# visualize
plot_graph(frmeteo)
plot_graph_custom3(frmeteo, e.size=1.3, v.size=6, vertex_color = frmeteo$temp.centered[,1], value="Temperature")

L.fr <- gasper::laplacian_mat(frmeteo$A)
N.fr <- nrow(frmeteo$xy)
val1 <- eigensort(L.fr)
evalues.fr <- val1$evalues
evectors.fr <- val1$evectors
# largest eigenvalue
lmax.fr <- max(evalues.fr)

# monthly avg
x <- rowMeans(frmeteo$temp.centered)
y <- rowMeans(frmeteo$humid.centered)

par(mfrow=c(1,1))
plot(x, type="l", main="temperature vs humidity", ylim=c(-10,10))
lines(y, col="red")
legend(legend=c("temperature", "humidity"), lty=c(1,1), col=c("black", "red"), "topright")


M.fr <- 40
cpsd.window.fr <- cpsd.graph(x=x, y=y, S=L.fr, g=g, M=M.fr, method="window")
psd.window.fr.x <- psd.graph(x=x, S=L.fr, g=g, M=M.fr, method="window")
psd.window.fr.y <- psd.graph(x=y, S=L.fr, g=g, M=M.fr, method="window")

# M is suitable?
(M.fr+1)*lmax.fr / M.fr^2 # sigma square
lmax.fr^2 / N.fr^2


tau.fr <- (M.fr+1)*lmax.fr / M.fr^2
par(mfrow=c(1,3), oma=c(0,1,0,0), mar=c(5,5,5,4)+0.1)
plot(tau.fr*c(1:M.fr), psd.window.fr.x, main="GPSD (Temperature)", pch=16, lwd=2,
     xlab="Graph frequency",ylab="Power spectral density", type="l", cex.lab=1.7, cex.axis=1.3, cex.main=1.7)
# abline(v=tau.fr*which.min(cpsd.window.fr[1:5]), lty=2, col="red", lwd=3)
# abline(v=tau.fr*which.max(psd.window.fr.y), lty=2, col="red", lwd=3)
abline(v=tau.fr*which.max(psd.window.fr.x), lty=2, col="red", lwd=3)
plot(tau.fr*c(1:M.fr), psd.window.fr.y, main="GPSD (Humidity)", pch=16, lwd=2,
     xlab="Graph frequency",ylab="Power spectral density", type="l", cex.lab=1.7, cex.axis=1.3, cex.main=1.7)
# abline(v=tau.fr*which.min(cpsd.window.fr[1:5]), lty=2, col="red", lwd=3)
abline(v=tau.fr*which.max(psd.window.fr.y), lty=2, col="blue", lwd=3)
# abline(v=tau.fr*which.max(psd.window.fr.x), lty=2, col="blue", lwd=3)
plot(tau.fr*c(1:M.fr), cpsd.window.fr, main="GCSD", pch=16, lwd=2,
     xlab="Graph frequency",ylab="Cross-spectral density", type="l", cex.lab=1.7, cex.axis=1.3, cex.main=1.7)
# abline(v=tau.fr*which.min(cpsd.window.fr[1:5]), lty=2, col="red", lwd=3)
abline(v=tau.fr*which.max(psd.window.fr.y), lty=2, col="blue", lwd=3)
abline(v=tau.fr*which.max(psd.window.fr.x), lty=2, col="red", lwd=3)


# visualize
g1 <- plot_graph_custom3(frmeteo, e.size=1.3, v.size=6, vertex_color = evectors.fr[,3], value="value", ratio=0.6,
                         min=-0.68, max=0.68)
g2 <- plot_graph_custom3(frmeteo, e.size=1.3, v.size=6, vertex_color = evectors.fr[,7], value="value", ratio=0.6,
                         min=-0.68, max=0.68)
g3 <- plot_graph_custom3(frmeteo, e.size=1.3, v.size=6, vertex_color = x, value="Temperature", ratio=0.6, 
                         min=-2, max=2)
g4 <- plot_graph_custom3(frmeteo, e.size=1.3, v.size=6, vertex_color = y, value="Humidity", ratio=0.6,
                         min=-5.1, max=5.1)
grid.arrange(g3,g4,g1,g2, nrow=2)
