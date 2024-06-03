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
library(expm)
library(Hmisc)

source("utils.R")

# window filterbank mother function for WGF-based estimation
g <- function(lambda, sigma.sq, m, tau){
  return(exp(-(lambda-m*tau)^2 / sigma.sq))
}

#############################
### 01. graph preparation ###
#############################


#########################
#### seoulmetro data ####
#########################

station.info = read.csv("Data/seoulmetro/master_station_info.csv")
station.distance = read.csv("Data/seoulmetro/master_station_distance.csv")

##############################
#### make seoulmetro data ####
##############################
seoulmetro <- list()

target_station <- unique(station.info$name) # total 243 stations
target_station_kor <- unique(station.info$name_kor)

# location
seoulmetro$xy <- sapply(station.info[!duplicated(station.info[,c("lon","lat")]),
                                     c("lon","lat")], as.numeric)
rownames(seoulmetro$xy) <- target_station


#########################################
#### edge weight matrix construction ####
#########################################
e.weight <- matrix(0,nrow=length(target_station), ncol=length(target_station))
colnames(e.weight) <- target_station
rownames(e.weight) <- target_station

for(i in 1:8){
  tmp <-station.distance[station.distance$line==i,]
  if(i==2){
    n <- 44 # circular line. 54th line : City Hall again
    e.weight["Seongsu", "Yongdap"] <- tmp$btwdist[45]
    e.weight["Yongdap", "Seongsu"] <- tmp$btwdist[45]
    for(j in (n+1):47){
      e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
      e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
    }
    e.weight["Sindorim", "Dorimcheon"] <- tmp$btwdist[49]
    e.weight["Dorimcheon", "Sindorim"] <- tmp$btwdist[49]
    for(j in 49:51){
      e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
      e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
    }
  } else if(i==5){
    n <- 46 # Hanam Pungsan is the end station of one line in Line5
    e.weight["Gangdong", "Dunchon-dong"] <- tmp$btwdist[47]
    e.weight["Dunchon-dong", "Gangdong"] <- tmp$btwdist[47]
    for(j in (n+1):52){
      e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
      e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
    }
  } else{
    n <- nrow(tmp)
  }
  for(j in 1:(n-1)){
    e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
    e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
  }
}

e.weight.old <- e.weight
e.weight.old[e.weight.old!=0] <- exp(-e.weight.old[e.weight.old!=0])
e.weight[e.weight!=0] <- exp(-(e.weight[e.weight!=0])^2 / mean(e.weight[e.weight!=0])^2)

e.sp.weight <- NULL
e.color <- c() # for line color
color.cand <- c("blue", "yellowgreen", "orangered", "cyan",
                "darkorchid", "chocolate3", "darkolivegreen", "hotpink")
for(i in 1:8){
  tmp <-station.distance[station.distance$line==i,]
  if(i==2){
    n <- 44 # circular line. 54th line : City Hall again
    e.sp.weight <- rbind(e.sp.weight, c("Seongsu", "Yongdap", tmp$btwdist[45]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight <- rbind(e.sp.weight, c("Yongdap", "Seongsu", tmp$btwdist[45]))
    e.color <- c(e.color, color.cand[i])
    for(j in (n+1):47){
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
    }
    e.sp.weight <- rbind(e.sp.weight, c("Sindorim", "Dorimcheon", tmp$btwdist[49]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight <- rbind(e.sp.weight, c("Dorimcheon", "Sindorim", tmp$btwdist[49]))
    e.color <- c(e.color, color.cand[i])
    
    for(j in 49:51){
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
    }
  } else if(i==5){
    n <- 46 # Hanam Pungsan is the end station of one line in Line5
    e.sp.weight<- rbind(e.sp.weight, c("Gangdong", "Dunchon-dong", tmp$btwdist[47]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight<- rbind(e.sp.weight, c("Dunchon-dong", "Gangdong", tmp$btwdist[47]))
    e.color <- c(e.color, color.cand[i])
    for(j in (n+1):52){
      e.sp.weight <- rbind(e.sp.weight, c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
      e.sp.weight <- rbind(e.sp.weight, c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
    }
  } else{
    n <- nrow(tmp)
  }
  for(j in 1:(n-1)){
    e.sp.weight <- rbind(e.sp.weight, 
                         c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight <- rbind(e.sp.weight, 
                         c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
    e.color <- c(e.color, color.cand[i])
  }
}

tmp <- as.data.frame(e.sp.weight)
tmp$V1 <- as.character(tmp$V1)
tmp$V2 <- as.character(tmp$V2)
tmp2 <- exp(-as.numeric(as.character(tmp$V3)))
tmp$V3 <- exp(-(as.numeric(as.character(tmp$V3)))^2 / (mean(as.numeric(as.character(tmp$V3))))^2)

e.sp.weight <- tmp
e.sp.weight.old <- e.sp.weight
e.sp.weight.old[,3] <- tmp2


# weight matrix
seoulmetro$A <- e.weight

# sparse weight matrix
seoulmetro$sA <- e.sp.weight

# dist matrix
distmat <- e.weight.old
distmat[distmat!=0] <- -log(e.weight.old[e.weight.old!=0])

tmp <- e.sp.weight.old
tmp[,3] <- -log(tmp[,3])

seoulmetro$dist <- distmat
seoulmetro$sdist <- tmp

plot_graph(seoulmetro)


L.metro <- laplacian_mat(seoulmetro$A) # laplacian matrix
val1 <- eigensort(L.metro)
evalues.metro <- val1$evalues
evectors.metro <- val1$evectors
# largest eigenvalue
lmax.metro <- max(evalues.metro)

N.metro <- nrow(L.metro)

#############################
### random sensor network ###
#############################
n <- 20
x <- rep(0:(n-1), n) # define coordinates of vertices
y <- rep(0:(n-1), rep(n,n))

set.seed(1)
x <- runif(n^2, min = 0, max = n)
y <- runif(n^2, min = 0, max = n)

irregular01 <- list()
irregular01$xy <- data.frame(x,y)

N.irregular01 <- nrow(irregular01$xy)
rownames(irregular01$xy) <- 1:N.irregular01

distmat.irregular01 <- as.matrix(dist(irregular01$xy))
A.irregular01 <- c()
for(i in 1:(nrow(distmat.irregular01)-1)){
  for(j in (i+1):ncol(distmat.irregular01)){
    val <- distmat.irregular01[i,j]
    A.irregular01 <- rbind(A.irregular01, c(i,j,val))
  }
}

# G.knn <- nng(dx=distmat.htemp, k=5, mutual=TRUE)
G.knn <- as.undirected(nng(dx=distmat.irregular01, k=5), mode="collapse")
edge.wt <- igraph::as_data_frame(G.knn, what="edges")
edge.wt <- sapply(edge.wt, as.numeric)
edge.wt <- cbind(edge.wt, 0)

for(i in 1:nrow(edge.wt)){
  edge.wt[i,3] <- distmat.irregular01[edge.wt[i,1], edge.wt[i,2]]
}  

wmat <- matrix(0, nrow=N.irregular01, ncol=N.irregular01)

colnames(wmat) <- 1:N.irregular01
rownames(wmat) <- 1:N.irregular01

for(i in 1:nrow(edge.wt)){
  wmat[edge.wt[i,1], edge.wt[i,2]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
  wmat[edge.wt[i,2], edge.wt[i,1]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
}  

sp.wmat <- c()
for(i in 1:nrow(edge.wt)){
  sp.wmat <- rbind(sp.wmat, c(edge.wt[i,1], edge.wt[i,2], 
                              wmat[edge.wt[i,1], edge.wt[i,2]]))
}


# sparse weight matrix
# weight matrix
irregular01$A <- wmat

# sparse weight matrix
irregular01$sA <- sp.wmat

irregular01$dist <- distmat.irregular01
irregular01$sdist <- A.irregular01

L.irregular01 <- gasper::laplacian_mat(irregular01$A)

plot_graph(irregular01)

val <- eigensort(L.irregular01)
lmax.irregular01 <- max(val$evalues)
evalues.irregular01 <- val$evalues
evectors.irregular01 <- val$evectors


###########################
### karate club network ###
###########################
data(karate)
V(karate)$name
E(karate)
edge.attributes(karate)
degree(karate)


N.karate <- length(V(karate)$name)
edge.wt <- igraph::as_data_frame(karate, what="edges")
for(i in 1:nrow(edge.wt)){
  edge.wt[i,1] <- which(V(karate)$name == edge.wt[i,1])
  edge.wt[i,2] <- which(V(karate)$name == edge.wt[i,2])
}
edge.wt <- sapply(edge.wt, as.numeric)
wmat.karate <- matrix(0, nrow=gorder(karate), ncol=gorder(karate))
# colnames(wmat.karate) <- V(karate)$name
# rownames(wmat.karate) <- V(karate)$name

colnames(wmat.karate) <- 1:N.karate
rownames(wmat.karate) <- 1:N.karate

for(i in 1:nrow(edge.wt)){
  wmat.karate[edge.wt[i,1], edge.wt[i,2]] <- edge.wt[i,3]
  wmat.karate[edge.wt[i,2], edge.wt[i,1]] <- edge.wt[i,3]
}

sp.wmat.karate <- c()
for(i in 1:nrow(edge.wt)){
  sp.wmat.karate <- rbind(sp.wmat.karate, c(edge.wt[i,1], edge.wt[i,2], 
                              wmat.karate[edge.wt[i,1], edge.wt[i,2]]))
}

zachary.karate <- list()
set.seed(1)
layout.karate <- layout.fruchterman.reingold(karate)
layout.karate <- layout.norm(layout.karate, -1, 1, -1, 1)
zachary.karate$xy <- data.frame(x=layout.karate[,1],
                                y=layout.karate[,2])

zachary.karate$A <- wmat.karate
zachary.karate$sA <- sp.wmat.karate

plot_graph(zachary.karate)

L.karate <- gasper::laplacian_mat(wmat.karate)
val <- eigensort(L.karate)
lmax.karate <- max(val$evalues)
evalues.karate <- val$evalues
evectors.karate <- val$evectors

eigenres.karate <- eigen(L.karate)


##############################
### Minnesota road network ###
##############################
Minnesota.network <- list()
Minnesota.network$xy <- data.frame(x=minnesota$xy[,1],
                                   y=minnesota$xy[,2])

Minnesota.network$A <- minnesota$A
Minnesota.network$A[Minnesota.network$A==2] <- 1
Minnesota.network$sA <- minnesota$sA
Minnesota.network$sA[Minnesota.network$sA[,3]==2,3] <- 1

Minnesota.network2 <- Minnesota.network
Minnesota.network2$xy <- Minnesota.network$xy[order(x),]
Minnesota.network2$A <- Minnesota.network$A[order(x),order(x)]
Minnesota.network2$sA <- Minnesota.network$sA
for(i in 1:nrow(Minnesota.network2$sA)){
  Minnesota.network2$sA[i,1] <- which(order(x) == Minnesota.network2$sA[i,1])
  Minnesota.network2$sA[i,2] <- which(order(x) == Minnesota.network2$sA[i,2])
}

plot_graph(Minnesota.network)
L.minnesota <- gasper::laplacian_mat(minnesota$A)
val <- eigensort(L.minnesota)
lmax.minnesota <- max(val$evalues)
evalues.minnesota <- val$evalues
evectors.minnesota <- val$evectors
N.minnesota <- nrow(L.minnesota)


##################################
### 02. estimation performance ###
##################################

# filter kernel들 polynomial이어야하지만, chebyshev polynomial로 근사 잘 되므로 okay.
filter_kernel <- function(lambda, lmax, tau1, tau2){
  return(sin(3*tau1*lambda/lmax)*exp(-tau2*lambda/lmax))
}

filter_kernel_highpass <- function(lambda, lmax, tau){
  return(lambda/lmax*exp(-tau*lambda/lmax))
}

# filter_kernel_lowpass <- function(lambda){
#   return(exp(-lambda^4))
# }

filter_kernel_heat <- function(lambda, lmax, tau){
  return(exp(-tau*lambda/lmax))
}

filter_kernel_Mexican <- function(lambda, lmax, tau1=5, tau2=25){
  return(tau1*lambda/lmax * exp(-tau2*lambda^2/lmax^2))
}
filter_poly <- function(lambda, h){
  res <- c()
  for(i in lambda){
    res <- c(res, sum(h*(i^c(0:(length(h)-1)))))
  }
  return(res)
}

tmp.heat <- polyApprox(function(l) filter_kernel_heat(l, lmax=lmax.karate, tau=10), 0, lmax.karate, 5)
tmp.Mexican <- polyApprox(function(l) filter_kernel_Mexican(l, lmax=lmax.karate,  tau1=5, tau2=25), 0, lmax.karate, 7)
tmp.sin.exp <- polyApprox(function(l) filter_kernel(l, lmax=lmax.karate, tau1=5, tau2=5), 0, lmax.karate, 9)
tmp.highpass <- polyApprox(function(l) filter_kernel_highpass(l, lmax=lmax.karate, tau=1), 0, lmax.karate, 4)

plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.heat$p)), type="l")
plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.Mexican$p)), type="l")
plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.sin.exp$p)), type="l")
plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.highpass$p)), type="l")

# R realizations case
R.list <- c(10,100,1000)
set.seed(1)
white.inputs <- t(mvrnorm(n=max(R.list), mu=rep(0,N.karate), Sigma=diag(1,N.karate))) # N.karate x R matrix
H.heat <- evectors.karate %*% diag(filter_poly(evalues.karate, h=rev(tmp.heat$p))) %*% t(evectors.karate)
H.Mexican <- evectors.karate %*% diag(filter_poly(evalues.karate, h=rev(tmp.Mexican$p))) %*% t(evectors.karate)
H.sin.exp <- evectors.karate %*% diag(filter_poly(evalues.karate, h=rev(tmp.sin.exp$p))) %*% t(evectors.karate)
H.highpass <- evectors.karate %*% diag(filter_poly(evalues.karate, h=rev(tmp.highpass$p))) %*% t(evectors.karate)

realizations_x <- H.heat %*% white.inputs
realizations_y <- H.Mexican %*% white.inputs

gcpsd.estimate.karate1 <- cpsd.graph(x=realizations_x[,1:R.list[1]], y=realizations_y[,1:R.list[1]], S=L.karate)
gcpsd.estimate.karate2 <- cpsd.graph(x=realizations_x[,1:R.list[2]], y=realizations_y[,1:R.list[2]], S=L.karate)
gcpsd.estimate.karate3 <- cpsd.graph(x=realizations_x[,1:R.list[3]], y=realizations_y[,1:R.list[3]], S=L.karate)

plot(evalues.karate/lmax.karate, 
     true.cpsd, type="l", ylim=c(-0.01, 0.35), lwd=5, cex.lab=3, cex.axis=1.8,
     xlab="Normalized graph frequency", ylab="Cross-spectal density")
lines(evalues.karate/lmax.karate, gcpsd.estimate.karate1, type="l", col="red", lwd=5, lty=2)
lines(evalues.karate/lmax.karate, gcpsd.estimate.karate2, type="l", col="blue", lwd=5, lty=2)
lines(evalues.karate/lmax.karate, gcpsd.estimate.karate3, type="l", col="green", lwd=5, lty=2)
legend(legend=c("True", paste("R =", R.list)), lty=c(1, rep(2, length(R.list))), lwd=5, 
       col=c("black","red", "blue", "green"), "topright", cex=2, text.font=2)

# window graph Fourier-based method
tmp.heat.minnesota <- polyApprox(function(l) filter_kernel_heat(l, lmax=lmax.minnesota, tau=10), 0, lmax.minnesota, 5)
tmp.Mexican.minnesota <- polyApprox(function(l) filter_kernel_Mexican(l, lmax=lmax.minnesota,  tau1=5, tau2=25), 0, lmax.minnesota, 7)
tmp.sin.exp.minnesota <- polyApprox(function(l) filter_kernel(l, lmax=lmax.minnesota, tau1=5, tau2=5), 0, lmax.minnesota, 9)
tmp.highpass.minnesota <- polyApprox(function(l) filter_kernel_highpass(l, lmax=lmax.minnesota, tau=1), 0, lmax.minnesota, 4)

H.heat.minnesota <- evectors.minnesota %*% diag(filter_poly(evalues.minnesota, h=rev(tmp.heat.minnesota$p))) %*% t(evectors.minnesota)
H.Mexican.minnesota <- evectors.minnesota %*% diag(filter_poly(evalues.minnesota, h=rev(tmp.Mexican.minnesota$p))) %*% t(evectors.minnesota)
H.sin.exp.minnesota <- evectors.minnesota %*% diag(filter_poly(evalues.minnesota, h=rev(tmp.sin.exp.minnesota$p))) %*% t(evectors.minnesota)
H.highpass.minnesota <- evectors.minnesota %*% diag(filter_poly(evalues.minnesota, h=rev(tmp.highpass.minnesota$p))) %*% t(evectors.minnesota)

# 1 realization case for windowed estimator
M.list <- c(10,100,1000)
set.seed(1)
white.inputs.window.minnesota <- t(mvrnorm(n=1, mu=rep(0,N.minnesota), Sigma=diag(1,N.minnesota))) # 1 x N.minnesota matrix
realizations_x.window.minnesota <- H.sin.exp.minnesota %*% t(white.inputs.window.minnesota)
realizations_y.window.minnesota <- H.highpass.minnesota %*% t(white.inputs.window.minnesota)

gcpsd.estimate.minnesota1.window <- cpsd.graph(x=realizations_x.window.minnesota, y=realizations_y.window.minnesota, S=L.minnesota, g=g, M=M.list[1], method="window")
gcpsd.estimate.minnesota2.window <- cpsd.graph(x=realizations_x.window.minnesota, y=realizations_y.window.minnesota, S=L.minnesota, g=g, M=M.list[2], method="window")
gcpsd.estimate.minnesota3.window <- cpsd.graph(x=realizations_x.window.minnesota, y=realizations_y.window.minnesota, S=L.minnesota, g=g, M=M.list[3], method="window")

true.cpsd.minnesota <- filter_poly(evalues.minnesota, h=rev(tmp.sin.exp.minnesota$p))*Conj(filter_poly(evalues.minnesota, h=rev(tmp.highpass.minnesota$p)))
plot(evalues.minnesota/lmax.minnesota, 
     true.cpsd.minnesota, type="l",
     ylim=c(-0.08, 0.1), lwd=5, cex.lab=3, cex.axis=1.8,
     xlab="Normalized graph frequency", ylab="Cross-spectal density")
lines(c(0:M.list[1])*(M.list[1]+1)/(M.list[1]^2), 
      c(approxExtrap(c(1:M.list[1])*(M.list[1]+1)/(M.list[1]^2),gcpsd.estimate.minnesota1.window, xout=0)$y,
        gcpsd.estimate.minnesota1.window), type="l", col="red", lwd=5, lty=2)
lines(c(0:M.list[2])*(M.list[2]+1)/(M.list[2]^2), 
      c(approxExtrap(c(1:M.list[2])*(M.list[2]+1)/(M.list[2]^2),gcpsd.estimate.minnesota2.window, xout=0)$y,
        gcpsd.estimate.minnesota2.window), type="l", col="green", lwd=5, lty=2)
lines(c(0:M.list[3])*(M.list[3]+1)/(M.list[3]^2), 
      c(approxExtrap(c(1:M.list[3])*(M.list[3]+1)/(M.list[3]^2),gcpsd.estimate.minnesota3.window, xout=0)$y,
        gcpsd.estimate.minnesota3.window), type="l", col="magenta", lwd=5, lty=2)
legend(legend=c("True", paste("K =", c(10,100,1000))), lty=c(1, rep(2, 3)), lwd=5, 
       col=c("black", "red", "green", "magenta"), "topright", cex=2, text.font=2)


#########################################
### 03. detection of common frequency ###
#########################################

## windowed average method ##
x <- 5*evectors.metro[,100] + 100*evectors.metro[,50]
y <- 5*evectors.metro[,100] + 100*evectors.metro[,150]

WB <- windowbank.random(N=N.metro, M=100, V=evectors.metro, sigma=0.5)
cpsd.datawindow.metro <- cpsd.graph(x=t(WB)*x, y=t(WB)*y, S=L.metro, g=g, M=M.metro)
psd.datawindow.metro.x <- cpsd.graph(x=t(WB)*x, y=t(WB)*x, S=L.metro, g=g, M=M.metro)
psd.datawindow.metro.y <- cpsd.graph(x=t(WB)*y, y=t(WB)*y, S=L.metro, g=g, M=M.metro)

x <- 5*evectors.irregular01[,300] + 100*evectors.irregular01[,100]
y <- 5*evectors.irregular01[,300] + 100*evectors.irregular01[,200]

WB <- windowbank.random(N=N.irregular01, M=100, V=evectors.irregular01, sigma=0.5)
cpsd.datawindow.irregular01 <- cpsd.graph(x=t(WB)*x, y=t(WB)*y, S=L.irregular01, g=g, M=M.irregular01)
psd.datawindow.irregular01.x <- cpsd.graph(x=t(WB)*x, y=t(WB)*x, S=L.irregular01, g=g, M=M.irregular01)
psd.datawindow.irregular01.y <- cpsd.graph(x=t(WB)*y, y=t(WB)*y, S=L.irregular01, g=g, M=M.irregular01)

x <- 5*evectors.minnesota[,500] + 100*evectors.minnesota[,1000]
y <- 5*evectors.minnesota[,500] + 100*evectors.minnesota[,2000]

WB <- windowbank.random(N=N.minnesota, M=100, V=evectors.minnesota, sigma=0.5)
cpsd.datawindow.minnesota <- cpsd.graph(x=t(WB)*x, y=t(WB)*y, S=L.minnesota, g=g, M=M.minnesota)
psd.datawindow.minnesota.x <- cpsd.graph(x=t(WB)*x, y=t(WB)*x, S=L.minnesota, g=g, M=M.minnesota)
psd.datawindow.minnesota.y <- cpsd.graph(x=t(WB)*y, y=t(WB)*y, S=L.minnesota, g=g, M=M.minnesota)

# visualize
par(mfrow=c(3,3), oma=c(0,1,0,0))
plot(evalues.metro, psd.datawindow.metro.x, main="GPSD (x)", pch=16, lwd=2,
     xlab="Graph frequency",ylab="Power spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.metro[50], lty=2, col="red", lwd=3)
plot(evalues.metro, psd.datawindow.metro.y, main="GPSD (y)", pch=16, lwd=2,
     xlab="Graph frequency",ylab="Power spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.metro[150], lty=2, col="red", lwd=3)
plot(evalues.metro, cpsd.datawindow.metro, main="GCSD", pch=16, lwd=2,
     xlab="Graph frequency", ylab="Cross-spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.metro[100], lty=2, col="red", lwd=3)

plot(evalues.irregular01, psd.datawindow.irregular01.x, main="GPSD (x)", pch=16, lwd=2,
     xlab="Graph frequency",ylab="Power spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.irregular01[100], lty=2, col="red", lwd=3)
plot(evalues.irregular01, psd.datawindow.irregular01.y, main="GPSD (y)", pch=16, lwd=2,
     xlab="Graph frequency",ylab="Power spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.irregular01[200], lty=2, col="red", lwd=3)
plot(evalues.irregular01, cpsd.datawindow.irregular01, main="GCSD", pch=16, lwd=2,
     xlab="Graph frequency", ylab="Cross-spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.irregular01[300], lty=2, col="red", lwd=3)

plot(evalues.minnesota, psd.datawindow.minnesota.x, main="GPSD (x)", pch=16, lwd=2,
     xlab="Graph frequency",ylab="Power spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.minnesota[1000], lty=2, col="red", lwd=3)
plot(evalues.minnesota, psd.datawindow.minnesota.y, main="GPSD (y)", pch=16, lwd=2,
     xlab="Graph frequency",ylab="Power spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.minnesota[2000], lty=2, col="red", lwd=3)
plot(evalues.minnesota, cpsd.datawindow.minnesota, main="GCSD", pch=16, lwd=2,
     xlab="Graph frequency", ylab="Cross-spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.minnesota[500], lty=2, col="red", lwd=3)


## window graph Fourier-based method ##
# metro
x <- 5*evectors.metro[,100] + 100*evectors.metro[,50]
y <- 5*evectors.metro[,100] + 100*evectors.metro[,150]
M.metro <- 100
cpsd.window.metro <- cpsd.graph(x=x, y=y, S=L.metro, g=g, M=M.metro, method="window")
psd.window.metro.x <- cpsd.graph(x=x, y=x, S=L.metro, g=g, M=M.metro, method="window")
psd.window.metro.y <- cpsd.graph(x=y, y=y, S=L.metro, g=g, M=M.metro, method="window")
tau.metro <- (M.metro+1)*lmax.metro / M.metro^2

# visualize
g0 <- plot_graph_custom4(seoulmetro, e.size=1.3, v.size=3, vertex_color = x, 
                         value="value", signal=FALSE)
g1 <- plot_graph_custom3(seoulmetro, e.size=1.3, v.size=7, vertex_color = x,
                         min=-52, max=52, value="value")
g2 <- plot_graph_custom3(seoulmetro, e.size=1.3, v.size=7, vertex_color = y,
                         min=-52, max=52, value="value")
g3 <- plot_graph_custom3(seoulmetro, e.size=1.3, v.size=7, vertex_color = evectors.metro[,100],
                         min=-0.52, max=0.52, value="value")
g4 <- plot_graph_custom3(seoulmetro, e.size=1.3, v.size=7, vertex_color = evectors.metro[,50],
                         min=-0.52, max=0.52, value="value")
g5 <- plot_graph_custom3(seoulmetro, e.size=1.3, v.size=7, vertex_color = evectors.metro[,150],
                         min=-0.52, max=0.52, value="value")


grid.arrange(g0,g1,g2,g3,g4,g5, nrow=2)

# irregular
x <- 5*evectors.irregular01[,300] + 100*evectors.irregular01[,100]
y <- 5*evectors.irregular01[,300] + 100*evectors.irregular01[,200]
M.irregular01 <- 50
cpsd.window.irregular01 <- cpsd.graph(x=x, y=y, S=L.irregular01, g=g, M=M.irregular01, method="window")
psd.window.irregular01.x <- cpsd.graph(x=x, y=x, S=L.irregular01, g=g, M=M.irregular01, method="window")
psd.window.irregular01.y <- cpsd.graph(x=y, y=y, S=L.irregular01, g=g, M=M.irregular01, method="window")
tau.irregular01<- (M.irregular01+1)*lmax.irregular01 / M.irregular01^2

# visualize
g0 <- plot_graph_custom4(irregular01, e.size=1.3, v.size=3, vertex_color = x, 
                         value="value", signal=FALSE)
g1 <- plot_graph_custom3(irregular01, e.size=1.3, v.size=7, vertex_color = x,
                         min=-57, max=57, value="value")
g2 <- plot_graph_custom3(irregular01, e.size=1.3, v.size=7, vertex_color = y,
                         min=-57, max=57, value="value")
g3 <- plot_graph_custom3(irregular01, e.size=1.3, v.size=7, vertex_color = evectors.irregular01[,300],
                         min=-0.73, max=0.73, value="value")
g4 <- plot_graph_custom3(irregular01, e.size=1.3, v.size=7, vertex_color = evectors.irregular01[,100],
                         min=-0.73, max=0.73, value="value")
g5 <- plot_graph_custom3(irregular01, e.size=1.3, v.size=7, vertex_color = evectors.irregular01[,200],
                         min=-0.73, max=0.73, value="value")


grid.arrange(g0,g1,g2,g3,g4,g5, nrow=2)

# Minnesota
x <- 5*evectors.minnesota[,500] + 100*evectors.minnesota[,1000]
y <- 5*evectors.minnesota[,500] + 100*evectors.minnesota[,2000]
M.minnesota <- 100
cpsd.window.minnesota <- cpsd.graph(x=x, y=y, S=L.minnesota, g=g, M=M.minnesota, method="window")
psd.window.minnesota.x <- cpsd.graph(x=x, y=x, S=L.minnesota, g=g, M=M.minnesota, method="window")
psd.window.minnesota.y <- cpsd.graph(x=y, y=y, S=L.minnesota, g=g, M=M.minnesota, method="window")
tau.minnesota<- (M.minnesota+1)*lmax.minnesota / M.minnesota^2

# visualize
g0 <- plot_graph_custom4(Minnesota.network2, e.size=1.3, v.size=3, vertex_color = x, 
                         value="value", signal=FALSE)
g1 <- plot_graph_custom3(Minnesota.network2, e.size=1.3, v.size=5, vertex_color = x[order(x)],
                         min=-38, max=38, value="value")
g2 <- plot_graph_custom3(Minnesota.network2, e.size=1.3, v.size=5, vertex_color = y[order(x)],
                         min=-38, max=38, value="value")
g3 <- plot_graph_custom3(Minnesota.network2, e.size=1.3, v.size=5, vertex_color = evectors.minnesota[order(x),500],
                         min=-0.44, max=0.44, value="value")
g4 <- plot_graph_custom3(Minnesota.network2, e.size=1.3, v.size=5, vertex_color = evectors.minnesota[order(x),1000],
                         min=-0.44, max=0.44, value="value")
g5 <- plot_graph_custom3(Minnesota.network2, e.size=1.3, v.size=5, vertex_color = evectors.minnesota[order(x),2000],
                         min=-0.44, max=0.44, value="value")


grid.arrange(g0,g1,g2,g3,g4,g5, nrow=2)



# estimation results
par(mfrow=c(3,3), oma=c(0,1,0,0))
plot(tau.metro*c(1:M.metro), psd.window.metro.x, main="GPSD (x)", pch=16, lwd=2,
     xlab="Graph frequency",ylab="Power spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.metro[50], lty=2, col="red", lwd=3)
plot(tau.metro*c(1:M.metro), psd.window.metro.y, main="GPSD (y)", pch=16, lwd=2,
     xlab="Graph frequency",ylab="Power spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.metro[150], lty=2, col="red", lwd=3)
plot(tau.metro*c(1:M.metro), cpsd.window.metro, main="GCSD", pch=16, lwd=2,
     xlab="Graph frequency", ylab="Cross-spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.metro[100], lty=2, col="red", lwd=3)

plot(tau.irregular01*c(1:M.irregular01), psd.window.irregular01.x, main="GPSD (x)", pch=16, lwd=2,
     xlab="Graph frequency",ylab="Power spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.irregular01[100], lty=2, col="red", lwd=3)
plot(tau.irregular01*c(1:M.irregular01), psd.window.irregular01.y, main="GPSD (y)", pch=16, lwd=2,
     xlab="Graph frequency",ylab="Power spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.irregular01[200], lty=2, col="red", lwd=3)
plot(tau.irregular01*c(1:M.irregular01), cpsd.window.irregular01, main="GCSD", pch=16, lwd=2,
     xlab="Graph frequency", ylab="Cross-spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.irregular01[300], lty=2, col="red", lwd=3)

plot(tau.minnesota*c(1:M.minnesota), psd.window.minnesota.x, main="GPSD (x)", pch=16, lwd=2,
     xlab="Graph frequency",ylab="Power spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.minnesota[1000], lty=2, col="red", lwd=3)
plot(tau.minnesota*c(1:M.minnesota), psd.window.minnesota.y, main="GPSD (y)", pch=16, lwd=2,
     xlab="Graph frequency",ylab="Power spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.minnesota[2000], lty=2, col="red", lwd=3)
plot(tau.minnesota*c(1:M.minnesota), cpsd.window.minnesota, main="GCSD", pch=16, lwd=2,
     xlab="Graph frequency", ylab="Cross-spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
abline(v=evalues.minnesota[500], lty=2, col="red", lwd=3)



#############################
### 04. robust estimators ###
#############################

# robust estimator for the GPSD
original.x <- 3*evectors.karate[,20]
outliers.x <- original.x
outliers.x[25] <- 4
plot(outliers.x)


WB <- windowbank.random(N=N.karate, M=100, V=evectors.karate, sigma=0.5)
psd.datawindow.karate.x.outlier <- cpsd.graph(x=t(WB)*outliers.x, y=t(WB)*outliers.x, S=L.karate, g=g, M=M.karate)
psd.datawindow.karate.x.original <- cpsd.graph(x=t(WB)*original.x, y=t(WB)*original.x, S=L.karate, g=g, M=M.karate)

# robust version
Sigma.window.karate <- (t(WB)*outliers.x) %*% t(t(WB)*outliers.x) / 100
sigma.window.karate <- as.vector(Sigma.window.karate)
G <- matrix(0, nrow=N.karate^2, ncol=N.karate)
for(i in 1:N.karate){
  G[,i] <- as.vector(as.matrix(evectors.karate[,i]) %*% t(as.matrix(evectors.karate[,i])))
}

# 50% quantile-based
psd.datawindow.karate.x.robust <- PIRLS(lambda=0, y=sigma.window.karate, tau=0.5, alpha=0.25, B=G,
                                        P=L.karate, gamma_init=rep(mean(outliers.x), length(outliers.x)), max_iterations=10000, tol=1e-6, check="Oh")


# robust estimator fot the GCSD
original.x <- 5*evectors.karate[,20] + 20*evectors.karate[,30]
original.y <- 5*evectors.karate[,20] + 20*evectors.karate[,10]
outliers.x <- original.x
outliers.x[25] <- -10
outliers.y <- original.y
outliers.y[15] <- 10

cpsd.datawindow.karate.original <- cpsd.graph(x=t(WB)*original.x, y=t(WB)*original.y, S=L.karate)
cpsd.datawindow.karate.outlier <- cpsd.graph(x=t(WB)*outliers.x, y=t(WB)*outliers.y, S=L.karate)

# robust version
Sigma.window.karate.cross <- (t(WB)*outliers.x) %*% t(t(WB)*outliers.y) / 100
sigma.window.karate.cross <- as.vector(Sigma.window.karate.cross)

cpsd.datawindow.karate.robust <- PIRLS(lambda=0, y=sigma.window.karate.cross, tau=0.5, alpha=0.25, B=G,
                                       P=L.karate, gamma_init=rep(mean(outliers.x), length(outliers.x)), max_iterations=10000, tol=1e-6, check="Oh")


# estimation results
par(mfrow=c(1,2), oma=c(0,1,0,0))
plot(evalues.karate, psd.datawindow.karate.x.original, main="GPSD", pch=16, lwd=3,
     xlab="Graph frequency",ylab="Power spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
lines(evalues.karate, psd.datawindow.karate.x.outlier, col="blue", lwd=3, lty=2)
lines(evalues.karate, psd.datawindow.karate.x.robust, col="green", lwd=3, lty=2)
abline(v=evalues.karate[20], lty=3, col="red", lwd=3)

plot(evalues.karate, cpsd.datawindow.karate.original, main="GCSD", pch=16, lwd=3, ylim=c(-2,60),
     xlab="Graph frequency",ylab="Cross-spectral density", type="l", cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
lines(evalues.karate, cpsd.datawindow.karate.outlier, col="blue", lwd=3, lty=2)
lines(evalues.karate, cpsd.datawindow.karate.robust, col="green", lwd=3, lty=2)
abline(v=evalues.karate[20], lty=3, col="red", lwd=3)



