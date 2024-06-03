# localization operator
localization.op <- function(evalues, evectors, g, i=NULL){ # support of T_i g is centered at node i
  if(is.null(i)){
    res <- evectors %*% (t(Conj(evectors))*g(evalues)) # ith column: T_i g
  } else{
    res <- as.vector(evectors %*% (Conj(evectors)[i,]*g(evalues))) # T_i g vector when i is specified
  }
  return(res)
}

# graph windowed Fourier transform
graph.window.FT <- function(x, S, g, M){
  # eigenres <- eigen(S)
  # evalues <- eigenres$values
  # evectors <- eigenres$vectors
  # lmax <- max(evalues)
  val <- eigensort(S)
  evalues <- val$evalues
  evectors <- val$evectors
  lmax <- max(evalues)
  
  C <- NULL
  normconst <- c()
  for(m in 1:M){
    gm <- function(lambda){
      return(g(lambda, sigma.sq=lmax*(M+1)/M^2, m, tau=lmax*(M+1)/M^2))
    }
    Tgm <- localization.op(evalues, evectors, gm)
    C <- cbind(C, t(Conj(Tgm))%*%x) # where C[i,m] = <x, Ti gm>
    normconst <- c(normconst, norm(Tgm, type="F")^2)
  }
  return(list(C=C, normconst=normconst))
}

# cpsd estimate
cpsd.graph <- function(x, y, S, g=NULL, M=NULL, method=NULL){
  if(is.null(method)){
    # x, y should be R realizations (N x R matrix)
    # eigenres <- eigen(S)
    # x.tilde <- t(Conj(eigenres$vectors))%*%x
    # y.tilde <- t(Conj(eigenres$vectors))%*%y
    val <- eigensort(S)
    evalues <- val$evalues
    evectors <- val$evectors
    x.tilde <- t(Conj(evectors))%*%x
    y.tilde <- t(Conj(evectors))%*%y
    cpsd <- rowMeans(x.tilde * Conj(y.tilde))
  } else if(method=="window"){
    # x, y should be one realization (N x 1 vector)
    C1 <- graph.window.FT(x, S, g, M)
    C2 <- graph.window.FT(y, S, g, M)
    cpsd <- colSums(C1$C * Conj(C2$C)) / C1$normconst
  }
  return(cpsd)
}

# psd estimate
psd.graph <- function(x, S, g=NULL, M=NULL, method=NULL){
  return(cpsd.graph(x=x, y=x, S=S, g=g, M=M, method=method))
}

# coherence estimate
coherence.graph <- function(x, y, S, g=NULL, M=NULL, method=NULL){
  if(is.null(method)){
    cpsd <- cpsd.graph(x=x, y=y, S=S, g=g, M=M)
    psd.x <- cpsd.graph(x=x, y=x, S=S, g=g, M=M)
    psd.y <- cpsd.graph(x=y, y=y, S=S, g=g, M=M)
  } else if(method=="window"){
    cpsd <- cpsd.graph(x=x, y=y, S=S, g=g, M=M, method=method)
    psd.x <- cpsd.graph(x=x, y=x, S=S, g=g, M=M, method=method)
    psd.y <- cpsd.graph(x=y, y=y, S=S, g=g, M=M, method=method)
  }
  return(cpsd*Conj(cpsd) / psd.x / psd.y)
}

# cross spectrum analysis All in one
cross.spectrum.graph <- function(x, y, S, g=NULL, M=NULL, method=NULL){
  cpsd <- cpsd.graph.window(x=x, y=y, S=S, g=g, M=M, method)
  psd.x <- cpsd.graph.window(x=x, y=x, S=S, g=g, M=M, method)
  psd.y <- cpsd.graph.window(x=y, y=y, S=S, g=g, M=M, method)
  return(list(cpsd=cpsd, psd.x = psd.x, psd.y = psd.y, coherence=cpsd*Conj(cpsd) / psd.x / psd.y))
}

graph.stationary.level <- function(cov, S){
  val <- eigensort(S)
  evalues <- val$evalues
  evectors <- val$evectors
  gamma <- t(Conj(evectors))%*%cov%*%evectors
  return(norm(diag(gamma), type=2) / norm(gamma, type="F"))
}

# robust cpsd estimate 
huberloss <- function(t, c){
  res <- c()
  for(i in t){
    if(abs(i) <= c){
      res <- c(res, i^2)
    } else{
      res <- c(res, 2*c*(abs(i)-c/2))
    }
  }
  return(res)
}

PIRLS <- function(lambda, y, tau, alpha, B, P, gamma_init, max_iterations, tol, check){
  n <- length(y)
  gamma <- gamma_init
  for (iteration in 1:max_iterations){
    gamma.old <- gamma
    if(check=="Zheng"){
      W <- diag(grad_checkft_by_Zheng(as.numeric(y-B%*%gamma), tau, alpha) / as.numeric((2*(y-B%*%gamma))))
    } else if(check=="Oh"){
      W <- diag(grad_checkft_by_Oh(as.numeric(y-B%*%gamma), tau, alpha) / as.numeric((2*(y-B%*%gamma))))
    }
    
    tmp <- solve(t(B)%*%W%*%B + n*lambda*P)%*%t(B)%*%W
    gamma <- tmp%*%y
    gammahat <- gamma
    if (norm(gamma-gamma.old, type="2") < tol) {
      cat("Iteration:", iteration, "\n")
      cat("Converged!\n")
      break
    }
  }
  
  res <- gammahat
  
  return(res)
}

grad_checkft_by_Oh <- Vectorize(function(x, tau, c){
  if((-c<=x) & (x<c)){
    return(tau*x/c*as.numeric(x >= 0) + (1-tau)*x/c*as.numeric(x < 0))
  }
  else{
    return(tau - as.numeric(x < -c))
  }
}, vectorize.args = "x")

# random window generation
windowbank.random <- function(N, M, V, sigma){
  res <- matrix(0, nrow=M, ncol=N)
  for(i in 1:M){
    W.tilde <- diag(N) + rnorm(N^2, 0, sigma)
    res[i,] <- diag(V %*% W.tilde %*% t(Conj(V)))
  }
  return(res)
}


# cpsd estimate
# robust.cpsd.graph.window <- function(x, y, S, g, M){
#   C1 <- graph.window.FT(x, S, g, M)
#   C2 <- graph.window.FT(y, S, g, M)
#   cpsd <- colSums(C1$C * Conj(C2$C)) / C1$normconst
#   return(cpsd)
# }



plot_graph_custom <- function (z, size = 0.75, edge_color, vertex_color=NULL) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$color <- factor(edge_color, levels = unique(edge_color))
  p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                  xend = y1, yend = y2, colour = color), 
                                              size = 2, data = df2) +
    scale_color_manual(values=color.cand, labels = paste("Line", 1:8, sep=""),
                       name = "Line number") + 
    geom_point(aes(fill=vertex_color), size = size, shape=21) + 
    scale_fill_gradient(low="white", high="black", na.value = "yellow", name = "People") +
    theme_void() +
    theme(legend.margin = margin(10,10,10,10))
  print(p1)
}

plot_graph_custom3 <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  
  if(signal==FALSE){
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2),  linewidth=e.size, color = "gray", data = df2) + 
      geom_point(size = v.size) + theme_void()+theme(aspect.ratio=ratio)
    print(p1)
  }
  
  else{
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="grey", high="black", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) +
      theme(legend.margin = margin(10,10,10,10), plot.margin = margin(10,10,10,10), 
            plot.title=element_text(size=30, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
}

plot_graph_custom4 <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  
  if(signal==FALSE){
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="gray", high="gray", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = "black", limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) +
      theme(legend.margin = margin(10,10,10,10), plot.margin = margin(10,10,10,10), 
            plot.title=element_text(size=30, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
  
  else{
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="grey", high="black", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) +
      theme(legend.margin = margin(10,10,10,10), plot.margin = margin(10,10,10,10), 
            plot.title=element_text(size=30, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
}