ref_matrix <- function(Nbask){
  ref <- matrix(0, nrow = Nbask, ncol = Nbask)
  ref[upper.tri(ref)] <- 1:choose(Nbask, 2)
  ref <- ref + t(ref)
  
  ref
}

j_matrix <- function(Nbask){
  jmtx <- matrix(0, ncol = Nbask - 1, nrow = Nbask)
  for(i in 1:Nbask){
    jmtx[i,] <- (1:Nbask)[1:Nbask != i]
  }
  
  jmtx
}

pi_success <- function(fit, pn){
  ps <- colMeans(as.data.frame(fit > pn))
  names(ps) <- sapply(1:Nbask, function(x){paste(c("pi[", x, "]"), collapse = "")})
  ps
}

JS_distance<-function(I,n,x,q1,q0,a,b){
  D=matrix(NA,I,I)
  for (i in 1:(I-1))
  {
    for (j in (i+1):I)
    {
      if ((x[i]/n[i])==(x[j]/n[j]))
      {
        x[j]=x[j]
      }
    }
  }
  for (i in 1:I)
  {
    for (j in i:I)
    {
      if (i == j)
      {
        D[i,j] = 0
      }
      if (i != j)
      {
        cdf1 = pbeta(seq(0,1,1/100),a+x[i],b+n[i]-x[i])
        freqs1 = sapply(1:(length(cdf1)-1),FUN=function(x){cdf1[x+1]-cdf1[x]}) + 0.0001
        cdf2 = pbeta(seq(0,1,1/100),a+x[j],b+n[j]-x[j])
        freqs2 = sapply(1:(length(cdf2)-1),FUN=function(x){cdf2[x+1]-cdf2[x]}) + 0.0001
        D[i,j] = D[j,i] = (KL.plugin(freqs1, freqs2)+KL.plugin(freqs2, freqs1))/2
      }
    }
  }
  D
}
