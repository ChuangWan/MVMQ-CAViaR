##############VaR & ES test###
require('rugarch')


DQTest <- function(y,VaR,theta,lag){
  n <- length(y)
  Hit <- ifelse(y<VaR,1,0)-theta
  Hit.mat <- matrix(NA,(n-lag),lag)
  Y <- Hit[(lag+1):n]
  for(i in 1:lag){
    Hit.mat[,i] <- Hit[i:(n-lag+i-1)]
  }
  X <- cbind(1,Hit.mat,VaR[(lag+1):n])
  ols <- solve((t(X) %*% X)) %*% t(X) %*% Y
  DQ <- (t(ols)%*%t(X) %*% X %*% ols)/(theta*(1-theta))
  df <- lag+2
  p.DQ <- 1-pchisq(DQ,df)
  return(p.DQ)
}

VaRtest <- function(y,VaR,theta,lag){
  rures <- VaRTest(theta, as.numeric(actual), as.numeric(VaR))
  UCstat <- rures$uc.LRstat
  CCstat <- rures$cc.LRstat
  hit <- rures$actual.exceed
  failrate <- hit/length(y)
  p.UC <- rures$uc.LRp
  p.CC <- rures$cc.LRp
  indstat <- UCstat+CCstat
  p.Ind <- 1-pchisq(indstat,2)
  p.DQ <- DQTest(y, VaR,theta,lag)
  Pvalue <- list(p.UC=p.UC,p.CC=p.CC,p.Ind=p.Ind,p.DQ=p.DQ)
  res <- list(Pvalue=Pvalue,hit=hit,rate=failrate)
  return(res)
}
