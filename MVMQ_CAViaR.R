




####mqRQobjectiveFunction###
mqRQobjective <- function(beta0,y,THETA,empiricalQuantile,OUT){
  T <- dim(y)[1] ##raw
  N <- dim(y)[2] ##columns
  q <- matrix(0,T,N)
  q[1,] <- empiricalQuantile
  A = t(beta0)
  c = t(A[,1]); a= t(A[,2:(N+1)]);b=t(A[,(N+2):(1+2*N)])
  for(i in 2:T){
    q[i,] <- c+abs(y[t-1,])*a+q[t-1,]*b
  }
  S = (y-q) %*% (THETA*matrix(1,T,N)-ifelse(y<q,1,0))
  S = apply(S,1,sum)
  S = mean(S)
  if(OUT=1) return(S)
  else return(q)
}

############mvmqCAViaR###
mvmqCAViaR <- function(y,THETA){
  T <- dim(y)[1]
  N <- dim(y)[2]
  REP =50
  WIN = 100
  ysort  = matrix(0,WIN,N)
  empiricalQuantile <- matrix(0,1,N)
  c = matrix(NA,3,N)
  for(i in 1:N){ysort[,i] <- sort(y[1:WIN,i])}
  for(i in 1:N){empiricalQuantile[1,i] = ysort[round(WIN*THETA),i]}
  for(i in 1:N){c[,i] = CAViaR_estim(y[,i],THETA)}
  w = t(t(c[1,]))
  a = diag(c[2,])
  b = diag(c[3,])
  Beta0 <- as.vector(cbind(w,a,b))
  Beta <- fminsearch(mqRQobjective,Beta0,y,THETA,empiricalQuantile,1)
  for
}


standarderror <- function(y,THETA,Beta,q){
  T <- dim(y)[1]
  N <- dim(y)[2]
  A <- t(Beta)
  end <- dim(A)[2]
  dA1 <- cbind(diag(x=1,N,N),matrix(0,N,2*N^2))
  dA2 <- cbind(matrix(0,N^2,N),diag(1,N^2,N^2),matrix(0,N^2,N^2))
  dA3 <- cbind(matrix(0,N^2,N),matrix(0,N^2,N^2),diag(1,N^2,N^2))
  dq <- array(N,N+2*N^2,T)
  for(t in 2:T){
    dq[,,t]=dA1+kron(abs(y[t-1,]),diag(1,N,N))*dA2+A[,(1+N+1):end]*dq[,,t-1]+kron(q[t-1,],diag(1,N,N))*dA3
  }
  eps = y-q
  kk = median(abs(eps[,1]-median(eps[,1])))
  hh = T^(-1/3)*(qnorm(1-0.05/2)^(2/3))*((1.5*(dnorm(qnorm(THETA)))^2)/(2*(qnorm(THETA))^2+1))^(1/3)
  c = kk*(qnorm(THETA+hh)-qnorm(THETA-hh))
  print(c)
  print(sum(eps[,1]<c))
  Q <- matrix(0,N+2*N^2,N+2*N^2)
  V = Q
  for(i in 1:T){
    psi = THETA - as.matrix(ifelse(eps[t,]<0,1,0))
    eta <- as.matrix(dq[,,t]) * (psi %*% matrix(1,1,N+2*N^2))
    eta <- apply(eta,2,sum)
    V = V + t(eta) %*% eta
    Qt = matrix(0,N+2*N^2,N+2*N^2)
    for(j in 1:N){
      Qt = Qt+ifelse(abs(eps[t,j])<c,1,0)*as.matrix(dq[j,,t]) %*% dq[j,,t]
    }
    Q = Q+Qt
  }
  V = V/T
  Q =Q/(2*c*T)
  VC = solve(Q) %*% V %*% solve(Q)/T
  se = sqrt(diag(VC))
  R=matrix(c(0,0,0,1,0,0,0,0,0,0,
             0,0,0,0,1,0,0,0,0,0,
             0,0,0,0,0,0,0,1,0,0,
             0,0,0,0,0,0,0,0,1,0),4,10,byrow=TRUE)
  TS = t(R %*% Beta) %*% solve(R %*% VC %*% t(R))%*% (R %*% Beta)
  pValueTS = 1- pchisq(TS,4)
  return(list(se=se,TS=TS,pValueTS=pValueTS,VC=VC))
}

##################CAViaR_estim
CAViaR_estim(rt,0.01)

library(pracma)
CAViaR_estim <- function(y,THETA){
  REP <- 30
  nInitialCond <- 10
  ysort <- sort(y[1:300])
  empiricalQuantile <- ysort[round(300*THETA)]
  initialTargetVectors = matrix(runif(3000,0,1),1000,3)
  RQfval <- rep(0,1000)
  for(i in 1:1000){
    RQfval[i] <- RQobjectiveFunction(initialTargetVectors[i,],1,y,THETA,empiricalQuantile)
  }
  Results <- cbind(RQfval,initialTargetVectors)
  SortedResults = sortrows(Results,1)
  BestInitialCond = SortedResults[1:nInitialCond,2:4]
  Beta  = matrix(NA,nInitialCond,3)
  fval = rep(NA,nInitialCond)
  for(i in 1:nInitialCond){
    Beta[i,] <- fminsearch(RQobjectiveFunction,BestInitialCond[i,],1,y,THETA,empiricalQuantile)$xval
    for(j in 1:REP){
      res <- fminsearch(RQobjectiveFunction,Beta[i,],1,y,THETA,empiricalQuantile)
      Beta[i,] <- res$xval
      fval[i] <- res$fval
    }
  }
  SortedFval = sortrows(cbind(fval,Beta),1)
  BetaHat <- SortedFval[1,2:4]
  return(BetaHat)
}

##################RQobjectiveFunction######
RQobjectiveFunction <- function(BETA,OUT,y,THETA,empiricalQuantile){
  q <- SAVloop(BETA,y,empiricalQuantile)
  Hit = (y<q)-THETA
  RQ <- sum(-Hit*(y-q))
  if(OUT == 1){
    return(RQ)
    }else{
    return(cbind(q,Hit))
  }
}

SAVloop <- function(BETA,y,empiricalQuantile){
  T = length(y)
  q = rep(0,T)
  q[1] <- empiricalQuantile
  for(t in 2:T){
    q[t] <- BETA[1]+BETA[2]*abs(y[t-1])+BETA[3]*q[t-1]
  }
  return(q)
}

#######################

