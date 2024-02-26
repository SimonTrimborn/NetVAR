.RemainInfo = function(b, Yx, Yr, TT) {
  store_testval = c()
  for (i in 1:length(b)) {
    res1 = Map('%*%', b[-i], Yx[-i])
    res2 = Reduce('+',res1)
    if (is.null(res2)) {YR = Yr - 0} else {YR = Yr - res2}
    M = Yx[[i]]
    testval = M %*% t(YR)
    testval2 = sqrt(sum(testval^2))
    testval3 = testval2
    store_testval = c(store_testval, testval3)
  }
  return(store_testval)
}

.SCAD = function(x, lambda, a.p) {
  Sfu = (sign(x)*pmax(abs(x) - lambda,0)) * (abs(x) <= 2*lambda) + (((a.p - 1) * x - sign(x) * a.p * lambda)/ (a.p - 2)) * (abs(x) > 2*lambda & abs(x) <= a.p*lambda) + x * (abs(x) > a.p*lambda)
  return(Sfu)
}

.LASSO = function(x, lambda) {
  Sfu = (sign(x)*pmax(abs(x) - lambda,0))
  return(Sfu)
}

.LambdaSequence = function(Ydata, Lags, N, sequence.l.1, sequence.l.2, sequence.l.3) {
  Yr = Ydata[[1]]
  Yx = Ydata[-1]
  dat1 = list()
  for (i in 1:Lags) {
    dat1[[i]] = Yx[[i]] %*% t(Yr)
  }
  
  lambdamax1 = c()
  for (i in 1:Lags) {
    lambdamax1[i] = sqrt(sum(dat1[[i]]^2))
  }
  
  lambdamax2 = c()
  lambdamax21 = c()
  for (i in 1:Lags) {
    lambdamax21 = c()
    dat2 = dat1[[i]]
    diag(dat2) = 0
    for (j in 1:N) {
      lambdamax21[j] = sqrt(sum(dat2[j,]^2)) 
    }
    lambdamax2[i] = max(lambdamax21)
  }
  
  lambda1 = max(lambdamax1)
  lambda2 = max(lambdamax2)
  lambda3 = max(abs(unlist(dat1)))
  
  k=1
  while (lambda1[k] > 0.0000001) {k=k+1; lambda1[k] = lambda1[k-1]*sequence.l.1} 
  k=1
  while (lambda2[k] > 0.0000001) {k=k+1; lambda2[k] = lambda2[k-1]*sequence.l.2}
  k=1
  while (lambda3[k] > 0.0000001) {k=k+1; lambda3[k] = lambda3[k-1]*sequence.l.3}
  lambda1 = c(lambda1, 0)
  lambda2 = c(lambda2, 0) 
  lambda3 = c(lambda3, 0)
  
  return(list(lambda1, lambda2, lambda3))
}

.FistaLassoK <- function(B, Y, Z, gam, eps) {
  OVF <- c()
  k <- nrow(Y)
  BOLD <- B
  BOLDOLD <- B
  tk <- 1/max(Mod(eigen(Z %*% t(Z), only.values = T)$values))
  q = 1
  for (i in 1:k) {
    threshold <- 10 * eps
    j = 1
    while (threshold > eps) {
      v <- matrix(BOLD[i,] + (j - 2)/(j + 1) * (BOLD[i,] - BOLDOLD[i,]), nrow = 1)
      B[i,] <- .LASSO(v + tk * (Y[i,] - v %*% Z) %*% t(Z), gam * tk) 
      threshold <- max(abs(B[i,] - v)/(1 + v))
      OVF[q] <- max(abs(B[i,] - v)/(1 + v))
      BOLDOLD[i,] <- BOLD[i,]
      BOLD[i,] <- B[i,]
      j = j + 1 
      q = q + 1 
    }
  }
  return(list(B = B, OV = OVF, iters = q))
}

.FistaSCADK <- function(B, Y, Z, gam, eps,a.p) {
  OVF <- c()
  k <- nrow(Y)
  BOLD <- B
  BOLDOLD <- B
  tk <- 1/max(Mod(eigen(Z %*% t(Z), only.values = T)$values))
  q = 1
  for (i in 1:k) {
    threshold <- 10 * eps
    j = 1
    while (threshold > eps) {
      v <- matrix(BOLD[i,] + (j - 2)/(j + 1) * (BOLD[i,] - BOLDOLD[i,]), nrow = 1)
      B[i,] <- .SCAD(v + tk * (Y[i,] - v %*% Z) %*% t(Z), gam * tk,a.p) 
      threshold <- max(abs(B[i,] - v)/(1 + v))
      OVF[q] <- max(abs(B[i,] - v)/(1 + v))
      BOLDOLD[i,] <- BOLD[i,]
      BOLD[i,] <- B[i,]
      j = j + 1 
      q = q + 1 
    }
  }
  return(list(B = B, OV = OVF, iters = q))
}

.FistaSCADK.mixed.pen <- function(B, Y, Z, gam, eps,a.p) {
  OVF <- c()
  k <- nrow(Y)
  BOLD <- B
  BOLDOLD <- B
  tk <- 1/max(Mod(eigen(Z %*% t(Z), only.values = T)$values))
  q = 1
  for (i in 1:k) {
    threshold <- 10 * eps
    j = 1
    while (threshold > eps) {
      v <- matrix(BOLD[i,] + (j - 2)/(j + 1) * (BOLD[i,] - BOLDOLD[i,]), nrow = 1)
      B[i,] <- .SCAD(v + tk * (Y[i,] - v %*% Z) %*% t(Z), gam[i,] * tk,a.p) 
      threshold <- max(abs(B[i,] - v)/(1 + v))
      OVF[q] <- max(abs(B[i,] - v)/(1 + v))
      BOLDOLD[i,] <- BOLD[i,]
      BOLD[i,] <- B[i,]
      j = j + 1 
      q = q + 1 
    }
  }
  return(list(B = B, OV = OVF, iters = q))
}

.FistaLASSOK.mixed.pen <- function(B, Y, Z, gam, eps) {
  OVF <- c()
  k <- nrow(Y)
  BOLD <- B
  BOLDOLD <- B
  tk <- 1/max(Mod(eigen(Z %*% t(Z), only.values = T)$values))
  q = 1
  for (i in 1:k) {
    threshold <- 10 * eps
    j = 1
    while (threshold > eps) {
      v <- matrix(BOLD[i,] + (j - 2)/(j + 1) * (BOLD[i,] - BOLDOLD[i,]), nrow = 1)
      B[i,] <- .LASSO(v + tk * (Y[i,] - v %*% Z) %*% t(Z), gam[i,] * tk) 
      threshold <- max(abs(B[i,] - v)/(1 + v))
      OVF[q] <- max(abs(B[i,] - v)/(1 + v))
      BOLDOLD[i,] <- BOLD[i,]
      BOLD[i,] <- B[i,]
      j = j + 1 
      q = q + 1 
    }
  }
  return(list(B = B, OV = OVF, iters = q))
}
