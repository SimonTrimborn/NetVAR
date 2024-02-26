.LASSO2 = function(N, Lags, Ydata, lambda3, eps1) {
  
  Yr = Ydata[[1]]
  Yx = Ydata[-1]
  store3 = list()
      for(j3 in 1:length(lambda3)) {
        store3[[j3]] = rep(list(matrix(0, N, N)), Lags) 
      }
  
    stepi = 0
        for(j3 in 1:length(lambda3)) {

          b = store3[[j3]]

          b1 = do.call(cbind, b)
          b1 = .FistaLassoK(B = b1, Y = Yr, Z = do.call(rbind, Yx), gam = lambda3[j3], eps = eps1)
          b = b1[[1]]
          b = lapply(seq_len(Lags), function(x) b[,(1:N)+N*(x-1)])
          store3[[j3]] = b
        }
  lambdas = list(lambda3)
  return(list(store3, lambdas))
}
