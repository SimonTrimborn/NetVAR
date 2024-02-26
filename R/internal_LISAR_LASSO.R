.LISAR_LASSO = function(N, TT, Lags, Ydata, lambda1, lambda2, lambda3, a.pen, eps1, alpha.pen, reshape) {

  pb = txtProgressBar(min = 0, max = length(lambda1) * length(lambda2)  * length(lambda3), initial = 0, style=3)
  
  Yr = Ydata[[1]]
  Yx = Ydata[-1]
  store3 = list()
  for(j1 in 1:length(lambda1)) {
    store2 = list()
    for(j2 in 1:length(lambda2)) {
      store1 = list()
      for(j3 in 1:length(lambda3)) {
        store1[[j3]] = rep(list(matrix(0, N, N)), Lags)
      }
      store2[[j2]] = store1
    }
    store3[[j1]] = store2
  }
  
  stepi = 0
  
  skip_all = rep(FALSE, length(lambda3))
  
    for(j1 in 1:length(lambda1)) {
      store2 = list()
      for(j2 in 1:length(lambda2)) {
        store1 = list()
        for(j3 in 1:length(lambda3)) {
          stepi = stepi + 1
          setTxtProgressBar(pb,stepi)
          
          if (j1 > 1 & skip_all[j3] == TRUE & reshape == TRUE) {
            store3[[j1]][[j2]][[j3]] = store3[[j1-1]][[j2]][[j3]]
            next
          } else if (skip_all[j3] == TRUE & reshape == FALSE) {
            store3[[j1]][[j2]][[j3]] = store3[[j1-1]][[j2]][[j3]]
            next
          }
          
          if (j1 > 1) {
            ## derive check values for lambda
            store_testvalue = .RemainInfo(b=store3[[j1-1]][[j2]][[j3]], Yx=Yx, Yr=Yr, TT=TT)
            eval_testvals1 = store_testvalue <= lambda1[j1-1]
            eval_testvals2 = store_testvalue <= lambda1[j1]
            
            ## setting values to old ones when lambda structure similar
            if (all(eval_testvals1 == eval_testvals2) & any(eval_testvals2 == FALSE)) {
              store3[[j1]][[j2]][[j3]] = store3[[j1-1]][[j2]][[j3]]
              next
            } else if (any(eval_testvals1 != eval_testvals2) & any(eval_testvals2 == FALSE)) {
              if (j2 > 1) {
                eval_testvals_group1 = c()
                eval_testvals_group2 = c()
                for (i in 1:Lags) {
                  YR1 = Map('%*%', store3[[j1]][[j2-1]][[j3]][-i], Yx[-i])
                  YR2 = Reduce('+',YR1)
                  if (is.null(YR2)) {Yr.diff = Yr-0} else {Yr.diff = Yr-YR2}
                  
                  a = store3[[j1]][[j2-1]][[j3]][[i]]
                  diag(a) = 0
                  YR = a %*% Yx[[i]] - Yr.diff 
                  M = matrix(0, ncol = dim(Yx[[i]])[1], nrow = dim(Yx[[i]])[1] * dim(Yx[[i]])[2])
                  for ( k in 1:dim(Yx[[i]])[2]) {M[(1 + dim(Yr.diff)[1]*(k-1)):(dim(Yr.diff)[1]*k),]=diag(Yx[[i]][,k])} 
                  testval = t(M) %*% matrix(YR,ncol=1)
                  
                  store_testval = sqrt(sum(testval^2))
                  for (j in 1:N) {
                    a = store3[[j1]][[j2-1]][[j3]][[i]]
                    a[-j,j] = 0
                    YR = (a %*% Yx[[i]] - Yr.diff)[-j,] 
                    M = Yx[[i]][j,] 
                    testval = M %*% t(YR)
                    store_testval = c(store_testval, sqrt(sum(testval^2)))
                  }
                  eval_testvals_group1 = c(eval_testvals_group1, store_testval <= lambda2[j2-1])
                  eval_testvals_group2 = c(eval_testvals_group2, store_testval <= lambda2[j2])
                }
                if (all(eval_testvals_group2 == eval_testvals_group1)) {
                  store3[[j1]][[j2]][[j3]] = store3[[j1]][[j2-1]][[j3]]
                  next
                }
              }
            } else if (all(eval_testvals1 == eval_testvals2) & all(eval_testvals2)) {
              store3[[j1]][[j2]][[j3]] = rep(list(matrix(0, N, N)), Lags)
              next
            }
          }
          
          b = store3[[j1]][[j2]][[j3]]
          if (j3 > 1 & reshape == TRUE) {
            # take previous group settings optimized result
            b = store3[[j1]][[j2]][[j3-1]]
          } else if (j1 > 1 & j2 > 1 & j3 > 1 & reshape == FALSE) { 
            # take previous group settings optimized result
            b = store3[[j1]][[j2]][[j3-1]]
          }
          
          store_testval = .RemainInfo(b=b, Yx=Yx, Yr=Yr, TT=TT)
          
          ## LAGS
          eval_testvals = store_testval <= lambda1[j1]
          b_order = order(store_testval, decreasing = TRUE)[0:sum(!eval_testvals)]
          store.k1 = rep(0, Lags)
          b_order_new = b_order
          b_order_new_store = c()
          iset_new_store = c()
          iset_new_store_list = list()
          
          k2 = 0
          while(k2 < Lags) {
            k2 = k2 + 1
            if (k2 > length(b_order_new)) {next}
            i = b_order_new[k2]
            b_order_new_store = c(b_order_new_store, i)
            
            YR1 = Map('%*%', b[-i], Yx[-i])
            YR2 = Reduce('+',YR1)
            if (is.null(YR2)) {Yr.diff = Yr-0} else {Yr.diff = Yr-YR2}
            store_testval = c()
            for (j in 1:N) {
              a = b[[i]]
              a[-j,j] = 0
              YR = (a %*% Yx[[i]] - Yr.diff)[-j,] 
              M = Yx[[i]][j,]
              testval = M %*% t(YR)
              store_testval = c(store_testval, sqrt(sum(testval^2)))
            }
            eval_testvals_groups = store_testval <= lambda2[j2]
            iset = order(store_testval, decreasing = TRUE)[0:sum(!eval_testvals_groups)]
            ##
            # Create a matrix of regularization parameters with alpha and 1/alpha
            # times lambda. The diagonal has currently no more penalty as it is 
            # usually nonzero when time series effects exist. Stronger penalty, 
            # 1/alpha * lambda for uninformative groups, weaker penalty, 
            # alpha * lambda for informative groups. 
            iset_new = iset
            iset_current = c()
            
            k1 = 0
            while(k1 < N) {
              k1 = k1 + 1
              if (k1 > length(iset_new)) {next}
              j = iset_new[k1]
              iset_current = c(iset_current, j)
              reg_para_matrix = matrix(1/alpha.pen * lambda3[j3], N, N * k2)
              if (k2 > 1) {
                which_elements = unlist(Map(function(x, y) relist(unlist(x) * y, skeleton = x), 1:length(iset_new_store_list), iset_new_store_list))
                reg_para_matrix[,sort(which_elements)] = alpha.pen * lambda3[j3]
              }
              reg_para_matrix[,sort(iset_current) + (N * (k2-1))] = alpha.pen * lambda3[j3]
              for (jjj in 1:k2) {
                diag(reg_para_matrix[,(1:N) + (N*(jjj-1))]) = alpha.pen * lambda3[j3]
              }
              
              b1 = .FistaLASSOK.mixed.pen(B = do.call(cbind, b[b_order_new[1:k2]]), 
                                        Y = Yr, 
                                        Z = do.call(rbind, Yx[b_order_new[1:k2]]), 
                                        gam = reg_para_matrix, eps = eps1)
              for (jjj in 1:k2) {
                b[[b_order_new[jjj]]] = b1[[1]][,(1:N) + (N * (jjj-1))]
              }

              store_testval = c()
              for (jj in 1:N) {
                a = b[[i]]
                a[-jj,jj] = 0
                YR = (a %*% Yx[[i]] - Yr.diff)[-jj,] 
                M = Yx[[i]][jj,]
                testval = M %*% t(YR)
                store_testval = c(store_testval, sqrt(sum(testval^2)))
              }
              # prevent current active groups from being filtered out of the new set of active groups
              store_testval[iset_current] = 0
              eval_testvals_groups = store_testval <= lambda2[j2]
              iset_new = unique(c(iset_current, order(store_testval, decreasing = TRUE)[0:sum(!eval_testvals_groups)]))
            }
            
            store_testval_new = .RemainInfo(b=b, Yx=Yx, Yr=Yr, TT=TT)
            # prevent current active groups from being filtered out of the new set of active groups
            store_testval_new[b_order_new_store] = 0
            eval_testvals_new = store_testval_new <= lambda1[j1]
            b_order_new = unique(c(b_order_new_store, order(store_testval_new, decreasing = TRUE)[0:sum(!eval_testvals_new)]))
            iset_new_store = c(iset_new_store, iset_new)
            iset_new_store_list[[k2]] = iset_new
          }
          
          store3[[j1]][[j2]][[j3]] = b
          
          # Check function. If all Lags and columns are active for this combination 
          # of lambda1 and lambda2 (j1, j2), then the solution is only determined
          # by lambda3 (j3). Hence we can skip all derivations from now on.
          if (length(b_order_new) == Lags & length(iset_new_store) == Lags * N) {
            skip_all[j3] = TRUE
          }
        }
      }
    }

  close(pb)
  lambdas = list(lambda1, lambda2, lambda3, alpha.pen)
  return(list(store3, lambdas))
}  
 

