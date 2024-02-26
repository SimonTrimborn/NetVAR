.LambdaSelect = function(LossFunction, test, Ydata_eval, Model) {
  test_eval = test[[1]]
  
  if (Model == "LISAR.Adap.LASSO") {
    lambdas_store = matrix(NA, 
                           nrow = length(test[[2]][[1]]) * length(test[[2]][[2]]) * length(test[[2]][[3]]) * length(test[[2]][[4]]), 
                           ncol = length(test[[2]]) - 1)
  } else {
    lambdas_store = matrix(NA, 
                           nrow = length(test[[2]][[1]]) * length(test[[2]][[2]]) * length(test[[2]][[3]]) * length(test[[2]][[4]]), 
                           ncol = length(test[[2]]))
  }
  
  TT.eval = dim(Ydata_eval[[1]])[2]
  counter = 0
  eval_save = c()
  
  for (j1 in 1:length(test[[2]][[1]])) {
    for (j2 in 1:length(test[[2]][[2]])) {
      for (j3 in 1:length(test[[2]][[3]])) {
        counter = counter + 1
        
        res1 = Map('%*%', test_eval[[j1]][[j2]][[j3]], Ydata_eval[-1])
        res2 = Reduce('+',res1)
        if (LossFunction == "BIC") {
          eval_save[counter] = TT.eval*log(det(((Ydata_eval[[1]] - res2) %*% t((Ydata_eval[[1]] - res2))) / TT.eval)) + sum(unlist(test_eval[[j1]][[j2]][[j3]]) != 0) * log(TT.eval)
          lambdas_store[counter,] = c(test[[2]][[1]][j1], test[[2]][[2]][j2], test[[2]][[3]][j3], test[[2]][[4]])
        } else if (LossFunction == "AIC") {
          eval_save[counter] = TT.eval*log(det(((Ydata_eval[[1]] - res2) %*% t((Ydata_eval[[1]] - res2))) / TT.eval)) + sum(unlist(test_eval[[j1]][[j2]][[j3]]) != 0) * 2
          lambdas_store[counter,] = c(test[[2]][[1]][j1], test[[2]][[2]][j2], test[[2]][[3]][j3], test[[2]][[4]])
        } else if (LossFunction == "MSFE") {
          eval_save[counter] = mean((Ydata_eval[[1]] - res2)^2)
          lambdas_store[counter,] = c(test[[2]][[1]][j1], test[[2]][[2]][j2], test[[2]][[3]][j3], test[[2]][[4]])
        }
      }
    }
  }
  
  lambdas_select = lambdas_store[head(which(eval_save == min(eval_save)), 1),]
  
  # lambda1
  lambda1_upper = max(which(test[[2]][[1]] == lambdas_select[1]) - 1, 1)
  lambda1_lower = min(which(test[[2]][[1]] == lambdas_select[1]) + 1, length(test[[2]][[1]]))
  
  # lambda2
  lambda2_upper = max(which(test[[2]][[2]] == lambdas_select[2]) - 1, 1)
  lambda2_lower = min(which(test[[2]][[2]] == lambdas_select[2]) + 1, length(test[[2]][[2]]))
  
  # lambda3
  lambda3_upper = max(which(test[[2]][[3]] == lambdas_select[3]) - 1, 1)
  lambda3_lower = min(which(test[[2]][[3]] == lambdas_select[3]) + 1, length(test[[2]][[3]]))
  
  # optimal model
  # sometimes the lambdas are no longer unique due to a too small sequence. Then [1] becomes relevant
  j1 = which(test[[2]][[1]] == lambdas_select[1])[1]
  j2 = which(test[[2]][[2]] == lambdas_select[2])[1]
  j3 = which(test[[2]][[3]] == lambdas_select[3])[1]
  opt_model = test_eval[[j1]][[j2]][[j3]]
  
  lambda_return = list(c(lambda1_upper, lambda1_lower), 
                       c(lambda2_upper, lambda2_lower), 
                       c(lambda3_upper, lambda3_lower))
  lambda_select = list(lambdas_select[1],lambdas_select[2],lambdas_select[3])
  
  return(list(lambda_return, opt_model, lambda_select))

}