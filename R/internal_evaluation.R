.EvaluateLossFunction = function(LossFunction, Model, model_res, alpha_vec, gamma_vec, Ydata, Ydata_oos, store.time2, store.time1) {
  
  if (Model == "LISAR.SCAD") {
    model_res_eval = model_res[[1]]
    lambdas_store = matrix(NA, 
                           nrow = length(model_res[[2]][[1]]) * length(model_res[[2]][[2]]) * length(model_res[[2]][[3]]) * length(alpha_vec), 
                           ncol = length(model_res[[2]]))
    
    TT.eval = dim(Ydata_oos[[1]])[2]
    counter = 0
    eval_save = c()
    
    for (j1 in 1:length(model_res[[2]][[1]])) {
      for (j2 in 1:length(model_res[[2]][[2]])) {
        for (j3 in 1:length(model_res[[2]][[3]])) {
          counter = counter + 1
          
          res1 = Map('%*%', model_res_eval[[j1]][[j2]][[j3]], Ydata_oos[-1])
          res2 = Reduce('+',res1)
          if (LossFunction == "BIC") {
            eval_save[counter] = TT.eval*log(det(((Ydata_oos[[1]] - res2) %*% t((Ydata_oos[[1]] - res2))) / TT.eval)) + sum(unlist(model_res_eval[[j1]][[j2]][[j3]]) != 0) * log(TT.eval)
            lambdas_store[counter,] = c(model_res[[2]][[1]][j1], model_res[[2]][[2]][j2], model_res[[2]][[3]][j3], model_res[[2]][[4]])
          } else if (LossFunction == "AIC") {
            eval_save[counter] = TT.eval*log(det(((Ydata_oos[[1]] - res2) %*% t((Ydata_oos[[1]] - res2))) / TT.eval)) + sum(unlist(model_res_eval[[j1]][[j2]][[j3]]) != 0) * 2
            lambdas_store[counter,] = c(model_res[[2]][[1]][j1], model_res[[2]][[2]][j2], model_res[[2]][[3]][j3], model_res[[2]][[4]])
          } else if (LossFunction == "MSFE") {
            eval_save[counter] = mean((Ydata_oos[[1]] - res2)^2)
            lambdas_store[counter,] = c(model_res[[2]][[1]][j1], model_res[[2]][[2]][j2], model_res[[2]][[3]][j3], model_res[[2]][[4]])
          }
        }
      }
    }
    lambdas_select = lambdas_store[tail(which(eval_save == min(eval_save)), 1),]
    names(lambdas_select) = c("lambda1", "lambda2", "lambda3", "alpha")
    j1 = which(model_res[[2]][[1]] == lambdas_select[1])
    j2 = which(model_res[[2]][[2]] == lambdas_select[2])
    j3 = which(model_res[[2]][[3]] == lambdas_select[3])
    opt_model = model_res_eval[[j1]][[j2]][[j3]]
  } else if (Model == "LISAR.LASSO") {
    model_res_eval = model_res[[1]]
    lambdas_store = matrix(NA, 
                           nrow = length(model_res[[2]][[1]]) * length(model_res[[2]][[2]]) * length(model_res[[2]][[3]]) * length(alpha_vec), 
                           ncol = length(model_res[[2]]))
    
    TT.eval = dim(Ydata_oos[[1]])[2]
    counter = 0
    eval_save = c()
    
    for (j1 in 1:length(model_res[[2]][[1]])) {
      for (j2 in 1:length(model_res[[2]][[2]])) {
        for (j3 in 1:length(model_res[[2]][[3]])) {
          counter = counter + 1
          
          res1 = Map('%*%', model_res_eval[[j1]][[j2]][[j3]], Ydata_oos[-1])
          res2 = Reduce('+',res1)
          if (LossFunction == "BIC") {
            eval_save[counter] = TT.eval*log(det(((Ydata_oos[[1]] - res2) %*% t((Ydata_oos[[1]] - res2))) / TT.eval)) + sum(unlist(model_res_eval[[j1]][[j2]][[j3]]) != 0) * log(TT.eval)
            lambdas_store[counter,] = c(model_res[[2]][[1]][j1], model_res[[2]][[2]][j2], model_res[[2]][[3]][j3], model_res[[2]][[4]])
          } else if (LossFunction == "AIC") {
            eval_save[counter] = TT.eval*log(det(((Ydata_oos[[1]] - res2) %*% t((Ydata_oos[[1]] - res2))) / TT.eval)) + sum(unlist(model_res_eval[[j1]][[j2]][[j3]]) != 0) * 2
            lambdas_store[counter,] = c(model_res[[2]][[1]][j1], model_res[[2]][[2]][j2], model_res[[2]][[3]][j3], model_res[[2]][[4]])
          } else if (LossFunction == "MSFE") {
            eval_save[counter] = mean((Ydata_oos[[1]] - res2)^2)
            lambdas_store[counter,] = c(model_res[[2]][[1]][j1], model_res[[2]][[2]][j2], model_res[[2]][[3]][j3], model_res[[2]][[4]])
          }
        }
      }
    }
    lambdas_select = lambdas_store[tail(which(eval_save == min(eval_save)), 1),]
    names(lambdas_select) = c("lambda1", "lambda2", "lambda3", "alpha")
    j1 = which(model_res[[2]][[1]] == lambdas_select[1])
    j2 = which(model_res[[2]][[2]] == lambdas_select[2])
    j3 = which(model_res[[2]][[3]] == lambdas_select[3])
    opt_model = model_res_eval[[j1]][[j2]][[j3]]
  } else if (Model == "LISAR.AdapLASSO") {
    model_res_eval = model_res
    lambdas_store = matrix(NA, 
                           nrow = length(model_res[[1]][[2]][[1]]) * length(model_res[[1]][[2]][[2]]) * length(model_res[[1]][[2]][[3]]) * length(alpha_vec) * length(gamma_vec), 
                           ncol = length(model_res[[1]][[2]]))
    
    TT.eval = dim(Ydata_oos[[1]])[2]
    counter = 0
    eval_save = c()
    
    for (j4 in 1:length(model_res)) {
      for (j1 in 1:length(model_res[[j4]][[2]][[1]])) {
        for (j2 in 1:length(model_res[[j4]][[2]][[2]])) {
          for (j3 in 1:length(model_res[[j4]][[2]][[3]])) {
            counter = counter + 1
            
            res1 = Map('%*%', model_res_eval[[j4]][[1]][[j1]][[j2]][[j3]], Ydata_oos[-1])
            res2 = Reduce('+',res1)
            if (LossFunction == "BIC") {
              eval_save[counter] = TT.eval*log(det(((Ydata_oos[[1]] - res2) %*% t((Ydata_oos[[1]] - res2))) / TT.eval)) + sum(unlist(model_res_eval[[j4]][[1]][[j1]][[j2]][[j3]]) != 0) * log(TT.eval)
              lambdas_store[counter,] = c(model_res[[j4]][[2]][[1]][j1], model_res[[j4]][[2]][[2]][j2], model_res[[j4]][[2]][[3]][j3], model_res[[j4]][[2]][[4]], model_res[[j4]][[2]][[5]])
            } else if (LossFunction == "AIC") {
              eval_save[counter] = TT.eval*log(det(((Ydata_oos[[1]] - res2) %*% t((Ydata_oos[[1]] - res2))) / TT.eval)) + sum(unlist(model_res_eval[[j4]][[1]][[j1]][[j2]][[j3]]) != 0) * 2
              lambdas_store[counter,] = c(model_res[[j4]][[2]][[1]][j1], model_res[[j4]][[2]][[2]][j2], model_res[[j4]][[2]][[3]][j3], model_res[[j4]][[2]][[4]], model_res[[j4]][[2]][[5]])
            } else if (LossFunction == "MSFE") {
              eval_save[counter] = mean((Ydata_oos[[1]] - res2)^2)
              lambdas_store[counter,] = c(model_res[[j4]][[2]][[1]][j1], model_res[[j4]][[2]][[2]][j2], model_res[[j4]][[2]][[3]][j3], model_res[[j4]][[2]][[4]], model_res[[j4]][[2]][[5]])
            }
          }
        }
      }
    }
    lambdas_select = lambdas_store[tail(which(eval_save == min(eval_save)), 1),]
    names(lambdas_select) = c("lambda1", "lambda2", "lambda3", "alpha", "gamma")
    
    gamma_model_res = c()
    for (j4 in 1:length(model_res)) {
      gamma_model_res = c(gamma_model_res, model_res[[j4]][[2]][[5]])
    }
    j4 = which(gamma_model_res == lambdas_select[5])
    j1 = which(model_res[[j4]][[2]][[1]] == lambdas_select[1])
    j2 = which(model_res[[j4]][[2]][[2]] == lambdas_select[2])
    j3 = which(model_res[[j4]][[2]][[3]] == lambdas_select[3])
    opt_model = model_res_eval[[j4]][[1]][[j1]][[j2]][[j3]]
  }

  LagsChosen = max(which(sapply(opt_model, function(x) !all(x == 0)) == TRUE),0)
  NumberIncludedPara = length(which(round(unlist(opt_model),4) != 0))
  MSE = mean((Ydata[[1]] - Reduce('+',Map('%*%', opt_model, Ydata[-1])))^2) 
  MSFE = mean((Ydata_oos[[1]] - Reduce('+',Map('%*%', opt_model, Ydata_oos[-1])))^2) 
  evals = c(LagsChosen, NumberIncludedPara, MSE, MSFE, as.numeric(store.time2) - as.numeric(store.time1))
  names(evals) = c("Chosen Lags", "Included parameters","MSE", "MSFE","estimation time")
  
  names(opt_model) = paste0("Lag", 1:length(opt_model))
  for (i in 1:length(opt_model)) {
    rownames(opt_model[[i]]) = rownames(Ydata[[1]])
    colnames(opt_model[[i]]) = paste0(rownames(Ydata[[1]]), ".", i)
  }
  
  return(list(Model = opt_model,  
              Model.evaluation = evals, 
              Model.regularizers = lambdas_select))
}