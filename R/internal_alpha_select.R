.AlphaSelect = function(LossFunction,new_opt_model_alpha,Ydata_eval,alpha_test) {
  
  new_opt_model_min_candidates = c()
  new_opt_model_which_candidates = c()
  for (k in 1:length(new_opt_model_alpha)) {
    test_eval = new_opt_model_alpha[[k]][[1]]
    
    TT.eval = dim(Ydata_eval[[1]])[2]
    counter = 0
    eval_save = c()
    for (j1 in 1:length(new_opt_model_alpha[[k]][[2]][[1]])) {
      for (j2 in 1:length(new_opt_model_alpha[[k]][[2]][[2]])) {
        for (j3 in 1:length(new_opt_model_alpha[[k]][[2]][[3]])) {
          counter = counter + 1
          
          res1 = Map('%*%', test_eval[[j1]][[j2]][[j3]], Ydata_eval[-1])
          res2 = Reduce('+',res1)
          if (LossFunction == "BIC") {
            eval_save[counter] = TT.eval*log(det(((Ydata_eval[[1]] - res2) %*% t((Ydata_eval[[1]] - res2))) / TT.eval)) + sum(unlist(test_eval[[j1]][[j2]][[j3]]) != 0) * log(TT.eval)
          } else if (LossFunction == "AIC") {
            eval_save[counter] = TT.eval*log(det(((Ydata_eval[[1]] - res2) %*% t((Ydata_eval[[1]] - res2))) / TT.eval)) + sum(unlist(test_eval[[j1]][[j2]][[j3]]) != 0) * 2
          } else if (LossFunction == "MSFE") {
            eval_save[counter] = mean((Ydata_eval[[1]] - res2)^2)
          }
        }
      }
    }
    new_opt_model_min_candidates = c(new_opt_model_min_candidates, min(eval_save))
    new_opt_model_which_candidates = c(new_opt_model_which_candidates, head(which(eval_save == min(eval_save)), 1))
  }
  eval_candidates = new_opt_model_min_candidates
  eval_candidates_which = head(which(eval_candidates == min(eval_candidates)), 1)
  alpha_which = alpha_test[eval_candidates_which]
  
  eval_list = new_opt_model_alpha
  return(list(alpha_which, eval_list[[eval_candidates_which]]))
}
