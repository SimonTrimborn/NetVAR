#' The LISAR model
#' 
#' \code{\link{LISAR}} contains the LISAR model described in Zhang and 
#' Trimborn (2023). The model is computed based on 3 lambda tuning parameters 
#' which are automatically determined by the package. A scaling parameter 
#' boosts/lowers the strength of the regularization on assets which values are 
#' determined less/more important than others (influencers). The function 
#' allows to compute the LISAR model under \code{"LASSO"}, \code{"SCAD"} and 
#' \code{"AdapLASSO"}. 
#' 
#' @param x A matrix containing the data with rows being observations and
#' columns being time series.
#' @param Model The LISAR model to use with either LASSO 
#' \code{"LISAR.LASSO"}, SCAD \code{"LISAR.SCAD"} or Adaptive LASSO 
#' \code{"LISAR.AdapLASSO"} regularization. Default \code{"LISAR.LASSO"}. 
#' @param eval.criteria The evaluation criteria to use to choose the 
#' best model under the regularization parameters. Can be either \code{"MSFE"}, 
#' \code{"AIC"} or \code{"BIC"}. Default \code{"MSFE"}. 
#' @param Lags The maximum number of lags to consider. Default 3. 
#' @param alpha.pens A number or vector specifying the boosting parameter 
#' increasing/decreasing the strength of regularization. Should be a number(s) 
#' between (0,1). See Zhang and Trimborn (2023) for details on the alpha 
#' parameters. 
#' @param gamma.pens A number or vector specifying the adaptive parameter 
#' for \code{"LISAR.AdapLASSO"}. Only required for \code{"LISAR.AdapLASSO"}. The number(s) 
#' should be larger than 0. 
#' @param lambda1_seq The factor by which the regularization sequence of 
#' lambda1, regularizing the lag structure, decreases towards 0. Should be a 
#' value between (0,1). 
#' @param lambda2_seq The factor by which the regularization sequence of 
#' lambda2, indicating which time series are more influential (influencers), 
#' decreases towards 0. Should be a value between (0,1). 
#' @param lambda3_seq The factor by which the regularization sequence of 
#' lambda3, regularizing the individual parameters, decreases towards 0. 
#' Should be a value between (0,1). 
#' @param a.pen The parameter specifying by which the SCAD penalty taperes off 
#' towards no regularization. Only required for \code{"LISAR.SCAD"}. 
#' The number should be larger than 0. Default \code{3.7}.
#' @param eps1 Control parameter for the inner optimization algorithm. The algorithm 
#' converged when between optimization steps the parameters change by less than 
#' \code{"eps1"}. Default \code{0.0001}.
#' @param eps2 Control parameter for the outer optimization algorithm. The algorithm 
#' converged when between optimization steps the parameters change by less than 
#' \code{"eps2"}. Default \code{0.0001}.
#' @param T1 A numeric stating the row of \code{"x"} where the training data 
#' end and the evaluation period starts. If \code{NULL}, then the first third 
#' of data are chosen as training data. Defaults to \code{NULL}.
#' @param T2 A numeric stating the row of \code{"x"} where the evaluation data 
#' end and the out-of-sample period starts. If \code{NULL}, then the second third 
#' of data are chosen as evaluation data. Defaults to \code{NULL}.
#' @param reoptim Logical. If TRUE, then the best model found under the 
#' initially derived lambda sequences is further optimized by a new lambda 
#' sequence around the previous best solution. Stops when a more granular 
#' lambda sequence no longer improves the model under \code{"eval.criteria"} 
#' criterion. Default \code{"FALSE"}. 
#' @return An object of the \code{class} netstruc with the components
#' \item{Model.optimal}{a list containing the optimal model, evaluation 
#' criteria and model regularizers}
#' \item{data.training}{the training data} \item{data.evaluation}{the evaluation 
#' data} \item{data.outofsample}{the out-of-sample data.} 
#' \item{Model.estimation}{estimation specifics for the model}
#' @references Kexin Zhang, Simon Trimborn (2023).
#' Influential assets in Large-Scale Vector Autoregressive Models
#' \emph{SSRN Working paper}.
#' \doi{http://dx.doi.org/10.2139/ssrn.4619531}\cr \cr 
#' @examples
#' 
#' \dontrun{
#' data(TradingData)
#' 
#' LISAR(TradingData, Model = "LISAR.LASSO")
#' }
#' @export LISAR
LISAR = function(x, Model = "LISAR.LASSO", eval.criteria = "MSFE", Lags = 3, 
                 alpha.pens = 0.5, gamma.pens = 0.5, lambda1_seq = 0.5, 
                 lambda2_seq = 0.5, lambda3_seq = 0.5, a.pen = 3.7, 
                 eps1 = 0.0001, eps2 = 0.0001, T1 = NULL, T2 = NULL, 
                 reoptim = FALSE ) {
  
  if (!is.matrix(x)) {
    stop("x has to be a matrix.")
  }
  if (!is.element(Model, c("LISAR.LASSO", "LISAR.SCAD", "LISAR.AdapLASSO"))) {
    stop("Model has to be either LISAR.LASSO, LISAR.SCAD, or LISAR.AdapLASSO")
  }
  if (!is.element(eval.criteria, c("MSFE", "AIC", "BIC"))) {
    stop("eval.criteria has to be either MSFE, AIC, or BIC")
  }
  if (Lags %% 1 != 0 || Lags < 0) {
    stop("Lags has to be a positive integer.")
  }
  if (any(alpha.pens <=0) || any(alpha.pens >= 1)) {
    stop("alpha.pens has to be a number or vector with every entry between 0 and 1.")
  }
  if (any(lambda1_seq <=0) || any(lambda1_seq >= 1)) {
    stop("lambda1_seq has to be a number between 0 and 1.")
  }
  if (any(lambda2_seq <=0) || any(lambda2_seq >= 1)) {
    stop("lambda2_seq has to be a number between 0 and 1.")
  }
  if (any(lambda3_seq <=0) || any(lambda3_seq >= 1)) {
    stop("lambda3_seq has to be a number between 0 and 1.")
  }
  if (length(lambda1_seq) != 1 || length(lambda2_seq) != 1 || length(lambda3_seq) != 1) {
    stop("lambda1_seq, lambda2_seq, and lambda3_seq have to be a number, not vector.")
  }
  if (length(a.pen) != 1 || a.pen < 0 ) {
    stop("a.pen has to be a number larger or equal 0.")
  }
  if (length(eps1) != 1 || eps1 < 0 ) {
    stop("eps1 has to be a number larger 0.")
  }
  if (length(eps2) != 1 || eps2 < 0 ) {
    stop("eps2 has to be a number larger 0.")
  }
  if (!is.null(T1) & is.null(T2) || !is.null(T2) & is.null(T1)) {
    stop("T1 and T2 have to be either both NULL or both whole numbers.")
  }
  if (!is.null(T1) & !is.null(T2)) {
    if (T1 < 2 || T1 > T2 || T1 > dim(x)[2] ) {
      stop("T1 has to be larger 1, smaller than T2, and cannot exceed the length of the time series.")
    }
    if (T2 < T1 || T2 > dim(x)[2] ) {
      stop("T2 has to be larger than T1 and cannot exceed the length of the time series.")
    }
  }
  if (!is.logical(reoptim)) {
    stop("reoptim has to be TRUE or FALSE.")
  }
  

  
  dat = scale(x)
  
  TT = nrow(dat)
  N = ncol(dat)
  
  tmp = t(as.matrix(dat))
  
  if (is.null(T1)) {
    T1 = dim(tmp)[2]/3
  }
  if (is.null(T2)) {
    T2 = 2*dim(tmp)[2]/3
  }
  
  # training data
  Ydata = list()
  for (i in 0:Lags) {
    Ydata[[i+1]] = tmp[, (Lags+1-i):(T1-i)] 
  }
  
  # evaluation data
  Ydata_eval = list()
  for (i in 0:Lags) {
    Ydata_eval[[i+1]] = tmp[, (T1+(1-i)):(T2-i)]
  }
  
  # out of sample data
  Ydata_oos = list()
  for (i in 0:Lags) {
    Ydata_oos[[i+1]] = tmp[, (T2+(1-i)):(dim(tmp)[2]-i)]
  }
  
  lambdas = .LambdaSequence(Ydata = Ydata, Lags = Lags, N = N, sequence.l.1 = lambda1_seq, sequence.l.2 = lambda2_seq, sequence.l.3 = lambda3_seq)
  
  lambda1 = lambdas[[1]] 
  lambda2 = lambdas[[2]] 
  lambda3 = lambdas[[3]] 
  
  ##########################
  ### LISAR SCAD
  ##########################
  if (Model == "LISAR.SCAD") {
    store.time1 = Sys.time()
    store_model_alpha = list()
    
    for (k in 1:length(alpha.pens)) {
      alpha_opt = alpha.pens[k]
      
      message(paste0("START estimation ", Model, " with alpha = ", alpha_opt))
      
      store_model_alpha[[k]] = .LISAR_SCAD(N=N, TT=TT, Lags = Lags, Ydata = Ydata, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, a.pen = a.pen, eps1 = eps1, alpha.pen=alpha_opt, reshape = FALSE)
      
      message(paste0("DONE estimation ", Model, " with alpha = ", alpha_opt))
      
      lambda1_new = lambda1
      lambda2_new = lambda2
      lambda3_new = lambda3
      lambdas_and_opt_model = .LambdaSelect(eval.criteria, store_model_alpha[[k]], Ydata_eval, Model)
      lambdas_upper_lower = lambdas_and_opt_model[[1]]
      new_opt_model = lambdas_and_opt_model[[2]]
      last_opt_model = lapply(new_opt_model, function(x) x + 1)
      lambda1_new = seq(lambda1_new[lambdas_upper_lower[[1]][1]], lambda1_new[lambdas_upper_lower[[1]][2]], length.out = 10)
      lambda2_new = seq(lambda2_new[lambdas_upper_lower[[2]][1]], lambda2_new[lambdas_upper_lower[[2]][2]], length.out = 10)
      lambda3_new = seq(lambda3_new[lambdas_upper_lower[[3]][1]], lambda3_new[lambdas_upper_lower[[3]][2]], length.out = 10)
      
      # Refine optimal lambda combination
      if (reoptim == TRUE) {
        message(paste0("START refined estimation ", Model, " with alpha = ", alpha_opt))
        while (any((unlist(last_opt_model) - unlist(new_opt_model)) > eps2)) {
          store_model_alpha[[k]] = .LISAR_SCAD(N=N, TT=TT, Lags = Lags, Ydata = Ydata, 
                                              lambda1 = lambda1_new, 
                                              lambda2 = lambda2_new, 
                                              lambda3 = lambda3_new, a.pen = a.pen, eps1 = eps1, 
                                              alpha.pen=alpha_opt, 
                                              reshape = TRUE)
          
          lambdas_and_opt_model = .LambdaSelect(eval.criteria, store_model_alpha[[k]], Ydata_eval, Model)
          lambdas_upper_lower = lambdas_and_opt_model[[1]]
          
          lambda1_select = lambdas_and_opt_model[[3]][[1]]
          lambda2_select = lambdas_and_opt_model[[3]][[2]]
          lambda3_select = lambdas_and_opt_model[[3]][[3]]
          
          lambda1_new = seq(lambda1_new[lambdas_upper_lower[[1]][1]], lambda1_new[lambdas_upper_lower[[1]][2]], length.out = 10)
          lambda2_new = seq(lambda2_new[lambdas_upper_lower[[2]][1]], lambda2_new[lambdas_upper_lower[[2]][2]], length.out = 10)
          lambda3_new = seq(lambda3_new[lambdas_upper_lower[[3]][1]], lambda3_new[lambdas_upper_lower[[3]][2]], length.out = 10)
          
          last_opt_model = new_opt_model
          new_opt_model = lambdas_and_opt_model[[2]]
        }
        message(paste0("DONE refined estimation ", Model, " with alpha = ", alpha_opt))
      }
    }
    
    alpha_new = .AlphaSelect(eval.criteria, store_model_alpha,Ydata_eval,alpha.pens)
    which_model = which(alpha.pens == alpha_new[[1]])
    
    store_model_alpha_select = store_model_alpha[[which_model]]
    
    store.time2 = Sys.time()
    store.time2 - store.time1
  }
  
  
  ##########################
  ### LISAR LASSO
  ##########################
  if (Model == "LISAR.LASSO") {
    store.time1 = Sys.time()
    store_model_alpha = list()
    
    for (k in 1:length(alpha.pens)) {
      alpha_opt = alpha.pens[k]
      
      message(paste0("START estimation ", Model, " with alpha = ", alpha_opt))
      
      store_model_alpha[[k]] = .LISAR_LASSO(N=N, TT=TT, Lags = Lags, Ydata = Ydata, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, a.pen = a.pen, eps1 = eps1, alpha.pen=alpha_opt, reshape = FALSE)
      
      message(paste0("DONE estimation ", Model, " with alpha = ", alpha_opt))
      
      lambda1_new = lambda1
      lambda2_new = lambda2
      lambda3_new = lambda3
      lambdas_and_opt_model = .LambdaSelect(eval.criteria, store_model_alpha[[k]], Ydata_eval, Model)
      lambdas_upper_lower = lambdas_and_opt_model[[1]]
      new_opt_model = lambdas_and_opt_model[[2]]
      last_opt_model = lapply(new_opt_model, function(x) x + 1)
      lambda1_new = seq(lambda1_new[lambdas_upper_lower[[1]][1]], lambda1_new[lambdas_upper_lower[[1]][2]], length.out = 10)
      lambda2_new = seq(lambda2_new[lambdas_upper_lower[[2]][1]], lambda2_new[lambdas_upper_lower[[2]][2]], length.out = 10)
      lambda3_new = seq(lambda3_new[lambdas_upper_lower[[3]][1]], lambda3_new[lambdas_upper_lower[[3]][2]], length.out = 10)
      
      # Refine optimal lambda combination
      if (reoptim == TRUE) {
        message(paste0("START refined estimation ", Model, " with alpha = ", alpha_opt))
        while (any((unlist(last_opt_model) - unlist(new_opt_model)) > eps2)) {
          store_model_alpha[[k]] = .LISAR_LASSO(N=N, TT=TT, Lags = Lags, Ydata = Ydata, 
                                               lambda1 = lambda1_new, 
                                               lambda2 = lambda2_new, 
                                               lambda3 = lambda3_new, a.pen = a.pen, eps1 = eps1, 
                                               alpha.pen=alpha_opt, 
                                               reshape = TRUE)
          
          lambdas_and_opt_model = .LambdaSelect(eval.criteria, store_model_alpha[[k]], Ydata_eval, Model)
          lambdas_upper_lower = lambdas_and_opt_model[[1]]
          
          lambda1_select = lambdas_and_opt_model[[3]][[1]]
          lambda2_select = lambdas_and_opt_model[[3]][[2]]
          lambda3_select = lambdas_and_opt_model[[3]][[3]]
          
          lambda1_new = seq(lambda1_new[lambdas_upper_lower[[1]][1]], lambda1_new[lambdas_upper_lower[[1]][2]], length.out = 10)
          lambda2_new = seq(lambda2_new[lambdas_upper_lower[[2]][1]], lambda2_new[lambdas_upper_lower[[2]][2]], length.out = 10)
          lambda3_new = seq(lambda3_new[lambdas_upper_lower[[3]][1]], lambda3_new[lambdas_upper_lower[[3]][2]], length.out = 10)
          last_opt_model = new_opt_model
          new_opt_model = lambdas_and_opt_model[[2]]
        }
        message(paste0("DONE refined estimation ", Model, " with alpha = ", alpha_opt))
      }
    }
    
    alpha_new = .AlphaSelect(eval.criteria, store_model_alpha,Ydata_eval,alpha.pens)
    which_model = which(alpha.pens == alpha_new[[1]])
    
    store_model_alpha_select = store_model_alpha[[which_model]]
    
    store.time2 = Sys.time()
    store.time2 - store.time1
  }
  
  
  
  ##########################
  ### LISAR Adaptive LASSO
  ##########################
  if (Model == "LISAR.Adap.LASSO") {
    store.time1 = Sys.time()
    store_model_alpha_select = list()
    
    for (adap.loop in 1:length(gamma.pens)) {
      gamma.pen = gamma.pens[adap.loop]
      store_model_alpha = list()
      for (k in 1:length(alpha.pens)) {
        alpha_opt = alpha.pens[k]
        
        message(paste0("START estimation ", Model, " with alpha = ", alpha_opt, " and gamma = ", gamma.pen))
        
        store_model_alpha[[k]] = .LISAR_AdapLASSO(N=N, TT=TT, Lags = Lags, Ydata = Ydata, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, gamma.pen = gamma.pen, eps1 = eps1, alpha.pen=alpha_opt, reshape = FALSE)
        
        message(paste0("DONE estimation ", Model, " with alpha = ", alpha_opt, " and gamma = ", gamma.pen))
        
        lambda1_new = lambda1
        lambda2_new = lambda2
        lambda3_new = lambda3
        lambdas_and_opt_model = .LambdaSelect(eval.criteria, store_model_alpha[[k]], Ydata_eval, Model)
        lambdas_upper_lower = lambdas_and_opt_model[[1]]
        new_opt_model = lambdas_and_opt_model[[2]]
        last_opt_model = lapply(new_opt_model, function(x) x + 1)
        lambda1_new = seq(lambda1_new[lambdas_upper_lower[[1]][1]], lambda1_new[lambdas_upper_lower[[1]][2]], length.out = 10)
        lambda2_new = seq(lambda2_new[lambdas_upper_lower[[2]][1]], lambda2_new[lambdas_upper_lower[[2]][2]], length.out = 10)
        lambda3_new = seq(lambda3_new[lambdas_upper_lower[[3]][1]], lambda3_new[lambdas_upper_lower[[3]][2]], length.out = 10)
        
        # Refine optimal lambda combination
        if (reoptim == TRUE) {
          message(paste0("START refined estimation ", Model, " with alpha = ", alpha_opt, " and gamma = ", gamma.pen))
          while (any((unlist(last_opt_model) - unlist(new_opt_model)) > eps2)) {
            
            store_model_alpha[[k]] = .LISAR_AdapLASSO(N=N, TT=TT, Lags = Lags, Ydata = Ydata, 
                                                     lambda1 = lambda1_new, 
                                                     lambda2 = lambda2_new, 
                                                     lambda3 = lambda3_new, gamma.pen = gamma.pen, eps1 = eps1, 
                                                     alpha.pen=alpha_opt, 
                                                     reshape = TRUE)
            
            lambdas_and_opt_model = .LambdaSelect(eval.criteria, store_model_alpha[[k]], Ydata_eval, Model)
            lambdas_upper_lower = lambdas_and_opt_model[[1]]
            
            lambda1_select = lambdas_and_opt_model[[3]][[1]]
            lambda2_select = lambdas_and_opt_model[[3]][[2]]
            lambda3_select = lambdas_and_opt_model[[3]][[3]]
            
            lambda1_new = seq(lambda1_new[lambdas_upper_lower[[1]][1]], lambda1_new[lambdas_upper_lower[[1]][2]], length.out = 10)
            lambda2_new = seq(lambda2_new[lambdas_upper_lower[[2]][1]], lambda2_new[lambdas_upper_lower[[2]][2]], length.out = 10)
            lambda3_new = seq(lambda3_new[lambdas_upper_lower[[3]][1]], lambda3_new[lambdas_upper_lower[[3]][2]], length.out = 10)
            
            last_opt_model = new_opt_model
            new_opt_model = lambdas_and_opt_model[[2]]
          }
          message(paste0("DONE refined estimation ", Model, " with alpha = ", alpha_opt, " and gamma = ", gamma.pen))
        }
      }
      
      alpha_new = .AlphaSelect(eval.criteria, store_model_alpha, Ydata_eval, alpha.pens)
      which_model = which(alpha.pens == alpha_new[[1]])
      store_model_alpha_select[[adap.loop]] = store_model_alpha[[which_model]]
    }
    store.time2 = Sys.time()
    store.time2 - store.time1
  }
  
  EvaluateModel = .EvaluateLossFunction(eval.criteria, Model, store_model_alpha_select, alpha.pens, gamma.pens, 
                                        Ydata, 
                                        Ydata_oos, 
                                        store.time2, 
                                        store.time1)
  
   
   # Assigning NetVAR structure
   res <- structure(
     class = "NetVAR",
     list(Model.optimal = EvaluateModel, 
          data.training = Ydata, 
          data.evaluation = Ydata_eval, 
          data.outofsample = Ydata_oos,
          Model.estimation = c(Model, eval.criteria))
   )
   
   return(res)
  
}