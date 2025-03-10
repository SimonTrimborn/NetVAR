#' Printing function for objects of class NetVAR
#' 
#' \code{\link{summary.NetVAR}} prints the values of an object of class NetVAR.
#' 
#' @method summary NetVAR
#' 
#' @param object An object of \code{class} NetVAR
#' @param ...  Further arguments to be passed over. Currently none.
#' @return None
#' @export
summary.NetVAR <- function(object, ...) {
  x = object
    cat(rep("-", getOption("width")), sep = "")
    cat("\n")
    cat(paste0("Model estimated with method ", 
               x$Model.estimation[1], 
    " under ", x$Model.estimation[2], 
    " as evaluation criteria. The best model under this evaluation criteria has ", 
    x$Model.optimal$Model.evaluation[1], 
    " Lags and ", x$Model.optimal$Model.evaluation[2], 
    " time dependent parameters.", collapse = NULL))

    cat("\n")
    cat(rep("-", getOption("width")), sep = "")
    cat("\n")
    cat("Evaluation:")
    cat("\n")
    out <- c(
      paste0("\t MSE: ", round(x$Model.optimal$Model.evaluation[3], 4)), 
      paste0("\t MSE eval: ", round(x$Model.optimal$Model.evaluation[4], 4)), 
      paste0("\t MSFE: ", round(x$Model.optimal$Model.evaluation[5], 4)),
      paste0("\t Runtime (in seconds): ", round(x$Model.optimal$Model.evaluation[6], 4))
      )
    cat(paste(out, collapse = "\n"))
    cat("\n")
    cat(rep("-", getOption("width")), sep = "")
    cat("\n")
    cat("Regularizers:")
    cat("\n")
    out_names = names(x$Model.optimal$Model.regularizers)
    out <- character()
    for (i in seq_along(x$Model.optimal$Model.regularizers)) {
      out <- c(out, paste0("\t", out_names[i], ": ", round(x$Model.optimal$Model.regularizers[i], 4)))
    }
    cat(paste(out, collapse = "\n"))
    cat("\n")
    cat(rep("-", getOption("width")), sep = "")
    cat("\n")

      cat("Parameters:")
      for (i in seq_along(x$Model.optimal$Model)) {
        cat("\n")
        cat("\n")
        cat(paste0("Lag ", i, ":"))
        cat("\n")
        out <- format(round(x$Model.optimal$Model[[i]], 2))
        out = cbind(rownames(x$Model.optimal$Model[[i]]), out)
        out = rbind(c("",colnames(x$Model.optimal$Model[[i]])), out)
        out = format(out, justify = "right")
        out2 = character()
        for (j in 1:dim(out)[1]) {
          out2 = c(out2, paste0(out[j,], collapse = " "))
        }
        cat(paste(out2, collapse = "\n"))
        cat("\n")
      }
      
      cat(
        "To plot the parameter matrices as heatmaps, call plot(x)."
      )
    
    cat("\n")
    cat("\n")
    
  invisible(x)
}
