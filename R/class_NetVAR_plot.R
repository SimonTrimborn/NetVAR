#' Plotting function for objects of class NetVAR
#' 
#' Plots an object of class NetVAR
#' 
#' 
#' @method plot NetVAR
#' 
#' @param x An object of \code{class} NetVAR
#' @param join_graphs Logical. If TRUE, then heatmaps for each lagged 
#' parameter matrix are returned in one plot. If FALSE, a single figure per 
#' heatmap is generated.
#' @param scale_range Logical. If TRUE, parameters larger 1 or smaller -1 are 
#' truncated to 1 and -1 purely for visualisation purposes. If FALSE, then the 
#' actual values are plotted in the heatmaps. 
#' @param legend Add a legend to plots. Default FALSE 
#' @param ...  Further arguments, currently none.
#' @return None
#' @export
plot.NetVAR <- function(x, join_graphs = TRUE, scale_range = TRUE, legend = FALSE, ...) {
  x1 = x$Model.optimal$Model
  if (join_graphs == TRUE) {
    layout(matrix(1:length(x1), length(x1), 1, byrow = TRUE))
  } else {
    layout(matrix(1, 1, 1, byrow = TRUE))
  }
  for (i in seq_along(x1)) {
    D = x1[[i]]
    rotate = function(x) t(apply(x, 2, rev)) # rotate clockwise 90 degrees
    floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
    ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)
    if (scale_range == TRUE) {
      color = colorRampPalette(c('blue2','white','red'))(201)
      cols_seq = round(seq(-1,1,0.01),2)
      D[which(D > 1)] = 1
      D[which(D < -1)] = -1
    } else {
      range_values = max(D) - min(D)
      color = colorRampPalette(c('blue2','white','red'))(range_values)
      cols_seq = round(seq(min(D),max(D),0.01),2)
    }
    
    pos_min = which(cols_seq == floor_dec(min(D),2))
    pos_max = which(cols_seq == ceiling_dec(max(D),2))
    
    if (legend == FALSE) {
      image(rotate(D), col= color[pos_min:pos_max],bty="o", cex.lab = 0.5,
            las=1, axes=F)
    } else if (legend == TRUE) {
      image.plot(rotate(D), col= color[pos_min:pos_max],bty="o", cex.lab = 0.5,
            las=1, axes=F)
    }
    
    nlabs = dim(D)[2]
    tck1=seq(0,1, length=nlabs)
    tck2=seq(1,0, length=nlabs)
    tick.labels1=colnames(D)
    tick.labels2=rownames(D)
    
    axis(1, at=tck1, labels=tick.labels1, tick=T, lwd=1, las = 3)
    axis(2, at=tck2, labels=tick.labels2, tick=T, lwd=1, las = 1)
    abline(v=tck1, h=tck2, lty=6, col = "grey")
    
    title(main = names(x1)[i], font.main = 4)
  }
#  dev.off()
}
