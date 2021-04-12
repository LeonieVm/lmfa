

#' Plot function for: A Generic Convex-Hull-Based Model Selection Method
#'
#' This function is based on the CHull function from the R package multichull
#'
#' @param x  Output from the main function step1() (must be of class lmfa_modelselection).
#' @param col Vector of \code{\link{colors}} used for plots.
#' @param pch Symbol used to indicate selected model(s).
#' @param ... Further arguments for the CHull plot function.


plot.CHull <-
  function(x, col = NULL, pch = NULL, ...){
    hull <- x$Hull
    Solution <- x$Solution
    data <- x$OrigData
    bound <- x$Bound
    if (is.null(col)){col <- 1:3}
    if (is.null(pch)){pch <- 19}
    
    label <- switch(bound,
                    upper="LL",
                    lower="Badness-of-fit")

    pos_vector <- rep(3, nrow(hull))
    pos_vector[1] <- 4
    pos_vector[nrow(hull)] <- 2
    
    plot(data[,1],data[,2],xlab="n_par",ylab=label,col=col[1],lty=1,type="p")
    points(hull[,1],hull[,2],col=col[2],type="b",lty=1)
    points(Solution[,1],Solution[,2],col=col[3],type="p",lty=1,pch=pch)
    text(hull[,1],hull[,2], labels=rownames(hull), cex= 1, pos=pos_vector,col=col[2])
  }