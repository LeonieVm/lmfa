#' A Generic Convex-Hull-Based Model Selection Method
#'
#' \code{chull_lmfa} is based on \code{CHull} from the R package multichull.
#'
#'
#' @param x Model-selection output from \code{step1} (must be of class lmfa_modelselection).
#' @param ... Further arguments for the CHull plot function.
#'
#'
#'
#' @return Returns CHull output.
#'
#' @examples
#' \dontrun{
#' chull_lmfa(x)
#' }
#' @export

chull_lmfa <- function(x, ...) {
  if (missing(x)) stop("argument data is missing, with no default")
  if (class(x) != "lmfa_modelselection") stop("x must be of class lmfa_modelselection")
  # if(!is.numeric(PercentageFit)) stop("PercentageFit must be a single scalar")
  # if(length(PercentageFit)>1) stop("PercentageFit must be a single scalar")
  PercentageFit <- 0.0
  object <- x
  # the following is just the summary for the model selection object. However, because of the added non-convergence note, it is not possible to simply extract the statistics. Thus, I compute them here as well.
  #---------------------------------------------------------------------------------------------------------#
  modelcomparison <- matrix(NA, nrow = length(object), ncol = 6)
  colnames(modelcomparison) <- c("LL", "BIC", "convergence", "n_par", "n_fact", "n_state")
  modelnames <- NULL
  for (i in 1:length(object)) {
    modelcomparison[i, 1] <- object[[i]]$LL
    modelcomparison[i, 2] <- object[[i]]$BIC
    modelcomparison[i, 3] <- object[[i]]$convergence
    modelcomparison[i, 4] <- object[[i]]$n_par
    modelcomparison[i, 5] <- sum(object[[i]]$n_fact)
    modelcomparison[i, 6] <- object[[i]]$n_state
    modelnames <- c(modelnames, paste("[", paste(object[[i]]$n_fact, collapse = ""), "]", sep = ""))
  }
  rownames(modelcomparison) <- modelnames
  modelcomparison <- modelcomparison[order(modelcomparison[, "BIC"]), ]

  #---------------------------------------------------------------------------------------------------------#

  objecModelselection <- (modelcomparison[, c(1:4)])

  CHullInput <- objecModelselection
  CHullInput <- CHullInput[CHullInput[, "convergence"] == 1, drop = FALSE, ]
  if (nrow(CHullInput) > 1) {
    fitCHull <- CHull(CHullInput[, c("n_par", "LL")], bound = "upper", PercentageFit = PercentageFit)
    if (length(fitCHull) > 1) {
      plot(fitCHull,
        col = c("black", "black", "red"), pch = 21,
        bg = "white", ylab = "LL", xlab = "n_par", ...
      )

      Hull <- fitCHull$Hull
      Solution <- fitCHull$Solution
      colnames(Hull) <- c("n_par", "LL", "st")
      colnames(Solution) <- c("n_par", "LL")

      cat("\n")
      cat(paste("Models on the upper boundary of the CHull:"), "\n")
      cat("\n")
      print(Hull)

      cat("\n")
      cat(paste("Selected model(s):"), "\n")
      cat("\n")
      print(Solution)


      # calculate numerators
      sumComplexModels <- 0
      if (nrow(Hull) > 2) {
        store <- matrix(NA, nrow = (nrow(Hull) - 2), ncol = 2)
        for (i in 2:(nrow(Hull) - 1)) {
          st <- Hull[i, 3]
          numerator <- ((Hull[i, 2] - Hull[i - 1, 2]) / (Hull[i, 1] - Hull[i - 1, 1]))
          store[i - 1, 1] <- st
          store[i - 1, 2] <- numerator
        }
        sumComplexModels <- sum(abs(store[, 1] - store[, 2]) < 10)
        if (sumComplexModels > 0) {
          cat("\n")
          cat(paste("Note 1: The least and most complex models cannot be selected."), "\n")
          cat(paste("  Therefore, it is advisable to also visually inspect the CHull plot."), "\n")

          cat("\n")
          cat(paste("Note 2: The st value(s) of the best model(s) might be artificially"), "\n")
          cat(paste("  inflated. Therefore, it is advisable to also consider the next best model(s)."), "\n")
        } else {
          cat("\n")
          cat(paste("Note: The least and most complex models cannot be selected."), "\n")
          cat(paste("  Therefore, it is advisable to also visually inspect the CHull plot."), "\n")
        }
      } else {
        cat("\n")
        cat(paste("Note: The least and most complex models cannot be selected."), "\n")
        cat(paste("  Therefore, it is advisable to also visually inspect the CHull plot."), "\n")
      }

      output <- (list(
        fitCHull = fitCHull,
        solution = Solution,
        chull = Hull,
        sumComplexModels = sumComplexModels
      ))

      class(output) <- "lmfa_chull"
      invisible(output)
    }
  } else {
    cat("\n")
    cat(paste("Not enough data points available to compute the convex hull"), "\n")
    cat("\n")
  }
}
