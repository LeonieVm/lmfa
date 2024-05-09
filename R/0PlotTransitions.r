#' Transition plot for a specific subject
#'
#'
#'
#'
#'
#'
#' @param x  Output from the main function step3() (must be of class lmfa_step3).
#' @param identifier The name of the column containing the subject identifiers (must be a single character).
#' @param id Identification number of the subject for which the transition pattern should be plotted.
#' @param ... Further arguments for the default S3 plot method.
#' @examples
#' \dontrun{
#' plot(transitionmodel)
#' }
#' @export
plot.lmfa_step3 <- function(x, identifier, id, ...) {
        if (missing(x)) stop("argument data is missing, with no default")
        if (missing(identifier)) stop("argument identifier is missing, with no default")
        if (missing(id)) stop("argument id is missing, with no default")
        dataPlot <- x$data
        if (is.element(id, unique(dataPlot[, identifier])) == FALSE) stop("id is not part of the specified identifier column")

        plot(dataPlot[dataPlot[, identifier] == id, "Modal"],
                type = "b",
                xlab = "Measurement Occasion",
                ylab = "State Membership",
                yaxt = "n",
                main = paste("Subject", id)
        )
        axis(side = 2, tick = seq(1, x$n_state, 1), at = seq(1, x$n_state, 1))
}