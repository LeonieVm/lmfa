#' Evaluating for which individuals invariance holds across all measurement occasions
#'
#' \code{invariance} evaluates for which individuals invariance holds across all measurement occasions based on the \code{step3} output.
#'
#'
#'
#'
#'
#'
#' @param model  Output from the main function step3() (must be of class lmfa_step3).
#' @param identifier The name of the column containing the subject identifiers (must be a single character).
#
#'
#'
#' @examples
#' \dontrun{
#' invariance(model, identifier)
#' }
#' @export
invariance <- function(model, identifier) {
    if (missing(model)) stop("argument data is missing, with no default")
    if (missing(identifier)) stop("argument identifier is missing, with no default")
    if (class(model) != "lmfa_step3") stop("model must be of class lmfa_step3")

    dataInvariance <- model$data
    invariance_time <- c()
    statemembership <- c()
    subjects <- unique(dataInvariance[, identifier])
    for (i in 1:length(subjects)) {
        subject_number <- unique(dataInvariance[, identifier])[i]

        # if within perosn invariance holds...
        if (length(unique(dataInvariance[dataInvariance[, identifier] == subject_number, "Modal"])) == 1) {
            invariance_time <- c(invariance_time, subject_number)
            statemembership <- c(statemembership, unique(dataInvariance[dataInvariance[, identifier] == subject_number, "Modal"]))
        }
    }
    # if it holds for anyone...
    if (!is.null(statemembership)) {
        statemembership <- as.matrix(statemembership, ncol = 1)
        rownames(statemembership) <- invariance_time


        invariance_subjects <- list()
        for (i in 1:length(sort(unique(statemembership)))) {
            invariance_subjects[[i]] <- as.numeric(row.names(statemembership)[which(statemembership == i)])
        }

        names(invariance_subjects) <- paste("S", sort(unique(statemembership)), sep = "")

        for (i in 1:length(invariance_subjects)) {
            cat(names(invariance_subjects)[i])
            cat("\n")

            cat(invariance_subjects[[i]])
            cat("\n")
            cat("\n")
        }

        # if it does not hold for anyone...
    } else {
        x_n_state <- 1:model$n_state
        invariance_subjects <- list()
        for (i in 1:length(x_n_state)) {
            invariance_subjects[[x_n_state[i]]] <- NA
        }

        names(invariance_subjects) <- paste("S", 1:model$n_state, sep = "")
        for (i in 1:length(invariance_subjects)) {
            cat(names(invariance_subjects)[i])
            cat("\n")

            cat(invariance_subjects[[i]])
            cat("\n")
            cat("\n")
        }
    }
}