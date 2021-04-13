#' Enumerate the Combinations of the Elements of a Vector
#'
#' This function was copied and adapted from R package gtools: https://github.com/cran/gtools/blob/master/R/combinations.R.
#'
#'
#' @param n Size of the source vector.
#' @param r Size of the target vectors.
#' @param v Source vector. Defaults to 1:n.
#'
#'
#' @return Returns matrix of unique combinations
#'
#' @examples
#' \dontrun{
#' 
#' }
#' @noRd

#Function adapted from R package combinations
#https://github.com/cran/gtools/blob/master/R/combinations.R
combinations_K_F_k <- function(n, r, v = 1:n)
{
  repeats.allowed <- TRUE
  set = TRUE
  if(mode(n) != "numeric" || length(n) != 1 
     || n < 1 || (n %% 1) != 0) stop("bad value of n") 
  if(mode(r) != "numeric" || length(r) != 1 
     || r < 1 || (r %% 1) != 0) stop("bad value of r") 
  if(!is.atomic(v) || length(v) < n) 
    stop("v is either non-atomic or too short")
  if( (r > n) & repeats.allowed==FALSE)
    stop("r > n and repeats.allowed=FALSE")
  if(set) {
    v <- unique(sort(v,decreasing = TRUE))
    if (length(v) < n) stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  ## Inner workhorse
  if(repeats.allowed)
    sub <- function(n, r, v)
    { 
      if(r == 0) v0 else
        if(r == 1) matrix(v, n, 1) else
          if(n == 1) matrix(v, 1, r) else
            rbind( cbind(v[1], Recall(n, r-1, v)),
                   Recall(n-1, r, v[-1]))
    }
  sub(n, r, v[1:n])
}