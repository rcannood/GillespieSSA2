#' Port GillespieSSA parameters to GillespieSSA2
#'
#' This is a helper function to tranform GillesieSSA-style
#' paramters to GillespieSSA2.
#'
#' @param x0 The `x0` parameter of [GillespieSSA::ssa()].
#' @param a The `a` parameter of [GillespieSSA::ssa()].
#' @param nu The `nu` parameter of [GillespieSSA::ssa()].
#'
#' @return A set of [reaction()]s to be used by [ssa()].
#'
#' @export
#'
#' @examples
#' x0  <- c(Y1 = 1000, Y2 = 1000)
#' a   <- c("c1*Y1","c2*Y1*Y2","c3*Y2")
#' nu  <- matrix(c(+1,-1,0,0,+1,-1),nrow=2,byrow=TRUE)
#' port_reactions(x0, a, nu)
port_reactions <- function(x0, a, nu) {
  map(
    seq_along(a),
    function(i) {
      ix <- which(nu[,i] != 0)
      effect <- set_names(nu[ix, i], names(x0)[ix])
      reaction(a[[i]], effect)
    }
  )
}
