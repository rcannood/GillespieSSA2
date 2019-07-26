#' Port GillespieSSA parameters to GillespieSSA2
#'
#' This is a helper function to tranform GillesieSSA-style
#' paramters to GillespieSSA2.
#'
#' @param a The `a` parameter of [GillespieSSA::ssa()].
#' @param nu The `nu` parameter of [GillespieSSA::ssa()].
#'
#' @export
#'
#' @examples
#' a   <- c("c1*Y1","c2*Y1*Y2","c3*Y2")
#' nu  <- matrix(c(+1,-1,0,0,+1,-1),nrow=2,byrow=TRUE)
#' port_reactions(a, nu)
port_reactions <- function(a, nu) {
  map(
    seq_along(a),
    function(i) {
      ix <- which(nu[,i] != 0)
      effect <- set_names(nu[ix, i], names(x0)[ix])
      reaction(a[[i]], effect)
    }
  )
}
