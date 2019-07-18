ssa_method <- function(name, params, factory) {
  l <- lst(
    name,
    params,
    factory
  )
  class(l) <- "gillespie::ssa_method"
  l
}
