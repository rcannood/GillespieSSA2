ssa_method <- function(
  name,
  class = paste0("SSA_", name),
  params
) {
  l <- lst(
    name,
    params
  )
  class(l) <- c(class, "SSA_method")
  l
}

configure_method <- function(
  method,
  simulation
) {
  UseMethod("configure_method")
}
