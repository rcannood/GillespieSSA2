ssa_method <- function(
  name,
  class = paste0("SSA_", name),
  params,
  factory
) {
  l <- list(
    name = name,
    params = params,
    factory = factory
  )
  class(l) <- c(class, "SSA_method")
  l
}
