ssa_method <- function(
  name,
  class = paste0("SSA_", name),
  params,
  factory
) {
  l <- lst(
    name,
    params,
    factory
  )
  class(l) <- c(class, "SSA_method")
  l
}
