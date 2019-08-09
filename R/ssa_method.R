# A wrapper object for creating SSA method descriptors
#
# @param name The name of the SSA method
# @param class The class name of the return object
# @param params The parameters to be passed to the algorithm
# @param factory A factory function that creates a pointer to a compiled version of the algorithm.
ssa_method <- function(
  name,
  class = paste0("SSA_", name),
  params,
  factory
) {
  l <- list(
    name = name,
    class = class,
    params = params,
    factory = function() do.call(factory, params)
  )
  class(l) <- c(class, "SSA_method")
  l
}

#' @rdname print_ssa
#' @export
print.SSA_method <- function(x, ...) {
  params <- x[["params"]]
  param_str <-
    if (length(params) > 0) {
      paste0(names(params), " = ", params, collapse = ", ")
    } else {
      ""
    }
  cat(
    x[["class"]], "(", param_str, ")\n",
    sep = ""
  )
}
