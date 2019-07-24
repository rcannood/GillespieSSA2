test_ssa_method_step <- function(
  ssa_method,
  state,
  propensity,
  nu
) {
  test_ssa_method_cpp(
    ssa_method$factory(),
    state,
    propensity,
    nu@i,
    nu@p,
    nu@x
  )
}
