#' Invoking the stochastic simulation algorithm
#'
#' Main interface function to the implemented \acronym{SSA} methods. Runs a
#' single realization of a predefined system. For a detailed explanation
#' on how to set up your first \acronym{SSA} system, check the vignette
#' available on [GitHub](https://github.com/dynverse/fastgssa/tree/master/vignettes/preparing_a_run.md)
#' or using `vignette("preparing_a_run", package = "fastgssa")`.
#'
#' Substantial improvements in speed and accuracy can be obtained by
#' adjusting the additional (and optional) `ssa` arguments. By default
#' `ssa` uses conservative parameters (o.a. [ssa_direct()]) which prioritise
#' computational accuracy over computational speed.
#'
#' Approximate methods ([ssa_etl()] and [ssa_btl()]) are not fool proof!
#' Some tweaking might be required for a stochastic model to run appropriately.
#'
#' @param initial_state `[named numeric vector]` The initial state to start the simulation with.
#' @param propensity_funs `[character vector]` A character representation of the propensity functions, written in C++.
#' @param nu `[(sparse) integer matrix]` The changes in number of individuals in a state (rows) caused
#'   by a single reaction (columns).
#' @param final_time `[numeric]` The final simulation time.
#' @param params `[named numeric vector]` Constant parameters to be used in the propensity functions.
#' @param method `[SSA]`] Which SSA algorithm to use. Must be one of: [ssa_direct()],
#'   [ssa_btl()], or [ssa_etl()].
#' @param census_interval `[numeric]` The approximate interval between recording the state of the system.
#'   Setting this parameter to `0` will cause each state to be recorded, and
#'   to `Inf` will cause only the end state to be recorded.
#' @param stop_on_neg_state `[logical]` Whether or not to stop the simulation when
#'   the a negative value in the state has occured. This can occur, for instance, in the [ssa_etl()]
#'   method.
#' @param max_walltime `[numeric]` The maximum duration (in seconds) that the
#'   simulation is allowed to run for before terminated.
#' @param hardcode_params `[logical]` Whether or not to hardcode the values of `params` in the compilation of the `propensity_funs`.
#'   Setting this to `TRUE` will result in a minor sacrifice in accuracy for a minor increase in performance.
#' @param verbose `[logical]` If `TRUE`, intermediary information pertaining to the simulation will be displayed.
#' @param console_interval `[numeric]` The approximate interval between intermediary information outputs.
#'
#' @return Returns a list object with the following elements:
#'
#' * `time`: `[numeric]` Simulatiom time for each recorded timepoint.
#' * `state`: `[numeric matrix]` The state values for each of the timepoints.
#' * `propensity`: `[numeric matrix]` The propensity values for each of the timepoints.
#' * `buffer`: `[numeric matrix]` The temporary calculation buffer used as part of the propensity functions.
#'
#' @seealso [fastgssa] for a high level explanation of the package
#'
#' @examples
#' initial_state <- c(x_prey = 1000, x_predators = 1000)
#' nu <- matrix(
#'   c(
#'     +1, -1, 0,
#'     0, +1, -1
#'   ),
#'   nrow = 2,
#'   byrow = TRUE,
#'   dimnames = list(
#'     c("x_prey", "x_predators"),
#'     c("p_prey_up", "p_predation", "p_pred_down")
#'   )
#' )
#' propensity_funs <- c(
#'   "p_prey_up = c1 * x_prey",
#'   "p_predation = c2 * x_prey * x_predators",
#'   "p_pred_down = c3 * x_predators"
#' )
#' params <- c(c1 = 10, c2 = 0.01, c3 = 10)
#'
#' out <-
#'   ssa(
#'     initial_state = initial_state,
#'     propensity_funs = propensity_funs,
#'     nu = nu,
#'     params = params,
#'     method = ssa_direct(),
#'     final_time = 5,
#'     census_interval = .001,
#'     verbose = TRUE
#'   )
#' ssa_plot(out)
#'
#' @export
#'
#' @importFrom methods is
#' @importFrom dynutils is_sparse
#' @importFrom tibble tibble
#' @importFrom Matrix Matrix
ssa <- function(
  initial_state,
  propensity_funs,
  nu,
  final_time,
  params = NULL,
  method = ssa_direct(),
  census_interval = 0,
  stop_on_neg_state = TRUE,
  max_walltime = Inf,
  hardcode_params = FALSE,
  verbose = FALSE,
  console_interval = 1
) {
  # implement these stats:
  # https://github.com/cran/GillespieSSA/blob/master/man/ssa.Rd#L71

  if (is.matrix(nu)) {
    nu <- Matrix::Matrix(nu, sparse = TRUE)
  }

  # check parameters
  assert_that(
    # initial state
    is.numeric(initial_state),
    !is.null(names(initial_state)),
    !any(c("state", "time", "params") %in% names(initial_state)),

    # nu
    dynutils::is_sparse(nu),
    length(initial_state) == nrow(nu),

    # params
    is.numeric(params),
    length(params) == 0 || !is.null(names(params)),
    !any(c("state", "time", "params") %in% names(params)),
    !any(duplicated(c(names(initial_state), names(params)))),

    # method
    is(method, "fastgssa::ssa_method"),

    # vector arguments
    is.numeric(final_time), length(final_time) == 1, final_time >= 0,
    is.numeric(census_interval), length(census_interval) == 1, census_interval >= 0,
    is.logical(stop_on_neg_state), length(stop_on_neg_state) == 1,
    is.numeric(max_walltime), length(max_walltime) == 1, max_walltime >= 0,
    is.logical(hardcode_params), length(hardcode_params) == 1,
    is.numeric(console_interval), length(console_interval) == 1, console_interval >= 0,
    is.logical(verbose), length(verbose) == 1
  )

  # compile propensity functions if this has not been done already
  comp_funs <-
    if (is.character(propensity_funs)) {
      assert_that(length(propensity_funs) == ncol(nu))
      compile_propensity_functions(
        propensity_funs = propensity_funs,
        state_ids = rownames(nu),
        params = params,
        hardcode_params = hardcode_params
      )
    } else {
      propensity_funs
    }

  # check propensity functions
  assert_that(
    is(comp_funs, "fastgssa::propensity_functions"),
    length(comp_funs$reaction_ids) == ncol(nu)
  )

  # run SSA
  output <- simulate(
    propensity_funs = comp_funs$functions_pointer,
    num_functions = comp_funs$num_functions,
    ssa_alg = method$factory(),
    initial_state = initial_state,
    params = params,
    nu_i = nu@i,
    nu_p = nu@p,
    nu_x = nu@x,
    final_time = final_time,
    census_interval = census_interval,
    buffer_size = comp_funs$buffer_size,
    max_walltime = max_walltime,
    stop_on_neg_state = stop_on_neg_state,
    verbose = verbose,
    console_interval = console_interval
  )

  # set colnames of objects
  colnames(output$state) <- names(initial_state)
  colnames(output$propensity) <- comp_funs$reaction_ids
  colnames(output$buffer) <- comp_funs$buffer_ids

  output
}
