#' Invoking the stochastic simulation algorithm
#'
#' Main interface function to the implemented \acronym{SSA} methods. Runs a
#' single realization of a predefined system. For a detailed explanation
#' on how to set up your first \acronym{SSA} system, check the vignette
#' available on [GitHub](https://github.com/dynverse/gillespie/tree/master/vignettes/preparing_a_run.md)
#' or using `vignette("preparing_a_run", package = "gillespie")`.
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
#' @param reactions A list of reactions, see [reaction()].
#' @param final_time `[numeric]` The final simulation time.
#' @param params `[named numeric vector]` Constant parameters to be used in the propensity functions.
#' @param method `[ssa_method]`] Which SSA algorithm to use. Must be one of: [ssa_direct()],
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
#' @param sim_name `[character]` An optional name for the simulation.
#'
#' @return Returns a list object with the following elements:
#'
#' * `time`: `[numeric]` Simulatiom time for each recorded timepoint.
#' * `state`: `[numeric matrix]` The state values for each of the timepoints.
#' * `propensity`: `[numeric matrix]` The propensity values for each of the timepoints.
#' * `buffer`: `[numeric matrix]` The temporary calculation buffer used as part of the propensity functions.
#'
#' @seealso [gillespie] for a high level explanation of the package
#'
#' @examples
#' initial_state <- c(prey = 1000, predators = 1000)
#' params <- c(c1 = 10, c2 = 0.01, c3 = 10)
#' reactions <- list(
#'   #        ↓ propensity function      ↓ effects                        ↓ name for reaction
#'   reaction(~c1 * prey,                c(prey = +1),                    name = "prey_up"),
#'   reaction(~c2 * prey * predators,    c(prey = -1, predators = +1),    name = "predation"),
#'   reaction(~c3 * predators,           c(predators = -1),               name = "pred_down")
#' )

#'
#' out <-
#'   ssa(
#'     initial_state = initial_state,
#'     reactions = reactions,
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
  reactions,
  final_time,
  params = NULL,
  method = ssa_direct(),
  census_interval = 0,
  stop_on_neg_state = TRUE,
  max_walltime = Inf,
  hardcode_params = FALSE,
  verbose = FALSE,
  console_interval = 1,
  sim_name = NA_character_
) {

  # check parameters
  assert_that(
    # initial state
    is.numeric(initial_state),
    !is.null(names(initial_state)),
    !any(c("state", "time", "params") %in% names(initial_state)),

    # params
    is.numeric(params),
    length(params) == 0 || !is.null(names(params)),
    !any(c("state", "time", "params") %in% names(params)),
    !any(duplicated(c(names(initial_state), names(params)))),

    # method
    is(method, "gillespie::ssa_method"),

    # vector arguments
    is.numeric(final_time), length(final_time) == 1, final_time >= 0,
    is.numeric(census_interval), length(census_interval) == 1, census_interval >= 0,
    is.logical(stop_on_neg_state), length(stop_on_neg_state) == 1,
    is.numeric(max_walltime), length(max_walltime) == 1, max_walltime >= 0,
    is.logical(hardcode_params), length(hardcode_params) == 1,
    is.numeric(console_interval), length(console_interval) == 1, console_interval >= 0,
    is.logical(verbose), length(verbose) == 1,
    is.character(sim_name)
  )

  # compile propensity functions if this has not been done already
  compiled_reactions <-
    if (is.list(reactions) && !is(reactions, "gillespie::reactions")) {
      compile_reactions(
        reactions = reactions,
        state_ids = names(initial_state),
        params = params,
        hardcode_params = hardcode_params
      )
    } else {
      reactions
    }

  assert_that(is(compiled_reactions, "gillespie::reactions"))

  # run SSA
  output <- simulate(
    propensity_funs = compiled_reactions$functions_pointer,
    num_functions = compiled_reactions$num_functions,
    ssa_alg = method$factory(),
    initial_state = initial_state,
    params = params,
    nu_i = compiled_reactions$state_change@i,
    nu_p = compiled_reactions$state_change@p,
    nu_x = compiled_reactions$state_change@x,
    final_time = final_time,
    census_interval = census_interval,
    buffer_size = compiled_reactions$buffer_size,
    max_walltime = max_walltime,
    stop_on_neg_state = stop_on_neg_state,
    verbose = verbose,
    console_interval = console_interval,
    sim_name = sim_name
  )

  # set colnames of objects
  colnames(output$state) <- names(initial_state)
  colnames(output$propensity) <- compiled_reactions$reaction_ids
  colnames(output$buffer) <- compiled_reactions$buffer_ids

  output
}
