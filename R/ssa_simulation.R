# boilerplate function but allows for parameter name auto completion
create_simulation <- function(
  compiled_reactions,
  params = NULL,
  method = ssa_exact(),
  initial_state,
  census_interval = 0,
  log_propensity = FALSE,
  log_firings = FALSE,
  log_buffer = FALSE,
  stop_on_neg_state = TRUE,
  final_time,
  max_walltime = Inf,
  sim_name = NA_character_,
  verbose = FALSE,
  console_interval = 1
) {
  # create new simulation object
  sim <- new(SSA_simulation)

  # set propensity functions
  sim$initialise(
    compiled_reactions$num_functions,
    compiled_reactions$functions_pointer,
    params,
    compiled_reactions$buffer_size,
    initial_state,
    compiled_reactions$state_change@i,
    compiled_reactions$state_change@p,
    compiled_reactions$state_change@x,
    census_interval,
    log_propensity,
    log_firings,
    log_buffer,
    stop_on_neg_state,
    final_time,
    max_walltime,
    sim_name,
    verbose,
    console_interval
  )
  configure_method(method, sim)
  sim$reset()

  sim
}


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
#' `ssa` uses conservative parameters (o.a. [ssa_exact()]) which prioritise
#' computational accuracy over computational speed.
#'
#' Approximate methods ([ssa_etl()] and [ssa_btl()]) are not fool proof!
#' Some tweaking might be required for a stochastic model to run appropriately.
#'
#' @param initial_state `[named numeric vector]` The initial state to start the simulation with.
#' @param reactions A list of reactions, see [reaction()].
#' @param final_time `[numeric]` The final simulation time.
#' @param params `[named numeric vector]` Constant parameters to be used in the propensity functions.
#' @param method `[ssa_method]`] Which SSA algorithm to use. Must be one of: [ssa_exact()],
#'   [ssa_btl()], or [ssa_etl()].
#' @param census_interval `[numeric]` The approximate interval between recording the state of the system.
#'   Setting this parameter to `0` will cause each state to be recorded, and
#'   to `Inf` will cause only the end state to be recorded.
#' @param stop_on_neg_state `[logical]` Whether or not to stop the simulation when
#'   the a negative value in the state has occured. This can occur, for instance, in the [ssa_etl()]
#'   method.
#' @param max_walltime `[numeric]` The maximum duration (in seconds) that the
#'   simulation is allowed to run for before terminated.
#' @param log_propensity `[logical]` Whether or not to store the propensity values at each census.
#' @param log_firings `[logical]` Whether or not to store number of firings of each reaction between censuses.
#' @param log_buffer `[logical]` Whether or not to store the buffer at each census.
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
#'     method = ssa_exact(),
#'     final_time = 5,
#'     census_interval = .001,
#'     verbose = TRUE
#'   )
#' autoplot.ssa(out)
#'
#' @export
#'
#' @importFrom methods is
#' @importFrom dynutils is_sparse inherit_default_params
#' @importFrom tibble tibble
#' @importFrom Matrix Matrix
#' @importFrom purrr is_scalar_double is_scalar_logical is_scalar_integer is_scalar_character
ssa <- dynutils::inherit_default_params(
  list(create_simulation),
  function(
    initial_state,
    reactions,
    final_time,
    params,
    method = ssa_exact(),
    census_interval,
    stop_on_neg_state,
    max_walltime,
    log_propensity,
    log_firings,
    log_buffer,
    verbose,
    console_interval,
    sim_name
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
      is(method, "SSA_method"),

      # vector arguments
      is_scalar_double(final_time), final_time >= 0,
      is_scalar_double(census_interval), census_interval >= 0,
      is_scalar_double(max_walltime), max_walltime >= 0,
      is_scalar_double(console_interval), console_interval >= 0,
      is_scalar_logical(log_propensity),
      is_scalar_logical(log_firings),
      is_scalar_logical(log_buffer),
      is_scalar_logical(verbose),
      is_scalar_character(sim_name)
    )

    # compile propensity functions if this has not been done already
    compiled_reactions <-
      if (is.list(reactions) && !is(reactions, "gillespie::reactions")) {
        compile_reactions(
          reactions = reactions,
          state_ids = names(initial_state),
          params = params
        )
      } else {
        reactions
      }

    assert_that(is(compiled_reactions, "gillespie::reactions"))

    # create new simulation object
    sim <- create_simulation(
      compiled_reactions = compiled_reactions,
      params = params,
      method = method,
      initial_state = initial_state,
      census_interval = census_interval,
      log_propensity = log_propensity,
      log_firings = log_firings,
      log_buffer = log_buffer,
      stop_on_neg_state = stop_on_neg_state,
      final_time = final_time,
      max_walltime = max_walltime,
      sim_name = sim_name,
      verbose = verbose,
      console_interval = console_interval
    )

    # run simulation
    sim$run()

    output <- list(
      time = sim$output_time,
      state = sim$output_state,
      propensity = sim$output_propensity,
      firings = sim$output_firings,
      buffer = sim$output_buffer,
      stats = sim$get_statistics(),
      name = sim$sim_name
    )

    rm(sim)

    # set colnames of objects
    colnames(output$state) <- names(initial_state)
    if (log_propensity) {
      colnames(output$propensity) <- compiled_reactions$reaction_ids
    }
    if (log_buffer) {
      colnames(output$buffer) <- compiled_reactions$buffer_ids
    }
    if (log_firings) {
      colnames(output$firings) <- compiled_reactions$reaction_ids
    }

    class(output) <- c("ssa", "list")

    output
  }
)
