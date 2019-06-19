#' Invoking the stochastic simulation algorithm
#'
#' Main interface function to the implemented \acronym{SSA} methods. Runs a
#' single realization of a predefined system.
#'
#' Substantial improvements in speed and accuracy can be obtained by
#' adjusting the additional (and optional) `ssa` arguments. By default
#' `ssa` (tries to) use conservative default values for the these
#' arguments, prioritizing computational accuracy over computational speed.
#' These default values are, however, **not** fool proof for the
#' approximate methods, and occasionally one will have to hand tweak them in
#' order for a stochastic model to run appropriately.
#'
#' @param initial_state `[named numeric vector]` The initial state to start the simulation with.
#' @param propensity_funs `[character vector]` A character representation of the propensity functions, written in C++.
#'   Check the sections below for more information on how to define the propensity functions.
#' @param nu `[numeric matrix]` The changes in number of individuals in a state (rows) caused
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
#' @param use_singular_optimisation `[logical]` An experimental optimisation. Do not use yet.
#'
#' @return Returns a list object with the following elements:
#'
#' * `time`: `[numeric]` Simulatiom time for each recorded timepoint.
#' * `state`: `[numeric matrix]` The state values for each of the timepoints.
#' * `propensity`: `[numeric matrix]` The propensity values for each of the timepoints.
#' * `buffer`: `[numeric matrix]` The temporary calculation buffer used as part of the propensity functions.
#'
#' @note Selecting the appropriate \acronym{SSA} method is a trade-off between
#'   computational speed, accuracy of the results, and which \acronym{SSA}
#'   actually works for a given scenario. This depends on the characteristics of
#'   the defined system (e.g. number of reaction channels, number of species, and
#'   the absolute and relative magnitude of the propensity functions).
#'   **Not all methods are appropriate for all models.** When selecting a
#'   \acronym{SSA} method all of these factors have to be taken into
#'   consideration. The various tau-leap methods accept a number of additional
#'   arguments. While the default values of these arguments may work for some
#'   scenarios they may have to be adjusted for others. The default values for
#'   the tau-leap methods are conservative in terms of computational speed and
#'   substantial increase in efficiency may be gained by optimizing their values
#'   for a specific system.
#' @section Preparing a run: In order to invoke \acronym{SSA} the stochastic
#'   model needs at least four components, the initial state vector (`initial.state`),
#'   state-change matrix (`nu`), propensity vector (`a`), and the final
#'   time of the simulation (`final.time`). The initial state vector defines the
#'   population sizes in all the states at \eqn{t = 0}, e.g. for a system with two
#'   species `X1` and `X2` where both have an initial population size
#'   of 1000 the initial state vector is defined as `initial.state <-
#'   c(X1 = 1000,X2 = 1000)`. The elements of the vector have to be labelled using
#'   the same notation as the state variables used in the propensity functions.
#'   The state-change matrix defines the change in the number of individuals in
#'   each state (rows) as caused by one reaction of a given type (columns). For
#'   example, the state-change matrix for system with the species \eqn{S_1}{S1}
#'   and \eqn{S_2}{S2} with two reactions \deqn{S_1
#'   \stackrel{c_1}{\longrightarrow} S_2}{S1 --c1--> S2} \deqn{S_2
#'   \stackrel{c_2}{\longrightarrow} 0}{S2 --c2--> 0}
#'
#'   is defined as `nu <- matrix(c(-1,0,+1,-1),nrow = 2,byrow = TRUE)` where
#'   \eqn{c_1}{c1} and \eqn{c_2}{c2} are the per capita reaction probabilities.
#'   The propensity vector, `a`, defines the probabilities that a particular
#'   reaction will occur over the next infinitesimal time interval \eqn{\left[
#'   t,t+dt \right]}{[t,t+dt]}. For example, in the previous example the
#'   propensity vector is defined as `a <- c("c1*X1","c2*X2")`. The
#'   propensity vector consists of character elements of each reaction's
#'   propensity function where each state variable requires the corresponding
#'   named element label in the initial state vector (`initial.state`).
#'
#' @seealso [package], [ssa_direct()], [ssa_etl()], [ssa_btl()]
#'
#' @examples
#' ## Irreversible isomerization
#' ## Large initial population size (X = 1000)
#' \dontrun{
#' params <- c(c = 0.5)
#' initial_state <- c(X = 1000)
#' propensity_funs <- c("A = c * X")
#' nu <- matrix(-1, dimnames = list(names(initial_state), "A"))
#' out <- ssa(
#'   initial_state = initial_state,
#'   propensity_funs = propensity_funs,
#'   nu = nu,
#'   params = params,
#'   final_time = 10,
#'   method = ssa_direct()
#' )
#' ssa_plot(out)
#' }
#'
#' ## Smaller initial population size (X = 100)
#' \dontrun{
#' initial_state <- c(X = 100)
#' out <- ssa(
#'   initial_state = initial_state,
#'   propensity_funs = propensity_funs,
#'   nu = nu,
#'   params = params,
#'   final_time = 10,
#'   method = ssa_direct()
#' )
#' ssa_plot(out)
#' }
#'
#' ## Small initial population size (X = 10)
#' \dontrun{
#' initial_state <- c(X = 10)
#' out <- ssa(
#'   initial_state = initial_state,
#'   propensity_funs = propensity_funs,
#'   nu = nu,
#'   params = params,
#'   final_time = 10,
#'   method = ssa_direct()
#' )
#' ssa_plot(out)
#' }
#'
#' ## Logistic growth
#' \dontrun{
#' params <- c(b = 2, d = 1, K = 1000)
#' initial_state  <- c(N = 500)
#' propensity_funs <- c(
#'   "A = b * N",
#'   "B = (d + (b - d) * N / K) * N"
#' )
#' nu  <- matrix(
#'   c(+1, -1),
#'   ncol = 2,
#'   dimnames = list(names(initial_state), c("A", "B"))
#' )
#' out <- ssa(
#'   initial_state = initial_state,
#'   propensity_funs = propensity_funs,
#'   nu = nu,
#'   params = params,
#'   final_time = 10,
#'   max_walltime = 5,
#'   method = ssa_direct()
#' )
#' ssa_plot(out)
#' }
#'
#' ## Kermack-McKendrick SIR model
#' \dontrun{
#' params <- c(beta = 0.001, gamma = 0.1)
#' initial_state  <- c(S = 500, I = 1, R = 0)
#' propensity_funs <- c(
#'   "A = beta * S * I",
#'   "B = gamma * I"
#' )
#' nu <- matrix(
#'   c(
#'     -1, 0,
#'     +1, -1,
#'     0, +1
#'   ),
#'   nrow = 3,
#'   byrow = TRUE,
#'   dimnames = list(names(initial_state), c("A", "B"))
#' )
#' out <- ssa(
#'   initial_state = initial_state,
#'   propensity_funs = propensity_funs,
#'   nu = nu,
#'   params = params,
#'   final_time = 100,
#'   method = ssa_direct(),
#'   census_interval = 0.01
#' )
#' ssa_plot(out)
#' }
#'
#' ## Lotka predator-prey model
#' \dontrun{
#' params <- c(c1 = 10, c2 = .01, c3 = 10)
#' initial_state <- c(Y1 = 1000, Y2 = 1000)
#' propensity_funs <- c(
#'   "A = c1 * Y1",
#'   "B = c2 * Y1 * Y2",
#'   "C = c3 * Y2"
#' )
#' nu <- matrix(
#'   c(
#'     +1, -1, 0,
#'     0, +1, -1
#'   ),
#'   nrow = 2,
#'   byrow = TRUE,
#'   dimnames = list(names(initial_state), c("A", "B", "C"))
#' )
#' out <- ssa(
#'   initial_state = initial_state,
#'   propensity_funs = propensity_funs,
#'   nu = nu,
#'   params = params,
#'   final_time = 100,
#'   method = ssa_direct(),
#'   census_interval = .01,
#'   verbose = TRUE
#' )
#' ssa_plot(out)
#' }
#'
#' @export
#'
#' @importFrom methods is
#' @importFrom dynutils is_sparse
#' @importFrom tibble tibble
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
  console_interval = 1,
  use_singular_optimisation = FALSE
) {
  # check parameters
  # TODO: check all params
  assert_that(
    is.numeric(initial_state),
    !is.null(names(initial_state)),
    is.matrix(nu) || dynutils::is_sparse(nu),
    is.numeric(final_time), final_time >= 0,
    is.numeric(console_interval), console_interval >= 0,
    is.numeric(census_interval), census_interval >= 0,
    is.logical(verbose),
    is(method, "fastgssa::ssa_method"),
    length(initial_state) == nrow(nu),
    is.numeric(params),
    length(params) == 0 || !is.null(names(params)),
    !any(c("state", "time", "params") %in% names(params)),
    !any(c("state", "time", "params") %in% names(initial_state)),
    !any(duplicated(c(names(initial_state), names(params))))
  )

  # compile propensity functions if this has not been done already
  comp_funs <-
    if (is.character(propensity_funs)) {
      assert_that(length(propensity_funs) == ncol(nu))
      compile_propensity_functions(
        propensity_funs = propensity_funs,
        reaction_ids = colnames(nu),
        state_ids = rownames(nu),
        params = params,
        hardcode_params = hardcode_params
      )
    } else {
      propensity_funs
    }

  # check propensity functions
  assert_that(
    is(comp_funs, "fastgssa::propensity_functions")
  )

  # create ssa algorithm instance
  ssa_alg <- method$factory()

  # run SSA
  output <- simulate(
    propensity_funs = comp_funs$pointer,
    ssa_alg = ssa_alg,
    initial_state = initial_state,
    params = params,
    nu = as.matrix(nu),
    final_time = final_time,
    census_interval = census_interval,
    buffer_size = comp_funs$buffer_size,
    max_walltime = max_walltime,
    stop_on_neg_state = stop_on_neg_state,
    verbose = verbose,
    console_interval = console_interval,
    use_singular_optimisation = use_singular_optimisation
  )

  # set colnames of objects
  colnames(output$state) <- rownames(nu)
  colnames(output$propensity) <- colnames(nu)
  colnames(output$buffer) <- comp_funs$buffer_ids

  output
}
