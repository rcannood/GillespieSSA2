#' Simple plotting of ssa output
#'
#' Provides basic functionally for simple and quick time series plot of simulation output from [ssa()].
#'
#' @param ssa_out Data object returned by [ssa()].
#' @param state Whether or not to plot the state values.
#' @param propensity Whether or not to plot the propensity values.
#' @param firings Whether or not to plot the reaction firings values.
#' @param buffer Whether or not to plot the buffer values.
#' @param geom Which geom to use, must be one of `"point"`, `"step"`.
#'
#' @importFrom tidyr gather
#' @importFrom dplyr bind_rows mutate
#' @export
plot_ssa <- function(ssa_out, state = TRUE, propensity = FALSE, buffer = FALSE, firings = FALSE, geom = c("point", "step")) {
  requireNamespace("ggplot2")
  geom <- match.arg(geom)

  var <- value <- time <- type <- NULL # satisfying r check

  # collect data to be plotted
  df <- NULL
  var_names <- c()

  if (state) {
    df <-
      data.frame(time = ssa_out$time, ssa_out$state, type = "state", stringsAsFactors = FALSE) %>%
      gather(var, value, -time, -type) %>%
      bind_rows(df)
    var_names <- c(var_names, colnames(ssa_out$state))
  }

  if (propensity) {
    df <-
      data.frame(time = ssa_out$time, ssa_out$propensity, type = "propensity", stringsAsFactors = FALSE) %>%
      gather(var, value, -time, -type) %>%
      bind_rows(df)
    var_names <- c(var_names, colnames(ssa_out$propensity))
  }

  if (firings) {
    df <-
      data.frame(time = ssa_out$time, ssa_out$firings, type = "firings", stringsAsFactors = FALSE) %>%
      gather(var, value, -time, -type) %>%
      bind_rows(df)
    var_names <- c(var_names, colnames(ssa_out$firings))
  }

  if (buffer && ncol(ssa_out$buffer) > 0) {
    df <-
      data.frame(time = ssa_out$time, ssa_out$buffer, type = "buffer", stringsAsFactors = FALSE) %>%
      gather(var, value, -time, -type) %>%
      bind_rows(df)
    var_names <- c(var_names, colnames(ssa_out$buffer))
  }

  # change levels of the var
  df <- df %>%
    mutate(
      type = factor(type, levels = c("state", "propensity", "firings", "buffer")),
      var = factor(var, levels = var_names)
    )

  # create plot
  g <-
    ggplot2::ggplot(df, ggplot2::aes_string("time", "value", colour = "var")) +
    ggplot2::facet_wrap(~type, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = as.character(ssa_out$stats$sim_name) %|% "SSA Simulation",
      subtitle = paste0(
        ssa_out$stats$method, ", ",
        round(ssa_out$stats$walltime_elapsed, 3), " sec, ",
        ssa_out$stats$num_steps, " steps"
      )
    ) +
    ggplot2::theme_bw()

  # add geom depending on parameter
  if (geom == "point") {
    g <- g + ggplot2::geom_point(size = .2)
  } else if (geom == "step") {
    g <- g + ggplot2::geom_step()
  }

  g
}
