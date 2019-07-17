#' Simple plotting of ssa output
#'
#' Provides basic functionally for simple and quick time series plot of simulation output from [ssa()].
#'
#' @param ssa_out Data object returned by [ssa()].
#' @param state Whether or not to plot the state values.
#' @param propensity Whether or not to plot the propensity values.
#' @param buffer Whether or not to plot the buffer values.
#'
#' @importFrom tidyr gather
#' @importFrom dplyr bind_rows
#' @export
ssa_plot <- function(ssa_out, state = TRUE, propensity = FALSE, buffer = FALSE) {
  requireNamespace("ggplot2")

  var <- value <- time <- type <- NULL # satisfying r check
  df <-
    bind_rows(
      if (state) {
        data.frame(time = ssa_out$time, ssa_out$state, type = "state", stringsAsFactors = FALSE) %>%
          gather(var, value, -time, -type)
      } else {
        NULL
      },
      if (propensity) {
        data.frame(time = ssa_out$time, ssa_out$propensity, type = "propensity", stringsAsFactors = FALSE) %>%
          gather(var, value, -time, -type)
      } else {
        NULL
      },
      if (buffer && ncol(ssa_out$buffer) > 0) {
        data.frame(time = ssa_out$time, ssa_out$buffer, type = "buffer", stringsAsFactors = FALSE) %>%
          gather(var, value, -time, -type)
      } else {
        NULL
      }
    )
  ggplot2::ggplot(df) +
    ggplot2::geom_path(ggplot2::aes_string("time", "value", colour = "var")) +
    ggplot2::facet_wrap(~type, ncol = 1, scales = "free_y") +
    ggplot2::labs(subtitle = paste0(
      ssa_out$stats$method, ", ",
      round(ssa_out$stats$walltime_elapsed, 2), " sec, ",
      ssa_out$stats$num_steps, " steps"
    ))
}
