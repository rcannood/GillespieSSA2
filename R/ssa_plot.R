#' Simple plotting of ssa output
#'
#' Provides basic functionally for simple and quick time series plot of simulation output from [ssa()].
#'
#' @param ssa_out Data object returned by [ssa()].
#'
#' @importFrom tidyr gather
#' @importFrom dplyr bind_rows
#' @export
ssa_plot <- function(ssa_out) {
  var <- value <- time <- type <- NULL # satisfying r check
  df <-
    bind_rows(
      data.frame(time = ssa_out$time, ssa_out$state, type = "state", stringsAsFactors = FALSE) %>%
        gather(var, value, -time, -type),
      data.frame(time = ssa_out$time, ssa_out$propensity, type = "propensity", stringsAsFactors = FALSE) %>%
        gather(var, value, -time, -type),
    )
  if (ncol(ssa_out$buffer) > 0) {
    df <-
      data.frame(time = ssa_out$time, ssa_out$buffer, type = "buffer", stringsAsFactors = FALSE) %>%
      gather(var, value, -time, -type) %>%
      bind_rows(df)
  }
  ggplot2::ggplot(df) +
    ggplot2::geom_path(ggplot2::aes_string("time", "value", colour = "var")) +
    ggplot2::facet_wrap(~type, ncol = 1, scales = "free_y")
}
