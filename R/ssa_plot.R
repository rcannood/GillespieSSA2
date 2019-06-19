#' @export
ssa_plot <- function(ssa_out) {
  df <-
    ssa_out$output %>%
    select(time, state) %>%
    mutate(state_name = map(state, names)) %>%
    unnest(state, state_name)
  ggplot(df) +
    geom_path(aes(time, state, colour = state_name))
}
