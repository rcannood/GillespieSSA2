context("ssa btl")

for (i in seq_len(10)) {
  test_that(paste0("ssa btl produces good results, seed ", i), {
    set.seed(i)
    mean_firings <- rnorm(1, mean = 10, sd = 2) %>% pmax(1)

    meth <- ssa_btl(mean_firings = mean_firings)
    expect_equal(meth$name, "BTL")
    expect_equal(meth$params, list(mean_firings = mean_firings))

    M <- sample.int(100, 1)
    N <- sample.int(26, 1)
    state <- sample.int(100, N) %>% set_names(sample(letters, N))
    propensity <- sample.int(1000, M) %>% pmax(20)
    nu <- Matrix::Matrix(
      rbinom(M * N, 100, .001) * sample(c(-1, +1), M * N, replace = TRUE, prob = c(.05, .95)),
      nrow = N,
      sparse = TRUE
    )

    out <- test_ssa_method_step(
      meth,
      state,
      propensity,
      nu
    )

    expect_length(out$firings, length(propensity))
    expect_is(out$firings, "numeric")
    expect_true(all(out$firings >= 0))

    expect_length(out$dstate, length(state))
    expect_is(out$dstate, "numeric")
    expected_dstate <- (nu %*% out$firings)[,1]
    expect_equal(out$dstate, expected_dstate)

    expect_length(out$dtime, 1)
    expect_is(out$dtime, "numeric")
    expected_tau <- mean_firings / sum(propensity)
    expect_equal(out$dtime, expected_tau)

    firings <- lapply(seq_len(1000), function(i) {
      out <- test_ssa_method_step(
        meth,
        state,
        propensity,
        nu
      )
      out$firings
    })

    avg_firings <- Reduce(`+`, firings) / length(firings)
    exp_firings <- propensity * expected_tau
    expect_equivalent(avg_firings, exp_firings, tolerance = 1)

    expect_true(all(out$dstate + state >= 0))
  })
}
