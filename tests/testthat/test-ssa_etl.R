context("ssa etl")

for (i in seq_len(10)) {
  test_that(paste0("ssa etl produces good results, seed ", i), {
    set.seed(i)

    tau <- runif(1, .001, .2)

    meth <- ssa_etl(tau = tau)
    expect_equal(meth$name, "ETL")
    expect_equal(meth$params, list(tau = tau))

    M <- sample.int(100, 1)
    N <- sample.int(26, 1)
    state <- sample.int(1000, N) %>% set_names(sample(letters, N))
    propensity <- sample.int(1000, M)
    nu <- Matrix::Matrix(
      sample(-5L:5L, M * N, replace = TRUE),
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
    expect_equal(out$dtime, tau)

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
    exp_firings <- propensity * tau
    expect_equivalent(avg_firings, exp_firings, tolerance = 1)
  })
}
