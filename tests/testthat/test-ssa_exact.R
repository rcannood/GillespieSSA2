context("ssa exact")

for (i in seq_len(10)) {
  test_that(paste0("ssa exact produces good results, seed ", i), {
    meth <- ssa_exact()
    expect_equal(meth$name, "exact")
    expect_equal(meth$params, list())

    set.seed(i)

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
    expect_equal(sum(out$firings == 1), 1L)
    expect_equal(sum(out$firings == 0), length(propensity) - 1L)

    expect_length(out$dstate, length(state))
    expect_is(out$dstate, "numeric")
    ix <- which(out$firings == 1)
    expect_equal(out$dstate, nu[, ix])
    expect_true(all(out$dstate + state > 0))

    expect_length(out$dtime, 1)
    expect_is(out$dtime, "numeric")
    expect_gt(out$dtime, 0)
    expect_lte(out$dtime, sum(propensity))
  })
}
