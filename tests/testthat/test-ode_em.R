context("ode em")

for (i in seq_len(10)) {
  test_that(paste0("ode em produces good results, seed ", i), {
    set.seed(i)

    tau <- runif(1, .001, .2)
    noise_strength <- rbinom(1, 10, .1)

    meth <- ode_em(tau = tau, noise_strength = noise_strength)
    expect_equal(meth$name, "EM")
    expect_equal(meth$params, list(tau = tau, noise_strength = noise_strength))

    M <- sample.int(100, 1)
    N <- sample.int(26, 1)
    state <- sample.int(1000, N) %>% set_names(sample(letters, N))
    propensity <- rnorm(M, mean = 50, sd = 10) %>% pmax(0)
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
    expect_lte(sum(abs(out$firings - propensity * tau)), .01)

    expect_length(out$dstate, length(state))
    expect_is(out$dstate, "numeric")
    expected_dstate <- (nu %*% out$firings)[,1]
    if (noise_strength > 0) {
      expect_true(any(abs(out$dstate - expected_dstate) > .1))
    } else {
      expect_true(sum(abs(out$dstate - expected_dstate)) < .1)
    }

    expect_length(out$dtime, 1)
    expect_is(out$dtime, "numeric")
    expect_equal(out$dtime, tau)

    if (noise_strength > 0) {
      noise <- sapply(seq_len(1000), function(i) {
        out <- test_ssa_method_step(
          meth,
          state,
          propensity,
          nu
        )
        (out$dstate - expected_dstate) / noise_strength / sqrt(abs(state))
      })
      expect_equal(mean(noise), 0, tolerance = .01)
      expect_equal(sd(noise), tau, tolerance = .01)
    }
  })
}

