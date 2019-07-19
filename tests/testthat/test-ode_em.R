context("ode em")


for (i in seq_len(10)) {
  test_that(paste0("ode em produces good results, seed ", i), {
    set.seed(i)

    tau <- runif(1, .001, .2)
    noise_strength <- rbinom(1, 10, .1)

    meth <- ode_em(tau = tau, noise_strength = noise_strength)
    expect_equal(meth$name, "EM")
    expect_equal(meth$params, list(tau = tau, noise_strength = noise_strength))

    ssa_alg <- meth$factory()

    M <- sample.int(100, 1)
    N <- sample.int(26, 1)
    state <- sample.int(1000, N) %>% set_names(sample(letters, N))
    propensity <- rnorm(M, mean = 50, sd = 10) %>% pmax(0)
    nu <- Matrix::Matrix(
      sample(-5L:5L, M * N, replace = TRUE),
      nrow = N,
      sparse = TRUE
    )

    out <- test_ssa_step(
      ssa_alg,
      state,
      propensity,
      nu_i = nu@i,
      nu_p = nu@p,
      nu_x = nu@x
    )
    expect_length(out$firings, length(propensity))
    expect_is(out$firings, "numeric")
    expect_lte(sum(abs(out$firings - propensity * tau)), .01)

    expect_length(out$dstate, length(state))
    expect_is(out$dstate, "numeric")
    expected_dstate <- (nu %*% out$firings)[,1]
    if (noise_strength > 0) {
      expect_true(any(abs(out$dstate - expected_dstate) > 0))
    } else {
      expect_true(sum(abs(out$dstate - expected_dstate)) < .1)
    }

    expect_length(out$dtime, 1)
    expect_is(out$dtime, "numeric")
    expect_equal(out$dtime, tau)

    if (noise_strength > 0) {
      dstates <- lapply(seq_len(1000), function(i) {
        out <- test_ssa_step(
          ssa_alg,
          state,
          propensity,
          nu_i = nu@i,
          nu_p = nu@p,
          nu_x = nu@x
        )
        out$dstate - expected_dstate
      })

      avg_dstates <- Reduce(`+`, dstates) / length(dstates)
      expect_true(all(abs(avg_dstates / sqrt(abs(state)) / noise_strength) - tau < 0))
    }
  })
}

