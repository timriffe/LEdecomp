test_that("Arriaga and Lopez-Ruzicka decompositions are additive", {
  a <- 0.001
  b <- 0.07
  x <- 0:100
  mx1 <- a * exp(x * b)
  mx2 <- a / 2 * exp(x * b)

  e0_1 <- mx_to_e0(mx1, age = x, sex = "t")
  e0_2 <- mx_to_e0(mx2, age = x, sex = "t")
  delta_e0 <- e0_2 - e0_1

  cc_arriaga <- arriaga(mx1, mx2, age = x, sex1 = "t", sex2 = "t")
  cc_lopez <- lopez_ruzicka(mx1, mx2, age = x, sex1 = "t", sex2 = "t")

  total_arriaga <- sum(cc_arriaga)
  total_lopez <- sum(cc_lopez)

  expect_equal(total_arriaga, delta_e0, tolerance = 1e-8,
               info = paste0("Arriaga method not additive: ",
                             total_arriaga, " vs ", delta_e0))

  expect_equal(total_lopez, delta_e0, tolerance = 1e-8,
               info = paste0("Lopez-Ruzicka method not additive: ",
                             total_lopez, " vs ", delta_e0))
})

test_that("All instantaneous methods (with opt = TRUE) give equal total effects", {
  # Setup example mortality schedule
  a <- 0.001
  b <- 0.07
  x <- 0:100
  mx1 <- a * exp(x * b)
  mx2 <- a / 2 * exp(x * b)

  # List all instantaneous method names
  inst_methods <- .get_registry()$method[.get_registry()$category == "opt_ok"]

  # Calculate reference Δe₀
  e0_1 <- mx_to_e0(mx1, age = x)
  e0_2 <- mx_to_e0(mx2, age = x)
  delta_e0 <- e0_2 - e0_1

  # Loop through each method and test that the sum of contributions ≈ Δe₀
  for (m in inst_methods) {
    result <- LEdecomp(mx1 = mx1, mx2 = mx2, age = x, method = m, opt = TRUE)
    contribution_sum <- sum(result$LEdecomp)
    expect_equal(contribution_sum, delta_e0, tolerance = 1e-8,
                 info = paste("Failed for method:", m))
  }
})

test_that("Sequential additivity holds in total contributions", {
  a <- 0.001
  b <- 0.07
  x <- 0:100
  mx1 <- a * exp(x * b)
  mx2 <- a / 2 * exp(x * b)
  mx_mid <- (mx1 + mx2) / 2

  methods_to_test <- .get_registry()$method[
    .get_registry()$category %in% c("direct", "direct_sen", "opt_ok")
  ]

  for (m in methods_to_test) {
    res_1 <- LEdecomp(mx1, mx_mid, age = x, method = m, opt = TRUE)
    res_2 <- LEdecomp(mx_mid, mx2, age = x, method = m, opt = TRUE)
    res_full <- LEdecomp(mx1, mx2, age = x, method = m, opt = TRUE)

    total_partial <- sum(res_1$LEdecomp) + sum(res_2$LEdecomp)
    total_full <- sum(res_full$LEdecomp)

    expect_equal(
      total_partial, total_full, tolerance = 1e-8,
      info = paste("Sequential additivity failed for method:", m)
    )
  }
})


test_that("Direct sensitivity methods give LEdecomp ≈ sensitivity × (mx2 - mx1)", {
  a <- 0.001
  b <- 0.07
  x <- 0:100
  mx1 <- a * exp(x * b)
  mx2 <- a / 2 * exp(x * b)

  direct_sens_methods <- .get_registry()$method[.get_registry()$category == "direct_sen"]

  for (m in direct_sens_methods) {
    result <- LEdecomp(mx1 = mx1, mx2 = mx2, age = x, method = m)
    expected_contrib <- result$sens * (mx2 - mx1)
    expect_equal(result$LEdecomp, expected_contrib, tolerance = 1e-10,
                 info = paste("Mismatch in LEdecomp vs sens * delta for method:", m))
  }
})
