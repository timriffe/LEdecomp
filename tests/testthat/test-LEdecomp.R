test_that("LEdecomp returns expected structure", {
  x <- 0:100
  mx1 <- 0.001 * exp(x * 0.07)
  mx2 <- 0.0005 * exp(x * 0.07)

  result <- LEdecomp(mx1, mx2, age = x, method = "arriaga")

  expect_s3_class(result, "LEdecomp")
  expect_named(result, c("mx1", "mx2", "age", "nx", "sex1", "sex2", "method",
                         "func", "closeout", "opt", "tol", "Num_Intervals",
                         "symmetrical", "direction", "perturb", "sens",
                         "LE1", "LE2", "LEdecomp","cause_names"))
  expect_type(result$LEdecomp, "double")
  expect_length(result$LEdecomp, length(mx1))
})


test_that("Decomposition sum matches LE difference", {
  x <- 0:100
  mx1 <- 0.001 * exp(x * 0.07)
  mx2 <- 0.0005 * exp(x * 0.07)

  out <- LEdecomp(mx1, mx2, age = x, method = "sen_lopez_ruzicka_inst2")
  delta <- out$LE2 - out$LE1
  approx_sum <- sum(out$LEdecomp)

  expect_lt(abs(delta - approx_sum), 1e-6)
})

test_that("All methods return output of correct length", {
  a <- 0.001
  b <- 0.07
  age <- 0:100
  mx1 <- a * exp(age * b)
  mx2 <- a / 2 * exp(age * b)

  methods <- method_registry$method
  n <- length(age)

  for (m in methods) {
    result <- LEdecomp(mx1 = mx1,
                       mx2 = mx2,
                       age = age,
                       method = m)
    output_len <- length(result$LEdecomp)
    expect_equal(output_len, n,
                 info = paste("Incorrect output length for method:", m))
  }
})
