test_that("Symmetric methods reverse sign when mx1 and mx2 are swapped", {
  # Example mortality schedules
  a <- 0.001
  b <- 0.07
  x <- 0:100
  mx1 <- a * exp(x * b)
  mx2 <- a / 2 * exp(x * b)

  # Methods known to be symmetric
  sym_methods <- c(
    "arriaga_sym",
    "lopez_ruzicka_sym",
    "chandrasekaran_ii",
    "chandrasekaran_iii",
    "stepwise",
    "horiuchi"
  )

  for (m in sym_methods) {
    d1 <- LEdecomp(mx1, mx2, age = x, method = m)
    d2 <- LEdecomp(mx2, mx1, age = x, method = m)

    expect_equal(d1$LEdecomp, -d2$LEdecomp, tolerance = 1e-10,
                 info = paste("Symmetry test failed for method:", m))
  }
})
