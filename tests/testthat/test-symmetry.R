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
    "andreev_sym",
    "sen_arriaga_sym",
    "sen_lopez_ruzicka_sym",
    "sen_andreev_sym",
    "chandrasekaran_ii",
    "chandrasekaran_iii",
    "stepwise",
    "horiuchi",
    "numerical",
    "lifetable",
    "sen_arriaga_sym_inst",
    "sen_andreev_sym_inst",
    "sen_lopez_ruzicka_sym_inst",
    "sen_arriaga_sym_inst2",
    "sen_andreev_sym_inst2",
    "sen_lopez_ruzicka_sym_inst2"
  )

  for (m in sym_methods) {
    d1 <- LEdecomp(mx1, mx2, age = x, method = m, tol=1e-14, perturb = 1e-6)
    d2 <- LEdecomp(mx2, mx1, age = x, method = m, tol=1e-14, perturb = 1e-6)

    expect_equal(d1$LEdecomp, -d2$LEdecomp, tolerance = 1e-7,
                 info = paste("Symmetry test failed for method:", m))
  }
})
