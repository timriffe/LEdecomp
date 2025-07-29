
test_that("Instantaneous methods are stable under small log-scale perturbations", {
  a <- 0.001
  b <- 0.07
  x <- 0:100
  mx <- a * exp(x * b)

  inst_methods <- method_registry$method[method_registry$category == "opt_ok"]
  inst_methods <- setdiff(inst_methods, c("lifetable","numerical"))  # exclude non-perturbation method

  perturb_vals <- c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7)

  for (method in inst_methods) {
    totals <- numeric(length(perturb_vals))
    for (i in seq_along(perturb_vals)) {
      result <- LEdecomp(mx1 = mx,
                         mx2 = mx,  # or use opt = FALSE and just pass mx
                         method = method,
                         age = x,
                         opt = FALSE,
                         perturb = perturb_vals[i])
      totals[i] <- sum(result$LEdecomp)
    }

    # Check that totals vary little as perturb shrinks
    diffs <- diff(totals)
    expect_true(all(abs(diffs) < 1e-8),
                info = paste("Instability detected for method:", method))
}
})

