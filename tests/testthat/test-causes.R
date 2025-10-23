test_that("Matrix input (causes of death) gives consistent total contributions", {
  a <- 0.001
  b <- 0.07
  age <- 0:100
  n_age <- length(age)

  # Create 3 synthetic causes by perturbing mx1/mx2 slightly
  base_mx1 <- a * exp(age * b)
  base_mx2 <- (a / 2) * exp(age * b)

  k <- 3
  set.seed(123)

  # Random positive weights for each cause (per age)
  weights1 <- matrix(runif(n_age * k, 0.9, 1.1), nrow = n_age)
  weights2 <- matrix(runif(n_age * k, 0.9, 1.1), nrow = n_age)

  # Normalize so that rowSums = 1
  weights1 <- weights1 / rowSums(weights1)
  weights2 <- weights2 / rowSums(weights2)

  # Apply weights to base all-cause mx to create cause-specific mx matrices
  mx1_mat <- base_mx1 * weights1
  mx2_mat <- base_mx2 * weights2
  colnames(mx2_mat) <-colnames(mx1_mat) <- LETTERS[1:k]

  # stepwise will get diff results guaranteed
  # horiuchi would need a very high Num_
  methods <- setdiff(.get_registry()$method, c("stepwise","horiuchi"))
  for (m in methods) {
    result_cause <- suppressWarnings(LEdecomp(mx1 = mx1_mat,
                             mx2 = mx2_mat,
                             age = age,
                             nx = rep(1,n_age),
                             method = m))
    result_all <- LEdecomp(mx1 = base_mx1,
                           mx2 = base_mx2,
                           age = age,
                           nx = rep(1,n_age),
                           method = m)

    # Compare total effect by summing cause-specific contributions
    sum_cause <- rowSums(as.matrix(result_cause$LEdecomp))
    total_all <- result_all$LEdecomp
    # 1e-6 is like 3-hour precision
    expect_equal(sum_cause, total_all, tolerance = 1e-6,
                 info = paste("Mismatch for method:", m))
  }
})

# explicitly assert that stepwise gives inconsistent results
test_that("Stepwise method gives different results for cause-specific decomposition", {
  a <- 0.001
  b <- 0.07
  age <- 0:100
  n_age <- length(age)

  # Create 3 synthetic causes by perturbing mx1/mx2 slightly
  base_mx1 <- a * exp(age * b)
  base_mx2 <- (a / 2) * exp(age * b)

  k <- 3
  set.seed(123)

  # Random positive weights for each cause (per age)
  weights1 <- matrix(runif(n_age * k, 0.9, 1.1), nrow = n_age)
  weights2 <- matrix(runif(n_age * k, 0.9, 1.1), nrow = n_age)

  # Normalize so that rowSums = 1
  weights1 <- weights1 / rowSums(weights1)
  weights2 <- weights2 / rowSums(weights2)

  # Apply weights to base all-cause mx to create cause-specific mx matrices
  mx1_mat <- base_mx1 * weights1
  mx2_mat <- base_mx2 * weights2
  colnames(mx2_mat) <-colnames(mx1_mat) <- LETTERS[1:k]

  # Get stepwise and reference (e.g., lifetable) decompositions
  result_stepwise_causes <- suppressMessages(suppressWarnings(
    LEdecomp(mx1 = mx1_mat, mx2 = mx2_mat,
             age = age, nx = rep(1,n_age), method = "stepwise")$LEdecomp
  ))

  result_stepwise_all <- LEdecomp(mx1 = base_mx1,
                               mx2 = base_mx2,
                               age = age,
                               nx = rep(1,n_age),
                               method = "stepwise")$LEdecomp

  # Compare total contributions
  diffs <- rowSums(result_stepwise_causes) - result_stepwise_all

  expect_gt(sum(abs(diffs)), 1e-5,
            label = "Cause-specific and all-cause stepwise decompositions should differ")
})
