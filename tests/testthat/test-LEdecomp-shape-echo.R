test_that("LEdecomp echoes input shape: vector in -> vector out (single cause)", {
  age <- 0:4
  mx1 <- c(0.10, 0.05, 0.03, 0.02, 0.08)
  mx2 <- mx1 * 1.05

  out <- LEdecomp(mx1 = mx1, mx2 = mx2, age = age, sex1 = "f", sex2 = "f",
                  method = "arriaga", closeout = TRUE)

  expect_true(is.list(out))
  expect_true("LEdecomp" %in% names(out))

  d <- out$LEdecomp
  expect_true(is.atomic(d))
  expect_null(dim(d))
  expect_equal(length(d), length(age))
})

test_that("LEdecomp echoes input shape: stacked vector in -> stacked vector out", {
  age <- 0:4
  n_causes <- 3L

  # stacked vector: age varies fastest within each cause (R column-major)
  mx1_mat <- matrix(
    c(0.10, 0.05, 0.03, 0.02, 0.08,
      0.02, 0.01, 0.01, 0.01, 0.02,
      0.03, 0.02, 0.02, 0.02, 0.03),
    nrow = length(age), ncol = n_causes
  )
  mx2_mat <- mx1_mat * 1.05

  mx1 <- as.vector(mx1_mat)
  mx2 <- as.vector(mx2_mat)

  out <- LEdecomp(mx1 = mx1, mx2 = mx2, age = age, n_causes = n_causes,
                  sex1 = "f", sex2 = "f", method = "arriaga", closeout = TRUE)

  d <- out$LEdecomp
  expect_true(is.atomic(d))
  expect_null(dim(d))
  expect_equal(length(d), length(age) * n_causes)
})

test_that("LEdecomp echoes input shape: wide matrix in -> wide matrix out", {
  age <- 0:4
  n_causes <- 3L

  mx1 <- matrix(
    c(0.10, 0.05, 0.03, 0.02, 0.08,
      0.02, 0.01, 0.01, 0.01, 0.02,
      0.03, 0.02, 0.02, 0.02, 0.03),
    nrow = length(age), ncol = n_causes
  )
  mx2 <- mx1 * 1.05

  out <- LEdecomp(mx1 = mx1, mx2 = mx2, age = age,
                  sex1 = "f", sex2 = "f", method = "arriaga", closeout = TRUE)

  d <- out$LEdecomp
  expect_true(is.matrix(d))
  expect_equal(dim(d), c(length(age), n_causes))
})

test_that("LEdecomp echoes input shape: wide data.frame in -> wide data.frame out", {
  age <- 0:4
  n_causes <- 3L

  mx1_mat <- matrix(
    c(0.10, 0.05, 0.03, 0.02, 0.08,
      0.02, 0.01, 0.01, 0.01, 0.02,
      0.03, 0.02, 0.02, 0.02, 0.03),
    nrow = length(age), ncol = n_causes
  )
  mx2_mat <- mx1_mat * 1.05

  mx1 <- data.frame(mx1_mat)
  mx2 <- data.frame(mx2_mat)

  out <- LEdecomp(mx1 = mx1, mx2 = mx2, age = age,
                  sex1 = "f", sex2 = "f", method = "arriaga", closeout = TRUE)

  d <- out$LEdecomp
  expect_true(is.data.frame(d))
  expect_equal(dim(d), c(length(age), n_causes))
})

test_that("LEdecomp echoes input shape: tibble in -> tibble out (optional)", {
  skip_if_not_installed("tibble")

  age <- 0:4
  n_causes <- 3L

  mx1_mat <- matrix(
    c(0.10, 0.05, 0.03, 0.02, 0.08,
      0.02, 0.01, 0.01, 0.01, 0.02,
      0.03, 0.02, 0.02, 0.02, 0.03),
    nrow = length(age), ncol = n_causes
  )
  mx2_mat <- mx1_mat * 1.05

  mx1 <- tibble::as_tibble(as.data.frame(mx1_mat))
  mx2 <- tibble::as_tibble(as.data.frame(mx2_mat))

  out <- LEdecomp(mx1 = mx1, mx2 = mx2, age = age,
                  sex1 = "f", sex2 = "f", method = "arriaga", closeout = TRUE)

  d <- out$LEdecomp
  expect_true(inherits(d, "tbl_df"))
  expect_equal(dim(d), c(length(age), n_causes))
})

test_that("stacked vector + sex1 != sex2 returns stacked vector (no matrix leak)", {
  age <- 0:100
  n_causes <- 18L

  set.seed(1)
  mx1_mat <- matrix(rexp(length(age) * n_causes, rate = 20000), nrow = length(age))
  mx2_mat <- mx1_mat * runif(1, 0.95, 1.05)

  mx1 <- as.vector(mx1_mat)
  mx2 <- as.vector(mx2_mat)

  out <- LEdecomp(
    mx1 = mx1, mx2 = mx2,
    age = age,
    n_causes = n_causes,
    sex1 = "m", sex2 = "f",
    method = "sen_arriaga_sym_inst",
    opt = TRUE
  )

  expect_null(dim(out$LEdecomp))
  expect_equal(length(out$LEdecomp), length(age) * n_causes)
})
