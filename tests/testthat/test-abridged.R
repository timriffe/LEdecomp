
test_that("Abridged ages are inferred for unlabeled stacked vectors and match matrix baseline", {
  skip_on_cran()
  set.seed(42)

  # Build an abridged grid: 0, 1, 5..110 by 5 (length 24)
  age_ab <- c(0L, 1L, seq.int(5L, 110L, by = 5L))
  n_age  <- length(age_ab)
  # nx inferred from age grid (repeat last width for open age group)
  nx_ab  <- c(diff(age_ab), tail(diff(age_ab), 1L))

  # Synthetic Gompertz rates on abridged lower bounds
  a <- 0.001
  b <- 0.07
  mx1_all <- a * exp(age_ab * b)
  mx2_all <- (a / 2) * exp(age_ab * b)

  # Make 3 causes by random weights per age, normalized to 1
  k <- 3
  w1 <- matrix(runif(n_age * k, 0.9, 1.1), nrow = n_age)
  w2 <- matrix(runif(n_age * k, 0.9, 1.1), nrow = n_age)
  w1 <- w1 / rowSums(w1)
  w2 <- w2 / rowSums(w2)

  mx1_mat <- mx1_all * w1
  mx2_mat <- mx2_all * w2

  # Matrix baseline (explicit age and nx)
  res_mat <- LEdecomp(
    mx1 = mx1_mat, mx2 = mx2_mat,
    age = age_ab, nx = nx_ab,
    sex1 = "t", sex2 = "t",
    method = "sen_arriaga", opt = TRUE
  )

  # Stack to unlabeled vectors (column-major; matches dim<- fill)
  mx1_stack <- as.vector(mx1_mat)
  mx2_stack <- as.vector(mx2_mat)

  # Stacked, unlabeled inputs: age = NULL, nx = NULL
  # Expect: age inferred as abridged grid; nx inferred from that age;
  #         output is stacked vector equal to matrix baseline when c()'d.
  res_vec <- LEdecomp(
    mx1 = mx1_stack, mx2 = mx2_stack,
    age = NULL, nx = NULL,
    sex1 = "t", sex2 = "t",
    method = "sen_arriaga", opt = TRUE
  )

  # 1) age should be the abridged grid
  expect_identical(res_vec$age, as.numeric(age_ab))

  # 2) stacked vector output equals matrix baseline when flattened
  expect_equal(res_vec$LEdecomp, c(res_mat$LEdecomp), tolerance = 1e-8)

  # 3) row sums across causes equal all-cause decomposition on abridged grid
  res_all <- LEdecomp(
    mx1 = mx1_all, mx2 = mx2_all,
    age = age_ab, nx = nx_ab,
    sex1 = "t", sex2 = "t",
    method = "sen_arriaga", opt = TRUE
  )
  expect_equal(rowSums(res_mat$LEdecomp), res_all$LEdecomp, tolerance = 1e-7)
})

test_that("Abridged: stacked vector with repeated explicit age infers n_causes and matches matrix", {
  skip_on_cran()
  set.seed(42)

  age_ab <- c(0L, 1L, seq.int(5L, 110L, by = 5L))
  n_age  <- length(age_ab)
  nx_ab  <- c(diff(age_ab), tail(diff(age_ab), 1L))

  a <- 0.001
  b <- 0.07
  mx1_all <- a * exp(age_ab * b)
  mx2_all <- (a / 2) * exp(age_ab * b)

  k <- 3
  w1 <- matrix(runif(n_age * k, 0.9, 1.1), nrow = n_age); w1 <- w1 / rowSums(w1)
  w2 <- matrix(runif(n_age * k, 0.9, 1.1), nrow = n_age); w2 <- w2 / rowSums(w2)

  mx1_mat <- mx1_all * w1
  mx2_mat <- mx2_all * w2

  res_mat <- LEdecomp(
    mx1 = mx1_mat, mx2 = mx2_mat,
    age = age_ab, nx = nx_ab,
    sex1 = "t", sex2 = "t",
    method = "sen_arriaga", opt = TRUE
  )

  mx1_stack <- as.vector(mx1_mat)
  mx2_stack <- as.vector(mx2_mat)
  rep_age   <- rep(age_ab, k)

  res_vec <- LEdecomp(
    mx1 = mx1_stack, mx2 = mx2_stack,
    age = rep_age, nx = NULL,
    sex1 = "t", sex2 = "t",
    method = "sen_arriaga", opt = TRUE
  )

  expect_identical(res_vec$age, as.numeric(age_ab))
  expect_equal(res_vec$LEdecomp, c(res_mat$LEdecomp), tolerance = 1e-8)

  res_all <- LEdecomp(
    mx1 = mx1_all, mx2 = mx2_all,
    age = age_ab, nx = nx_ab,
    sex1 = "t", sex2 = "t",
    method = "sen_arriaga", opt = TRUE
  )
  expect_equal(rowSums(res_mat$LEdecomp), res_all$LEdecomp, tolerance = 1e-7)
})

test_that("Abridged: unlabeled single-cause vector infers abridged ages and matches explicit", {
  skip_on_cran()

  age_ab <- c(0L, 1L, seq.int(5L, 110L, by = 5L))
  n_age  <- length(age_ab)
  nx_ab  <- c(diff(age_ab), tail(diff(age_ab), 1L))

  a <- 0.001
  b <- 0.07
  mx1 <- a * exp(age_ab * b)
  mx2 <- (a / 2) * exp(age_ab * b)

  # explicit baseline
  res_exp <- LEdecomp(
    mx1 = mx1, mx2 = mx2,
    age = age_ab, nx = nx_ab,
    sex1 = "t", sex2 = "t",
    method = "sen_arriaga", opt = TRUE
  )

  # unlabeled (age = NULL, nx = NULL), length equals abridged grid
  res_inf <- LEdecomp(
    mx1 = mx1, mx2 = mx2,
    age = NULL, nx = NULL,
    sex1 = "t", sex2 = "t",
    method = "sen_arriaga", opt = TRUE
  )

  expect_identical(res_inf$age, as.numeric(age_ab))
  expect_equal(res_inf$LEdecomp, res_exp$LEdecomp, tolerance = 1e-8)
})


test_that("Matrix rownames parsed as ages override mismatched age argument", {
  skip_on_cran()

  a <- 0.001
  b <- 0.07
  age <- 0:100
  n_age <- length(age)
  mx1 <- a * exp(age * b)
  mx2 <- (a/2) * exp(age * b)

  # 3-cause matrices with rownames set to ages
  mx1_mat <- cbind(mx1/2, mx1/3, mx1/6)
  mx2_mat <- cbind(mx2/4, mx2/2, mx2/4)
  rownames(mx1_mat) <- rownames(mx2_mat) <- as.character(age)
  colnames(mx1_mat) <- colnames(mx2_mat) <- c("c1","c2","c3")

  # baseline with correct explicit age
  res_base <- LEdecomp(
    mx1 = mx1_mat, mx2 = mx2_mat,
    age = age, nx = rep(1, n_age),
    sex1 = "t", sex2 = "t",
    method = "sen_arriaga", opt = TRUE
  )

  # call with intentionally wrong age; expect a warning and rownames override
  wrong_age <- age + 10
  expect_warning(
    res_rn <- LEdecomp(
      mx1 = mx1_mat, mx2 = mx2_mat,
      age = wrong_age, nx = rep(1, n_age),
      sex1 = "t", sex2 = "t",
      method = "sen_arriaga", opt = TRUE
    ),
    regexp = "Age argument differs from rownames; using rownames"
  )

  # ages should equal parsed rownames (0:100), not wrong_age
  expect_identical(res_rn$age, as.numeric(age))
  # decomposition identical to baseline
  expect_equal(res_rn$LEdecomp, res_base$LEdecomp, tolerance = 1e-8)
})


test_that("Data.frame with an age column is handled and matches matrix baseline", {
  skip_on_cran()

  a <- 0.001
  b <- 0.07
  age <- 0:100
  n_age <- length(age)
  mx1 <- a * exp(age * b)
  mx2 <- (a/2) * exp(age * b)

  # 3-cause matrices
  mx1_mat <- cbind(mx1/2, mx1/3, mx1/6)
  mx2_mat <- cbind(mx2/4, mx2/2, mx2/4)
  colnames(mx1_mat) <- colnames(mx2_mat) <- c("c1","c2","c3")

  # matrix baseline
  res_mat <- LEdecomp(
    mx1 = mx1_mat, mx2 = mx2_mat,
    age = age, nx = rep(1, n_age),
    sex1 = "t", sex2 = "t",
    method = "sen_arriaga", opt = TRUE
  )

  # build data.frames with an age column up front
  df1 <- data.frame(age = age,
                    c1 = mx1_mat[,1],
                    c2 = mx1_mat[,2],
                    c3 = mx1_mat[,3],
                    check.names = FALSE)
  df2 <- data.frame(age = age,
                    c1 = mx2_mat[,1],
                    c2 = mx2_mat[,2],
                    c3 = mx2_mat[,3],
                    check.names = FALSE)

  # call passing the data.frames; age column should be ripped internally
  res_df <- LEdecomp(
    mx1 = df1, mx2 = df2,
    age = NULL, nx = rep(1, n_age),
    sex1 = "t", sex2 = "t",
    method = "sen_arriaga", opt = TRUE
  )

  expect_identical(res_df$age, as.numeric(age))
  expect_equal(as.matrix(res_df$LEdecomp), res_mat$LEdecomp, tolerance = 1e-8)
})


test_that("Abridged: wide data.frame with age column matches matrix baseline", {
  skip_on_cran()
  set.seed(123)

  age_ab <- c(0L, 1L, seq.int(5L, 110L, by = 5L))
  n_age  <- length(age_ab)
  nx_ab  <- c(diff(age_ab), tail(diff(age_ab), 1L))

  a <- 0.001; b <- 0.07
  mx1_all <- a * exp(age_ab * b)
  mx2_all <- (a / 2) * exp(age_ab * b)

  k <- 3
  w1 <- matrix(runif(n_age * k, 0.9, 1.1), nrow = n_age); w1 <- w1 / rowSums(w1)
  w2 <- matrix(runif(n_age * k, 0.9, 1.1), nrow = n_age); w2 <- w2 / rowSums(w2)

  mx1_mat <- mx1_all * w1
  mx2_mat <- mx2_all * w2
  colnames(mx1_mat) <- colnames(mx2_mat) <- c("c1","c2","c3")

  # matrix baseline
  res_mat <- LEdecomp(
    mx1 = mx1_mat, mx2 = mx2_mat,
    age = age_ab, nx = nx_ab,
    sex1 = "t", sex2 = "t",
    method = "sen_arriaga", opt = TRUE
  )

  # wide df: one row per age, 3 cause columns, plus age column
  df1 <- data.frame(age = age_ab,
                    c1 = mx1_mat[,1],
                    c2 = mx1_mat[,2],
                    c3 = mx1_mat[,3],
                    check.names = FALSE)
  df2 <- data.frame(age = age_ab,
                    c1 = mx2_mat[,1],
                    c2 = mx2_mat[,2],
                    c3 = mx2_mat[,3],
                    check.names = FALSE)

  res_df <- LEdecomp(
    mx1 = df1, mx2 = df2,
    age = NULL, nx = nx_ab,
    sex1 = "t", sex2 = "t",
    method = "sen_arriaga", opt = TRUE
  )

  expect_identical(res_df$age, as.numeric(age_ab))
  expect_true(is.data.frame(res_df$LEdecomp))
  expect_identical(dim(res_df$LEdecomp), dim(res_mat$LEdecomp))
  expect_equal(as.matrix(res_df$LEdecomp), res_mat$LEdecomp, tolerance = 1e-8)
})

test_that("Abridged: tidy/long df (repeated age) works when passed as stacked vectors", {
  skip_on_cran()
  set.seed(123)

  age_ab <- c(0L, 1L, seq.int(5L, 110L, by = 5L))
  n_age  <- length(age_ab)
  nx_ab  <- c(diff(age_ab), tail(diff(age_ab), 1L))

  a <- 0.001; b <- 0.07
  mx1_all <- a * exp(age_ab * b)
  mx2_all <- (a / 2) * exp(age_ab * b)

  k <- 3
  w1 <- matrix(runif(n_age * k, 0.9, 1.1), nrow = n_age); w1 <- w1 / rowSums(w1)
  w2 <- matrix(runif(n_age * k, 0.9, 1.1), nrow = n_age); w2 <- w2 / rowSums(w2)

  mx1_mat <- mx1_all * w1
  mx2_mat <- mx2_all * w2
  colnames(mx1_mat) <- colnames(mx2_mat) <- c("c1","c2","c3")

  # matrix baseline
  res_mat <- LEdecomp(
    mx1 = mx1_mat, mx2 = mx2_mat,
    age = age_ab, nx = nx_ab,
    sex1 = "t", sex2 = "t",
    method = "sen_arriaga", opt = TRUE
  )

  # "tidy" long-like vectors: stacked by cause with repeated age column
  mx1_stack <- as.vector(mx1_mat)
  mx2_stack <- as.vector(mx2_mat)
  age_rep   <- rep(age_ab, k)

  res_vec <- LEdecomp(
    mx1 = mx1_stack, mx2 = mx2_stack,
    age = age_rep, nx = NULL,
    sex1 = "t", sex2 = "t",
    method = "sen_arriaga", opt = TRUE
  )

  expect_identical(res_vec$age, as.numeric(age_ab))
  expect_equal(res_vec$LEdecomp, c(res_mat$LEdecomp), tolerance = 1e-8)
})

