
test_that("Arriaga, Chandrasekaran II/III, and Lopez-Ruzicka decompositions give identical results", {
  a <- 0.001
  b <- 0.07
  x <- 0:100
  mx1 <- a * exp(x * b)
  mx2 <- a / 2 * exp(x * b)

  # Compute decompositions
  arr     <- arriaga_sym(mx1, mx2, age = x)
  ch_ii   <- chandrasekaran_II(mx1, mx2, age = x)
  ch_iii  <- chandrasekaran_III(mx1, mx2, age = x)
  lr      <- lopez_ruzicka_sym(mx1, mx2, age = x)

  # Compare pairwise
  expect_equal(arr, ch_ii, tolerance = 1e-10)
  expect_equal(arr, ch_iii, tolerance = 1e-10)
  expect_equal(arr, lr, tolerance = 1e-10)
})

test_that("Symmetrical sensitivity variants give identical results", {
  a <- 0.001
  b <- 0.07
  x <- 0:100
  mx1 <- a * exp(x * b)
  mx2 <- a / 2 * exp(x * b)

  s1 <- sen_arriaga_sym(mx1, mx2, age = x)
  s2 <- sen_chandrasekaran_II(mx1, mx2, age = x)
  s3 <- sen_chandrasekaran_III(mx1, mx2, age = x)
  s4 <- sen_lopez_ruzicka_sym(mx1, mx2, age = x)

  expect_equal(s1, s2, tolerance = 1e-10)
  expect_equal(s1, s3, tolerance = 1e-10)
  expect_equal(s1, s4, tolerance = 1e-10)
})


test_that("Directional Arriaga and Lopez-Ruzicka decompositions give identical results", {
  a <- 0.001
  b <- 0.07
  x <- 0:100
  mx1 <- a * exp(x * b)
  mx2 <- a / 2 * exp(x * b)

  a1 <- arriaga(mx1, mx2, age = x)
  a2 <- lopez_ruzicka(mx1, mx2, age = x)

  expect_equal(a1, a2, tolerance = 1e-10)
})


