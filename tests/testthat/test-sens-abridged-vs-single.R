test_that("compare sensitivity magnitudes: single-year vs abridged", {
  skip_on_cran()

  # need the registry
  methods <- LEdecomp::available_methods()

  ## 1. single-year setup
  age_sy <- 0:100
  nx_sy  <- rep(1, length(age_sy))
  a  <- 0.0001
  b  <- 0.07
  mx1_sy <- a * exp(age_sy * b)
  mx2_sy <- (a / 1.2) * exp(age_sy * b)  # roughly +1 year LE
  mx_to_e0(mx2_sy,age=age_sy, sex="t", nx = nx_sy)-mx_to_e0(mx1_sy,age=age_sy, sex="t", nx = nx_sy)
  ## 2. abridged setup (same underlying shape)
  age_ab <- c(0L, 1L, seq.int(5L, 100L, by = 5L))
  nx_ab  <- c(diff(age_ab), tail(diff(age_ab), 1L))

  # collapse single-year mx to abridged by taking the value at the lower bound
  # (for a smoother version you could average over the interval)
  mx1_ab <- abridge_mx(mx1_sy,age=age_sy,sex="t")
  mx2_ab <- abridge_mx(mx2_sy,age=age_sy,sex="t")

  out_list <- vector("list", length(methods))

  for (i in seq_along(methods)) {
    m <- methods[i]

    # single-year run
    res_sy <- LEdecomp::LEdecomp(
      mx1 = mx1_sy,
      mx2 = mx2_sy,
      age = age_sy,
      nx  = nx_sy,
      sex1 = "t",
      method = m,
      opt = TRUE
    )

    # abridged run
    res_ab <- suppressWarnings(LEdecomp::LEdecomp(
      mx1 = mx1_ab,
      mx2 = mx2_ab,
      age = age_ab,
      nx  = nx_ab,
      sex1 = "t",
      method = m,
      opt = TRUE
    ))

    # pull sensitivities (always a vector in current LEdecomp)
    sens_sy <- res_sy$sens
    sens_ab <- res_ab$sens

    out_list[[i]] <- data.frame(
      method      = m,
      max_abs_sy  = max(abs(sens_sy), na.rm = TRUE),
      max_abs_ab  = max(abs(sens_ab), na.rm = TRUE),
      ratio_ab_sy = max(abs(sens_ab), na.rm = TRUE) /
        max(abs(sens_sy),  na.rm = TRUE)
    )
  }

  res_tbl <- do.call(rbind, out_list)

  # sanity: no infinities, no NaNs
  expect_true(all(is.finite(res_tbl$max_abs_sy)))
  expect_true(all(is.finite(res_tbl$max_abs_ab)))
  expect_true(all(is.finite(res_tbl$ratio_ab_sy)))

  c("lopez_ruzicka",
  "sen_lopez_ruzicka",
  "arriaga",
  "sen_arriaga")
})


test_that("sensitivity magnitudes are consistent within method blocks", {
  skip_on_cran()

  # rebuild the table just like before -------------------------------
  methods <- LEdecomp::available_methods()

  age_sy <- 0:100
  nx_sy  <- rep(1, length(age_sy))

  a  <- 0.0001
  b  <- 0.07
  mx1_sy <- a * exp(age_sy * b)
  mx2_sy <- (a / 1.2) * exp(age_sy * b)  # about 3-year LE gap

  # abridge via helper (the one you just added)
  mx1_ab <- abridge_mx(mx1_sy, age_sy)
  mx2_ab <- abridge_mx(mx2_sy, age_sy)
  age_ab <- c(0L, 1L, seq.int(5L, max(age_sy), by = 5L))
  nx_ab  <- c(diff(age_ab), tail(diff(age_ab), 1L))

  out_list <- vector("list", length(methods))
options(warn=2)
  for (i in seq_along(methods)) {
    m <- methods[i]

    res_sy <- LEdecomp::LEdecomp(
      mx1 = mx1_sy,
      mx2 = mx2_sy,
      age = age_sy,
      nx  = nx_sy,
      sex1 = "t",
      method = m,
      opt = TRUE
    )

    res_ab <- LEdecomp::LEdecomp(
      mx1 = mx1_ab,
      mx2 = mx2_ab,
      age = age_ab,
      nx  = nx_ab,
      sex1 = "t",
      method = m,
      opt = TRUE
    )

    sens_sy <- res_sy$sens
    sens_ab <- res_ab$sens

    out_list[[i]] <- data.frame(
      method      = m,
      max_abs_sy  = max(abs(sens_sy), na.rm = TRUE),
      max_abs_ab  = max(abs(sens_ab), na.rm = TRUE),
      ratio_ab_sy = max(abs(sens_ab), na.rm = TRUE) /
        max(abs(sens_sy),  na.rm = TRUE)
    )
  }

  res_tbl <- do.call(rbind, out_list)

  # ------------------------------------------------------------------
  # define blocks
  block1 <- c("lopez_ruzicka",
              "sen_lopez_ruzicka",
              "arriaga",
              "sen_arriaga",
              "andreev","sen_andreev")

  # whatever else is in the registry but not in block1
  block2 <- setdiff(res_tbl$method, block1)

  # function to test one block

  # run checks
  cols_to_check <- c("max_abs_sy", "max_abs_ab", "ratio_ab_sy")

  # block 1
  df_b1 <- res_tbl[res_tbl$method %in% block1, , drop = FALSE]
  expect_equal(df_b1$max_abs_sy, rep(89.45379,6), tolerance = 1e-5)
  expect_equal(df_b1$max_abs_ab, rep(410.3238,6), tolerance = 1e-5)

   # block 2
  df_b2 <- res_tbl[res_tbl$method %in% block2, , drop = FALSE]
  expect_true(sd(df_b2$max_abs_sy) / mean(df_b2$max_abs_sy) < 1e-3)
  expect_true(sd(df_b2$max_abs_ab) / mean(df_b2$max_abs_ab) < 1e-3)
})
