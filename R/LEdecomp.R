#' Function for applying different Life-Expectancy decomposition and sensitivity methods
#' @description A variety of exact or asymptotically exact life expectancy decomposition methods are implemented. Also, several life-expectancy decomposition sensitivity methods are implemented to answer how each age will change with an increase/decrease in life expectancy. See the package README and references for details.
#'
#' @param mx1 numeric. Age-structured mortality rates for population 1 (vector, matrix, or data.frame).
#' @param mx2 numeric. Age-structured mortality rates for population 2 (same shape as `mx1`).
#' @param age integer. Lower bound of each age group. If `NULL`, it will be inferred from data (see Details).
#' @param nx integer vector of age intervals (defaults to 1 when missing).
#' @param n_causes integer or `NULL`. If provided with stacked vectors, forces the number of causes (columns).
#' @param sex1 character. `"m"`,`"f"`, or `"t"`, affects a0 treatment.
#' @param sex2 character. `"m"`,`"f"`, or `"t"`, affects a0 treatment.
#' @param method character. One of the methods in `method_registry$method`.
#' @param closeout logical. Close out at top age (`TRUE`) or assume closed final age group (`FALSE`).
#' @param opt logical. For lifetable, numerical, and instantaneous sensitivity-based methods, optimize rate averaging
#'   to eliminate the decomposition residual?
#' @param tol numeric. Tolerance for rate-averaging optimization.
#' @param Num_Intervals integer. For methods that discretize an integral (e.g., Horiuchi).
#' @param symmetrical logical. For stepwise replacement only: average 1 to 2 and 2 to 1?
#' @param direction character. For stepwise replacement: "up", "down", or "both".
#' @param perturb numeric. Small perturbation for numerical derivatives.
#' @param ... optional arguments passed to `numDeriv::grad()` or other internals
#'
#' @return An object of class `"LEdecomp"`:
#' \itemize{
#'   \item `mx1`, `mx2`, `age`, `sex1`, `sex2`, `method`, `closeout`, `opt`, `tol`, `Num_Intervals`,
#'         `symmetrical`, `direction`, `perturb`
#'   \item `sens`: vector/matrix of sensitivities (same dimensions as inputs)
#'   \item `LE1`, `LE2`: life expectancy for `mx1` and `mx2`
#'   \item `LEdecomp`: vector/matrix of contributions (same shape as inputs)
#' }
#'
#' @seealso [LEdecomp::sen_e0_mx_lt()], [LEdecomp::arriaga()], [LEdecomp::arriaga_sym()],
#'   [LEdecomp::sen_arriaga()], [LEdecomp::sen_arriaga_sym()]
#'
#' @references
#' \insertRef{arriaga1984measuring}{LEdecomp}
#' \insertRef{Chandrasekaran1986}{LEdecomp}
#' \insertRef{preston2000demography}{LEdecomp}
#' \insertRef{Ponnapalli2005}{LEdecomp}
#'
#' @importFrom DemoDecomp horiuchi
#' @importFrom DemoDecomp stepwise_replacement
#' @importFrom numDeriv grad
#' @export
LEdecomp <- function(mx1,
                     mx2,
                     age = NULL,
                     nx = NULL,
                     n_causes = NULL,
                     sex1 = "t",
                     sex2 = sex1,
                     method = c("lifetable",
                                "arriaga", "arriaga_sym",
                                "sen_arriaga", "sen_arriaga_sym",
                                "sen_arriaga_inst", "sen_arriaga_inst2",
                                "sen_arriaga_sym_inst", "sen_arriaga_sym_inst2",
                                "chandrasekaran_ii",
                                "sen_chandrasekaran_ii", "sen_chandrasekaran_ii_inst",
                                "sen_chandrasekaran_ii_inst2",
                                "chandrasekaran_iii",
                                "sen_chandrasekaran_iii", "sen_chandrasekaran_iii_inst",
                                "sen_chandrasekaran_iii_inst2",
                                "lopez_ruzicka", "lopez_ruzicka_sym",
                                "sen_lopez_ruzicka", "sen_lopez_ruzicka_sym",
                                "sen_lopez_ruzicka_inst", "sen_lopez_ruzicka_inst2",
                                "horiuchi", "stepwise", "numerical"),
                     closeout = TRUE,
                     opt = TRUE,
                     tol = 1e-10,
                     Num_Intervals = 20,
                     symmetrical = TRUE,
                     direction = "both",
                     perturb = 1e-6,
                     ...) {

  if (is.null(mx1) || is.null(mx2)) {
    stop("Arguments 'mx1' and 'mx2' must be provided.")
  }
  if (length(mx1) != length(mx2)) {
    stop("Arguments 'mx1' and 'mx2' must have the same length (prior to shaping).")
  }


  # --- Normalize all inputs (shape + ages) in one place ---------------------
  norm <- normalize_inputs(mx1 = mx1, mx2 = mx2, age = age, n_causes = n_causes)
  mx1        <- norm$mx1
  mx2        <- norm$mx2
  age        <- norm$age
  nages      <- norm$nages
  n_causes   <- norm$n_causes
  deez_dims  <- norm$deez_dims

  if (is.null(nx)) {
    nx <- default_nx_from_age(age)
  }
  stopifnot(length(nx) == nages)
  # If different sexes requested, do two runs (male vs male; female vs female), then average
  if (sex1 != sex2) {
    d1 <- LEdecomp(mx1 = mx1, mx2 = mx2, age = age, nx = nx,
                   sex1 = sex1, sex2 = sex1,
                   method = method, closeout = closeout, opt = opt, tol = tol,
                   Num_Intervals = Num_Intervals, symmetrical = symmetrical,
                   direction = direction, perturb = perturb, ...)
    d2 <- LEdecomp(mx1 = mx1, mx2 = mx2, age = age, nx = nx,
                   sex1 = sex2, sex2 = sex2,
                   method = method, closeout = closeout, opt = opt, tol = tol,
                   Num_Intervals = Num_Intervals, symmetrical = symmetrical,
                   direction = direction, perturb = perturb, ...)
    decomp <- (d1$LEdecomp + d2$LEdecomp) / 2
    sen    <- (d1$sens + d2$sens) / 2

    LE2 <- mx_to_e0(rowSums(as.matrix(mx2)), age = age, nx = nx, sex = sex1, closeout = closeout)
    LE1 <- mx_to_e0(rowSums(as.matrix(mx1)), age = age, nx = nx, sex = sex1, closeout = closeout)

    out <- list(mx1 = mx1, mx2 = mx2, age = age,
                sex1 = sex1, sex2 = sex2, method = tolower(method),
                func = get_dec_fun(tolower(match_method(method))),
                closeout = closeout, opt = opt, tol = tol,
                Num_Intervals = Num_Intervals, symmetrical = symmetrical,
                direction = direction, perturb = perturb,
                sens = sen, LE1 = LE1, LE2 = LE2, LEdecomp = decomp)
    class(out) <- "LEdecomp"
    return(out)
  }

  # --- Method selection ------------------------------------------------------
  method <- tolower(match_method(method))
  if (!(method %in% .get_registry()$method)) {
    stop("Method '", method, "' not found in method_registry.")
  }
  dec_fun <- get_dec_fun(method)

  # --- General methods (horiuchi, stepwise, numerical) ----------------------
  gen_methods <- .get_registry()$method[.get_registry()$category == "general"]
  if (method %in% gen_methods) {

    mx_to_e0_vec <- function(mx, n_causes, age, nx, sex, closeout, ...) {
      if (!is.null(n_causes)) {
        dim(mx) <- c(length(mx) / n_causes, n_causes)
        mx <- rowSums(mx)
      }
      mx_to_e0(mx, age = age, nx = nx, sex = sex, closeout = closeout)
    }

    if (method == "stepwise") {
      if (!is.null(n_causes)) {
        message("\nNote: stepwise replacement with multiple causes does not arrange all cause-orderings; results may differ from other methods.\n")
      }
      decomp <- DemoDecomp::stepwise_replacement(
        func = mx_to_e0_vec,
        pars1 = c(mx1), pars2 = c(mx2),
        symmetrical = symmetrical, direction = direction,
        n_causes = n_causes, age = age, nx = nx, sex = sex1, closeout = closeout, ...
      )
    } else if (method == "horiuchi") {
      decomp <- DemoDecomp::horiuchi(
        func = mx_to_e0_vec,
        pars1 = c(mx1), pars2 = c(mx2),
        age = age, nx = nx, n_causes = n_causes,
        sex = sex1, N = Num_Intervals, closeout = closeout, ...
      )
    } else if (method == "numerical") {
      # numerical falls into "general" in your registry; routed via horiuchi-style func
      decomp <- DemoDecomp::horiuchi(
        func = mx_to_e0_vec,
        pars1 = c(mx1), pars2 = c(mx2),
        age = age, nx = nx, n_causes = n_causes,
        sex = sex1, N = Num_Intervals, closeout = closeout, ...
      )
    }

    dim(decomp) <- deez_dims
    delta <- rowSums(as.matrix(mx2)) - rowSums(as.matrix(mx1))
    sen   <- rowSums(as.matrix(decomp)) / delta
  }

  # --- Direct methods (arriaga, lopez_ruzicka, chandrasekaran_ii/iii, etc.) --
  dir_methods <- .get_registry()$method[.get_registry()$category == "direct"]
  if (method %in% dir_methods) {
    if (!is.null(n_causes)) {
      mx1_all     <- rowSums(mx1)
      mx2_all     <- rowSums(mx2)
      delta_all   <- mx2_all - mx1_all
      delta_split <- (mx2 - mx1)
      delta_ratio <- delta_split / matrix(delta_all, nrow = nages, ncol = n_causes)

      decomp_all  <- dec_fun(mx1 = mx1_all, mx2 = mx2_all,
                             age = age, nx = nx, sex1 = sex1, sex2 = sex1, closeout = closeout)
      decomp      <- delta_ratio * matrix(decomp_all, nrow = nages, ncol = n_causes)
      sen         <- matrix(decomp_all / delta_all, nrow = nages, ncol = n_causes)

      # stability warnings
      if (any(abs(delta_all) < 1e-6, na.rm = TRUE)) {
        warning("\nAll-cause rate difference < 1e-6 for at least one age; cause partitioning may be unstable.")
      }
      if (any(rowSums(abs(decomp), na.rm = TRUE) / abs(rowSums(decomp, na.rm = TRUE)) > 3, na.rm = TRUE)) {
        warning("\nAbsolute sum of cause contributions > 3 times total age contribution for at least one age; consider sensitivity-based methods.")
      }
    } else {
      decomp <- dec_fun(mx1 = mx1, mx2 = mx2,
                        age = age, nx = nx, sex1 = sex1, sex2 = sex1, closeout = closeout)
      delta  <- mx2 - mx1
      sen    <- decomp / delta
    }
  }

  # --- Sensitivity methods with a single schedule (opt OK) -------------------
  opt_methods <- .get_registry()$method[.get_registry()$category == "opt_ok"]
  if (method %in% opt_methods) {
    if (!is.null(n_causes)) {
      mx1_all <- rowSums(mx1)
      mx2_all <- rowSums(mx2)
      delta   <- mx2 - mx1  # matrix (cause-specific)

      if (opt) {
        sen_all <- sen_min(mx1 = mx1_all, mx2 = mx2_all,
                           age = age, nx = nx, sex1 = sex1, sex2 = sex1,
                           closeout = closeout, sen_fun = dec_fun, tol = tol)
      } else {
        mx_avg  <- (mx1_all + mx2_all) / 2
        sen_all <- dec_fun(mx = mx_avg, age = age, nx = nx, sex = sex1, closeout = closeout)
      }

      sen    <- matrix(sen_all, nrow = nages, ncol = n_causes)
      decomp <- sen * delta
    } else {
      delta <- mx2 - mx1
      if (opt) {
        sen <- sen_min(mx1 = mx1, mx2 = mx2,
                       age = age, nx = nx, sex1 = sex1, sex2 = sex1,
                       closeout = closeout, sen_fun = dec_fun, tol = tol)
      } else {
        mx_avg <- (mx1 + mx2) / 2
        sen    <- dec_fun(mx = mx_avg, age = age, nx = nx, sex = sex1, closeout = closeout)
      }
      decomp <- sen * delta
    }

    # residual check against e0 difference
    Delta <-
      mx_to_e0(rowSums(as.matrix(mx2)), age = age, nx = nx, sex = sex1, closeout = closeout) -
      mx_to_e0(rowSums(as.matrix(mx1)), age = age, nx = nx, sex = sex1, closeout = closeout)
    if (abs(sum(decomp) - Delta) > 0.1) {
      msg <- paste0("\nSensitivity method (", method, ") ",
                    if (opt) "with opt = TRUE" else "at midpoint",
                    " residual = ", round(sum(decomp) - Delta, 3),
                    ". Consider comparing with other methods.\n")
      warning(msg)
    }
  }

  # --- Sensitivity methods using both mx1 and mx2 (direct_sen) ---------------
  ds_methods <- .get_registry()$method[.get_registry()$category == "direct_sen"]
  if (method %in% ds_methods) {
    delta <- mx2 - mx1
    if (!is.null(n_causes)) {
      mx1_all <- rowSums(mx1)
      mx2_all <- rowSums(mx2)
      sen_all <- dec_fun(mx1 = mx1_all, mx2 = mx2_all,
                         age = age, nx = nx, sex1 = sex1, sex2 = sex1, closeout = closeout)
      sen     <- matrix(sen_all, nrow = nages, ncol = n_causes)
    } else {
      sen <- dec_fun(mx1 = mx1, mx2 = mx2,
                     age = age, nx = nx, sex1 = sex1, sex2 = sex1, closeout = closeout)
    }
    decomp <- sen * delta
  }

  # --- Final outputs ---------------------------------------------------------
  LE2 <- mx_to_e0(rowSums(as.matrix(mx2)), age = age, nx = nx, sex = sex1, closeout = closeout)
  LE1 <- mx_to_e0(rowSums(as.matrix(mx1)), age = age, nx = nx, sex = sex1, closeout = closeout)

  if (norm$return_as == "stacked_vector") {
    # collapse multi-cause matrices back to stacked vectors (cause-major)
    if (is.matrix(decomp)) decomp <- c(decomp)
    if (is.matrix(sen))    sen    <- c(sen)
    if (is.matrix(mx1))    mx1    <- c(mx1)
    if (is.matrix(mx2))    mx2    <- c(mx2)
  } else if (norm$return_as == "matrix") {
    # make sure they are matrices with expected dims
    if (!is.null(norm$deez_dims)) {
      if (!is.null(dim(decomp))) dim(decomp) <- norm$deez_dims
      if (!is.null(dim(sen)))    dim(sen)    <- norm$deez_dims
    }
  } else { # "vector"
    # drop any accidental dims
    decomp <- c(decomp)
    sen    <- c(sen)
  }
  out <- list(mx1 = mx1, mx2 = mx2, age = age,
              sex1 = sex1, sex2 = sex2, method = method, func = dec_fun,
              closeout = closeout, opt = opt, tol = tol,
              Num_Intervals = Num_Intervals, symmetrical = symmetrical,
              direction = direction, perturb = perturb,
              sens = sen, LE1 = LE1, LE2 = LE2, LEdecomp = decomp)
  class(out) <- "LEdecomp"
  out
}

#' @export
print.LEdecomp <- function(x, ...) {
  if (is.null(x) || !"LEdecomp" %in% class(x)) {
    stop("The object is not of class 'LEdecomp'.")
  }
  if (!is.list(x)) stop("The object is not a list. Use 'LEdecomp()' first.")

  is_mat <- is.matrix(x$LEdecomp)
  m <- x$method

  if (!is_mat) {
    if (m %in% c("arriaga", "arriaga_sym", "chandrasekaran_ii",
                 "chandrasekaran_iii", "lopez_ruzicka",
                 "lopez_ruzicka_sym", "horiuchi", "stepwise", "numerical")) {
      cat(paste("Estimated the", m, "Life-Expectancy decomposition method."))
    } else {
      cat("Estimated a sensitivity Life-Expectancy decomposition method.")
    }
  } else {
    if (m %in% c("arriaga", "arriaga_sym", "chandrasekaran_ii", "lopez_ruzicka",
                 "lopez_ruzicka_sym", "horiuchi", "stepwise")) {
      cat(paste("Estimated the", m, "cause-of-death Life-Expectancy decomposition method."))
    } else {
      cat("Estimated the cause-of-death sensitivity Life-Expectancy decomposition method.")
    }
  }
  cat(paste("\nThe total difference explained is:", round(sum(x$LEdecomp), 4)))
}

# -----------------------------------------------------------------------------
# Helpers (no documentation needed; kept lightweight and package-scope)
# -----------------------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b
is_num <- function(x) is.numeric(x) && all(is.finite(x))

parse_age_labels <- function(lbl) {
  if (is.null(lbl)) return(NULL)
  lbl <- as.character(lbl)
  x <- suppressWarnings(as.integer(sub(".*?(\\d+).*", "\\1", lbl, perl = TRUE)))
  if (anyNA(x)) return(NULL)
  x
}

looks_abridged <- function(age) {
  if (length(age) < 3L) return(FALSE)
  age <- as.integer(age)
  length(age) <= 25L && age[1] == 0L && age[2] == 1L && all(diff(age[-(1:2)]) == 5L)
}

strip_age_col <- function(df) {
  if (!is.data.frame(df)) return(list(x = df, age = NULL))
  cn <- tolower(names(df))
  hit <- which(cn == "age")
  if (length(hit) == 1L) {
    age_col <- df[[hit]]
    df <- df[-hit]
    return(list(x = df, age = age_col))
  }
  list(x = df, age = NULL)
}

df_to_matrix <- function(x) {
  if (is.data.frame(x)) return(data.matrix(x))
  x
}

# Abridged grids commonly used (<= 25 age groups)
.abrg_grid <- function(oag) c(0L, 1L, seq.int(5L, oag, by = 5L))
.ABRIDGED_OAG <- c(80L, 85L, 90L, 95L, 100L, 105L, 110L)
.ABRIDGED_CANDIDATES <- lapply(.ABRIDGED_OAG, .abrg_grid)
.ABRIDGED_CANDIDATE_LENGTHS <- vapply(.ABRIDGED_CANDIDATES, length, 1L)

# Inference of the age vector from whatever labels or shapes we can find
infer_age <- function(mx_like, explicit_age = NULL, ripped_age = NULL, prefer_abridged = TRUE) {
  # ----- Case 1: explicit age provided --------------------------------------
  if (!is.null(explicit_age)) return(as.numeric(explicit_age))

  # ----- Case 2: 'age' column ripped from a data.frame ----------------------
  if (!is.null(ripped_age))   return(as.numeric(ripped_age))

  # ----- Case 3: matrix/data.frame rownames look like ages ------------------
  rn <- tryCatch(rownames(mx_like), error = function(e) NULL)
  a  <- parse_age_labels(rn)
  if (!is.null(a) && is_num(a)) return(as.numeric(a))

  # ----- Case 4: vector names look like ages --------------------------------
  if (!is.matrix(mx_like) && !is.data.frame(mx_like)) {
    a <- parse_age_labels(names(mx_like))
    if (!is.null(a) && is_num(a)) return(as.numeric(a))
  }

  # ----- Case 5: matrix/data.frame with no labels -> use nrow ---------------
  if (is.matrix(mx_like) || is.data.frame(mx_like)) {
    nr <- nrow(mx_like)
    if (!is.null(nr) && nr > 0L) return(seq_len(nr) - 1L)
  }

  # ----- Case 6: unlabeled vector (possibly stacked causes) -----------------
  # Try to deduce an age grid length that divides the total length.
  # Priority: (6a) abridged candidates (if prefer_abridged), else (6b) single-year.
  nmx <- length(mx_like)

  # (6a) Abridged candidates: use the longest candidate that divides nmx
  if (prefer_abridged) {
    idx <- which((nmx %% .ABRIDGED_CANDIDATE_LENGTHS) == 0L)
    if (length(idx) > 0L) {
      nages_guess <- max(.ABRIDGED_CANDIDATE_LENGTHS[idx])
      return(as.numeric(seq_len(nages_guess) - 1L))
    }
  }

  # (6b) Single-year candidates: typical max ages (0:max_age) -> lengths are max_age + 1
  # Adjust the range if your data routinely go higher/lower.
  single_year_lengths <- (60L:130L) + 1L  # lengths 61..131 correspond to 0:60 .. 0:130
  div <- single_year_lengths[nmx %% single_year_lengths == 0L]
  if (length(div) > 0L) {
    nages_guess <- max(div)  # prefer the finest grid (e.g., 101 over 81)
    return(as.numeric(seq_len(nages_guess) - 1L))
  }

  # ----- Case 7: final fallback (rare): assume single-year 0..(n-1) ---------
  seq_len(length(mx_like)) - 1L
}

normalize_shape <- function(x, age, n_causes = NULL, what = "mx") {
  sa <- strip_age_col(x); x <- sa$x
  age_override <- sa$age

  x <- df_to_matrix(x)
  if (!is.numeric(x)) stop(sprintf("'%s' must be numeric or coercible to numeric.", what))

  nages <- length(age)

  if (is.matrix(x)) {
    if (nrow(x) != nages) {
      rn_age <- parse_age_labels(rownames(x))
      if (!is.null(rn_age) && length(rn_age) == nrow(x) && is_num(rn_age)) {
        age_override <- as.numeric(rn_age)
        nages <- length(age_override)
      }
    }
    if (nrow(x) != nages) {
      stop(sprintf("For '%s' matrix, nrow(%s) = %d but inferred length(age) = %d.",
                   what, what, nrow(x), nages))
    }
    return(list(x = x[, , drop = FALSE], age_override = age_override))
  }

  # vector (possibly stacked)
  nmx <- length(x)
  if (is.null(n_causes)) {
    if (nmx %% nages != 0L) {
      stop(sprintf("Length of '%s' (%d) is not divisible by length(age) (%d). Supply 'age' and/or 'n_causes'.",
                   what, nmx, nages))
    }
    nc <- as.integer(nmx / nages)
    if (nc == 1L) {
      list(x = as.numeric(x), age_override = age_override)
    } else {
      list(x = matrix(x, nrow = nages, ncol = nc), age_override = age_override)
    }
  } else {
    stopifnot(n_causes >= 1L)
    expected <- nages * n_causes
    if (nmx != expected) {
      stop(sprintf("Given n_causes = %d and length(age) = %d, expected length(%s) = %d but got %d.",
                   n_causes, nages, what, expected, nmx))
    }
    if (n_causes == 1L) {
      list(x = as.numeric(x), age_override = age_override)
    } else {
      list(x = matrix(x, nrow = nages, ncol = n_causes), age_override = age_override)
    }
  }
}

normalize_inputs <- function(mx1, mx2, age = NULL, n_causes = NULL, prefer_abridged = TRUE) {
  # record original forms (before ripping/ coercing)
  orig_form1 <- if (is.data.frame(mx1) && any(tolower(names(mx1)) == "age")) "df_age" else
    if (is.matrix(mx1)) "matrix" else "vector"
  orig_form2 <- if (is.data.frame(mx2) && any(tolower(names(mx2)) == "age")) "df_age" else
    if (is.matrix(mx2)) "matrix" else "vector"

  # peek, infer, normalize (your existing code) ...
  deez_dims_orig <- if (is.matrix(mx1) || is.data.frame(mx1)) dim(mx1) else length(mx1)
  rip1 <- strip_age_col(mx1); rip2 <- strip_age_col(mx2)
  age_inferred <- infer_age(rip1$x, explicit_age = age, ripped_age = rip1$age %||% rip2$age,
                            prefer_abridged = prefer_abridged)

  n1 <- normalize_shape(rip1$x, age_inferred, n_causes, what = "mx1")
  age_final <- n1$age_override %||% age_inferred
  mx1_n <- n1$x

  n2 <- normalize_shape(rip2$x, age_final, n_causes, what = "mx2")
  age_final <- n2$age_override %||% age_final
  mx2_n <- n2$x

  # align benign 1-col differences (existing)

  if (is.matrix(mx1_n) != is.matrix(mx2_n) ||
      (is.matrix(mx1_n) && !identical(dim(mx1_n), dim(mx2_n)))) {
    stop("After normalization, 'mx1' and 'mx2' have incompatible shapes.")
  }

  nages    <- length(age_final)
  n_causes_out <- if (is.matrix(mx1_n)) if (ncol(mx1_n) > 1L) ncol(mx1_n) else NULL else NULL
  deez_dims <- if (is.matrix(mx1_n)) dim(mx1_n) else length(mx1_n)

  # >>> decide how to shape the *output* <<<
  # - If multi-cause and BOTH inputs were vectors originally -> return stacked vector
  # - If any input was matrix or data.frame-with-age -> return matrix
  # - If single-cause -> return vector
  return_as <-
    if (!is.null(n_causes_out)) {
      if (orig_form1 == "vector" && orig_form2 == "vector") "stacked_vector" else "matrix"
    } else {
      "vector"
    }

  list(mx1 = mx1_n, mx2 = mx2_n, age = as.numeric(age_final),
       nages = nages, n_causes = n_causes_out,
       deez_dims = deez_dims, deez_dims_orig = deez_dims_orig,
       return_as = return_as)
}


# tiny helper to accept a vector default for 'method' but force a single match
match_method <- function(method) {
  # method can be a vector default from the signature; make it a single, valid choice
  match.arg(tolower(method), choices = tolower(.get_registry()$method))
}

#' @importFrom utils tail
default_nx_from_age <- function(age) {
  age <- as.numeric(age)
  if (length(age) <= 1L) return(1)              # degenerate case
  # infer widths from the grid; repeat the last width for the terminal group
  dx <- diff(age)
  # common abridged: age == c(0,1,5,10,15,...) -> dx starts 1,4,5,5,...
  # single-year: dx all 1
  c(dx, tail(dx, 1L))
}


# plausible single-year grid lengths to try when x is an unlabeled stacked vector
.single_year_lengths <- function(min_age = 60L, max_age = 130L) {
  # lengths are (max_age + 1) when ages are 0:max_age
  seq.int(min_age + 1L, max_age + 1L)
}

# try to guess "nages" for an unlabeled vector by divisibility
guess_nages_from_length <- function(nxm,
                                    abridged_lengths = .ABRIDGED_CANDIDATE_LENGTHS,
                                    single_year_lengths = .single_year_lengths()) {
  cand <- c(abridged_lengths, single_year_lengths)
  cand <- sort(unique(cand))
  div <- cand[nxm %% cand == 0L]
  if (!length(div)) return(NA_integer_)
  # prefer the largest (finest age grid) that divides nxm
  max(div)
}
