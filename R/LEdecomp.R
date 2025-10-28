#' Function for applying different Life-Expectancy decomposition and sensitivity methods
#' @description A variety of exact or asymptotically exact life expectancy decomposition methods are implemented. Also, several life-expectancy decomposition sensitivity methods are implemented to answer how each age will change with an increase/decrease in life expectancy. See the package README and references for details.
#'
#' @param mx1 numeric. Age-structured mortality rates for population 1 (vector, matrix, or data.frame). See Details section for more info.
#' @param mx2 numeric. Age-structured mortality rates for population 2 (same shape as `mx1`).
#' @param age integer. Lower bound of each age group. If `NULL`, it will be inferred from data (see Details).
#' @param nx integer vector of age intervals (defaults to 1 when missing).
#' @param n_causes integer or `NULL`. If provided with stacked vectors, forces the number of causes (columns).
#' @param cause_names optional character vector of length `n_causes` giving labels for causes. Alternatively detected from `colnames(mx1)` in case given as a `matrix` or `data.frame`
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
#' #' @details
#' Input dimensions are flexible to accommodate different coding styles and data layouts:
#'
#' **Accepted forms of `mx1` and `mx2`:**
#'
#' * **Vector:** A single all-cause mortality schedule, one value per age.
#'   In this case `age` must be the same length, or can be omitted and will default to
#'   `0:(length(mx1)-1)`, unless we detect you might be using abridged age groups.
#'
#' * **Matrix:** Rows represent ages, columns represent causes of death.
#'   Row names, if any, and if numeric, are interpreted as ages and override the supplied `age` argument or our inferences. Column names are retained in the output for clarity.
#'
#' * **Data frame (wide):** Same layout as a matrix, with an optional column named `age`.
#'
#' * **Stacked vector:** A long, concatenated vector representing causes stacked on top of
#'   each other (i.e., all ages for cause 1, then all ages for cause 2, and so on).
#'   If you don't specify age, we try to detect age and the number of causes. But please specify age- it could be stacked also, or not! For example, when used inside a tidy pipeline you might write `mutate(LEdecomp(mx1, mx2, age))` where `age` is repeated for each cause, i.e. the code might look the same as if you were dealing with all-cause data. But in that case be careful data are ordered consistently.
#'
#' **Age detection and inference:**
#'
#' * If `age` is supplied explicitly, it is used as given.
#' * If missing, `LEdecomp()` attempts to infer it from (in order):
#'   row names, names of the input vector, a column named `"age"` in a data frame,
#'   or heuristics for single-year (0:100) or abridged (0,1,5,10,…) schedules.
#' * If `age` is repeated (e.g., `c(0:100, 0:100, 0:100)`), the function assumes
#'   a stacked structure and collapses `age` to its unique sorted values.
#'   The number of repetitions becomes `n_causes`.
#'
#' **Return shape:**
#'
#' The output mirrors the input form:
#' * If the inputs were vectors, outputs are vectors.
#' * If inputs were matrices or data frames (wide), outputs are matrices.
#' * If inputs were stacked vectors, outputs are stacked vectors in the same order.
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
#' @examples
#' ## Simple reproducible setup
#' set.seed(123)
#' a <- 0.001
#' b <- 0.07
#'
#' ## 1) Vector (single cause), single-year ages
#' age <- 0:50
#' mx1 <- a * exp(age * b)
#' mx2 <- (a / 2) * exp(age * b)
#'
#' res_vec <- LEdecomp(
#'   mx1 = mx1, mx2 = mx2,
#'   age = age, nx = rep(1, length(age)),
#'   sex1 = "t", method = "sen_arriaga", opt = TRUE
#' )
#' round(sum(res_vec$LEdecomp), 4)
#'
#' ## 2) Matrix (multiple causes): rows = age, cols = causes
#' ##    Build 3 causes with random positive weights per age
#' k <- 3
#' w1 <- matrix(runif(length(age) * k, 0.9, 1.1), nrow = length(age)); w1 <- w1 / rowSums(w1)
#' w2 <- matrix(runif(length(age) * k, 0.9, 1.1), nrow = length(age)); w2 <- w2 / rowSums(w2)
#' mx1_mat <- (mx1) * w1
#' mx2_mat <- (mx2) * w2
#' colnames(mx1_mat) <- colnames(mx2_mat) <- paste0("c", 1:k)
#' rownames(mx1_mat) <- rownames(mx2_mat) <- as.character(age)
#'
#' res_mat <- LEdecomp(
#'   mx1 = mx1_mat, mx2 = mx2_mat,
#'   age = age, nx = rep(1, length(age)),
#'   sex1 = "t", method = "sen_arriaga", opt = TRUE
#' )
#' ## Check: row-summed cause contributions equal all-cause result
#' res_all <- LEdecomp(
#'   mx1 = mx1, mx2 = mx2,
#'   age = age, nx = rep(1, length(age)),
#'   sex1 = "t", method = "sen_arriaga", opt = TRUE
#' )
#' all.equal(rowSums(res_mat$LEdecomp), res_all$LEdecomp, tolerance = 1e-7)
#'
#' ## 3) Data frame (wide): same as matrix but with an 'age' column
#' df1 <- data.frame(age = age, mx1_mat, check.names = FALSE)
#' df2 <- data.frame(age = age, mx2_mat, check.names = FALSE)
#' res_df <- LEdecomp(
#'   mx1 = df1, mx2 = df2,
#'   age = NULL, nx = rep(1, length(age)),
#'   sex1 = "t", method = "sen_arriaga", opt = TRUE
#' )
#' all.equal(res_df$LEdecomp, res_mat$LEdecomp, tolerance = 1e-8)
#'
#' ## 4) Stacked vector (long/concatenated): all ages for cause 1, then cause 2, etc.
#' ##    If 'age' is repeated per cause, LEdecomp infers n_causes and collapses age.
#' mx1_stack <- as.vector(mx1_mat)  # column-major flattening
#' mx2_stack <- as.vector(mx2_mat)
#' age_rep   <- rep(age, k)         # typical tidy pipeline: age repeated per cause
#'
#' res_stack <- LEdecomp(
#'   mx1 = mx1_stack, mx2 = mx2_stack,
#'   age = age_rep, nx = NULL,
#'   sex1 = "t", method = "sen_arriaga", opt = TRUE
#' )
#' ## Output is a stacked vector matching the matrix baseline when flattened
#' all.equal(res_stack$LEdecomp, c(res_mat$LEdecomp), tolerance = 1e-8)
#'
#' ## 5) Abridged ages (0,1,5,10,...,110): inference when labels are missing
#' age_ab <- c(0L, 1L, seq.int(5L, 110L, by = 5L))
#' nx_ab  <- c(diff(age_ab), tail(diff(age_ab), 1L))
#' mx1_ab <- a * exp(age_ab * b)
#' mx2_ab <- (a / 2) * exp(age_ab * b)
#'
#' ## Explicit abridged example
#' res_ab_explicit <- LEdecomp(
#'   mx1 = mx1_ab, mx2 = mx2_ab,
#'   age = age_ab, nx = nx_ab,
#'   sex1 = "t", method = "sen_arriaga", opt = TRUE
#' )
#'
#' ## Unlabeled abridged vector of the same length: age and nx inferred
#' res_ab_infer <- LEdecomp(
#'   mx1 = mx1_ab, mx2 = mx2_ab,
#'   age = NULL, nx = NULL,
#'   sex1 = "t", method = "sen_arriaga", opt = TRUE
#' )
#' all.equal(res_ab_infer$age, as.numeric(age_ab))
#' all.equal(res_ab_infer$LEdecomp, res_ab_explicit$LEdecomp, tolerance = 1e-8)
#'
#' ## 6) Rownames override age when they look like ages
#' ##    Here we give the wrong 'age' but set rownames to "0","1",...,"50".
#' wrong_age <- age + 10
#' mx1_rn <- mx1_mat; mx2_rn <- mx2_mat
#' rownames(mx1_rn) <- rownames(mx2_rn) <- as.character(age)
#' res_rn <- suppressWarnings(LEdecomp(
#'   mx1 = mx1_rn, mx2 = mx2_rn,
#'   age = wrong_age, nx = rep(1, length(age)),
#'   sex1 = "t", method = "sen_arriaga", opt = TRUE
#' ))
#' all.equal(res_rn$age, as.numeric(age))
#'
#' ## 7) List available methods
#' available_methods()
LEdecomp <- function(mx1,
                     mx2,
                     age = NULL,
                     nx = NULL,
                     n_causes = NULL,
                     cause_names = NULL,          # NEW ARG
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

  # ---------------------------------------------------------------------------
  # TOP: cause_names pre-processing (non-intrusive)
  # ---------------------------------------------------------------------------
  # We will collect an optional 'cn_levels_pref' to carry the intended cause order.
  cn_levels_pref <- NULL

  if (!is.null(cause_names)) {
    # If cause_names length equals length(mx1), treat as stacked labels.
    if (length(cause_names) == length(mx1)) {
      if (is.factor(cause_names)) {
        # Respect factor order
        cn_levels_pref <- levels(cause_names)
      } else {
        # Preserve first-appearance order
        cn_levels_pref <- unique(as.character(cause_names))
      }
      # If n_causes not given, infer it from label levels
      if (is.null(n_causes)) n_causes <- length(cn_levels_pref)
    } else {
      # Else, could be length n_causes; we do not yet know n_causes. Will handle later.
      # If it is a factor, remember its levels order, otherwise the given order.
      if (is.factor(cause_names)) {
        cn_levels_pref <- levels(cause_names)
      } else {
        # Keep as-given; will validate length after normalization.
        cn_levels_pref <- as.character(cause_names)
      }
      # do not set n_causes here; we may not know it yet
    }
  }

  # --- Normalize all inputs (shape + ages) in one place ---------------------
  norm <- normalize_inputs(mx1 = mx1, mx2 = mx2, age = age, n_causes = n_causes)
  mx1        <- norm$mx1
  mx2        <- norm$mx2
  age        <- norm$age
  nages      <- norm$nages
  n_causes   <- norm$n_causes
  deez_dims  <- norm$deez_dims

  # Apply preferred cause names to shaped matrices when possible
  # (matrix case or stacked-reshaped inside normalize_inputs)
  if (!is.null(n_causes) && n_causes >= 1L) {
    final_cols <- NULL

    # 1) If user supplied a vector of length n_causes, take it as-is (respect factor order)
    if (!is.null(cause_names) && length(cause_names) == n_causes) {
      if (is.factor(cause_names)) {
        final_cols <- levels(cause_names)
      } else {
        final_cols <- as.character(cause_names)
      }
    }
    # 2) Else, if we had levels from a long vector or a factor, prefer those
    if (is.null(final_cols) && !is.null(cn_levels_pref)) {
      if (length(cn_levels_pref) == n_causes) {
        final_cols <- as.character(cn_levels_pref)
      }
    }
    # 3) Else, try to recover from existing colnames
    if (is.null(final_cols)) {
      if (is.matrix(mx1) && !is.null(colnames(mx1))) {
        final_cols <- colnames(mx1)
      } else if (is.matrix(mx2) && !is.null(colnames(mx2))) {
        final_cols <- colnames(mx2)
      }
    }
    # 4) Else, synthesize
    if (is.null(final_cols)) {
      final_cols <- paste("cause", seq_len(n_causes), sep = "_")
    }

    # Set on matrices if present
    if (is.matrix(mx1)) colnames(mx1) <- final_cols
    if (is.matrix(mx2)) colnames(mx2) <- final_cols

    # Keep for output regardless of shape
    cause_names_out <- final_cols
  } else {
    cause_names_out <- NULL
  }

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
                sens = sen, LE1 = LE1, LE2 = LE2, LEdecomp = decomp,
                cause_names = cause_names_out)  # NEW
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

  # --- Direct methods --------------------------------------------------------
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

  # --- Sensitivity methods (opt_ok) -----------------------------------------
  opt_methods <- .get_registry()$method[.get_registry()$category == "opt_ok"]
  if (method %in% opt_methods) {
    if (!is.null(n_causes)) {
      mx1_all <- rowSums(mx1)
      mx2_all <- rowSums(mx2)
      delta   <- mx2 - mx1

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

  # --- Direct_sen methods ----------------------------------------------------
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

  LE2 <- mx_to_e0(rowSums(as.matrix(mx2)), age = age, nx = nx, sex = sex1, closeout = closeout)
  LE1 <- mx_to_e0(rowSums(as.matrix(mx1)), age = age, nx = nx, sex = sex1, closeout = closeout)

  # ---------------------------------------------------------------------------
  # BOTTOM: echo input shapes and attach cause_names to output
  # ---------------------------------------------------------------------------
  shape_policy <- norm$return_as %||% {
    if (is.matrix(mx1)) "matrix" else if (!is.null(n_causes) && n_causes > 1L) "stacked_vector" else "vector"
  }

  if (identical(shape_policy, "stacked_vector")) {
    if (is.matrix(decomp)) decomp <- c(decomp)
    if (is.matrix(sen))    sen    <- c(sen)
    if (is.matrix(mx1))    mx1    <- c(mx1)
    if (is.matrix(mx2))    mx2    <- c(mx2)
  } else if (identical(shape_policy, "matrix")) {
    if (!is.null(norm$deez_dims)) {
      if (!is.null(dim(decomp))) dim(decomp) <- norm$deez_dims
      if (!is.null(dim(sen)))    dim(sen)    <- norm$deez_dims
    }
  } else {
    decomp <- c(decomp)
    sen    <- c(sen)
  }

  out <- list("mx1" = mx1,
              "mx2" = mx2,
              "age" = age,
              "sex1" = sex1,
              "sex2" = sex2,
              "method" = method,
              "func" = dec_fun,
              "closeout" = closeout,
              "opt" = opt,
              "tol" = tol,
              "Num_Intervals" = Num_Intervals,
              "symmetrical" = symmetrical,
              "direction" = direction,
              "perturb" = perturb,
              "sens" = sen,
              "LE1" = LE1,
              "LE2" = LE2,
              "LEdecomp" = decomp,
              "cause_names" = if (!is.null(n_causes) && n_causes >= 1L) cause_names_out else NULL)
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
  # (6a) abridged candidates: use the longest candidate that divides nmx
  if (prefer_abridged) {
    idx <- which((nmx %% .ABRIDGED_CANDIDATE_LENGTHS) == 0L)
    if (length(idx) > 0L) {
      # pick the one with the largest length, then return the actual grid
      best_len_idx <- idx[which.max(.ABRIDGED_CANDIDATE_LENGTHS[idx])]
      best_grid <- .ABRIDGED_CANDIDATES[[best_len_idx]]
      return(as.numeric(best_grid))
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

normalize_shape <- function(x, age, n_causes, what = "mx") {
  sa <- strip_age_col(x); x <- sa$x
  age_override <- sa$age

  x <- df_to_matrix(x)
  if (!is.numeric(x)) stop(sprintf("'%s' must be numeric or coercible to numeric.", what))

  nages <- length(age)

  if (is.matrix(x)) {
    # Try to parse rownames as ages
    rn_age <- parse_age_labels(rownames(x))
    has_rn_age <- !is.null(rn_age) && length(rn_age) == nrow(x) && is_num(rn_age)

    if (has_rn_age) {
      # If provided age is NULL OR differs in any value, prefer rownames
      age_differs <- is.null(age) ||
        length(age) != length(rn_age) ||
        any(as.numeric(age) != as.numeric(rn_age))
      if (age_differs) {
        warning(sprintf("Age argument differs from rownames; using rownames for '%s'.", what), call. = FALSE)
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

  # ... vector branch unchanged ...
  nmx <- length(x)
  if (is.null(n_causes)) {
    if (nmx %% nages != 0L) {
      stop(sprintf(
        "Length of '%s' (%d) is not divisible by length(age) (%d). Supply 'age' and/or 'n_causes'.",
        what, nmx, nages
      ))
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
      stop(sprintf(
        "Given n_causes = %d and length(age) = %d, expected length(%s) = %d but got %d.",
        n_causes, nages, what, expected, nmx
      ))
    }
    if (n_causes == 1L) {
      list(x = as.numeric(x), age_override = age_override)
    } else {
      list(x = matrix(x, nrow = nages, ncol = n_causes), age_override = age_override)
    }
  }
}

normalize_inputs <- function(mx1,
                             mx2,
                             age = NULL,
                             n_causes = NULL,
                             prefer_abridged = TRUE) {
  # --- record original types/dims for later echoing -------------------------
  deez_dims_orig   <- if (is.matrix(mx1) || is.data.frame(mx1)) dim(mx1) else length(mx1)
  orig_is_matrix1  <- is.matrix(mx1) || is.data.frame(mx1)
  orig_is_matrix2  <- is.matrix(mx2) || is.data.frame(mx2)
  orig_is_vector1  <- is.atomic(mx1) && is.null(dim(mx1))
  orig_is_vector2  <- is.atomic(mx2) && is.null(dim(mx2))

  # --- explicit 'age' may be repeated (e.g., 0:100,0:100,...); harmonize ----
  if (!is.null(age)) {
    # If mx1 is a vector, we can check coherence with its length
    x_len <- if (!orig_is_matrix1) length(mx1) else NULL
    h <- harmonize_repeated_age(age, x_len = x_len)
    age      <- h$age
    n_causes <- n_causes %||% h$n_causes
  }

  # examine age columns before coercion so infer_age can use them ----------
  rip1 <- strip_age_col(mx1)  # list(x = <mx without age col>, age = <age or NULL>)
  rip2 <- strip_age_col(mx2)

  # --- infer age (priority: explicit > age col > rownames > names > nrow > heuristics)
  age_inferred <- infer_age(
    mx_like        = rip1$x,
    explicit_age   = age,
    ripped_age     = rip1$age %||% rip2$age,
    prefer_abridged = prefer_abridged
  )

  # --- normalize shapes; allow age override by df age col / rownames ----------
  n1 <- normalize_shape(rip1$x, age_inferred, n_causes, what = "mx1")
  age_final <- n1$age_override %||% age_inferred
  mx1_n     <- n1$x

  n2 <- normalize_shape(rip2$x, age_final, n_causes, what = "mx2")
  age_final <- n2$age_override %||% age_final
  mx2_n     <- n2$x

  # --- align benign differences: vector vs 1-col matrix ----------------------
  if (is.matrix(mx1_n) && !is.matrix(mx2_n)) {
    if (ncol(mx1_n) == 1L && length(mx2_n) == length(age_final)) {
      mx2_n <- matrix(mx2_n, nrow = length(age_final), ncol = 1L)
    }
  } else if (!is.matrix(mx1_n) && is.matrix(mx2_n)) {
    if (ncol(mx2_n) == 1L && length(mx1_n) == length(age_final)) {
      mx1_n <- matrix(mx1_n, nrow = length(age_final), ncol = 1L)
    }
  }

  # --- final compatibility guard ---------------------------------------------
  if (is.matrix(mx1_n) != is.matrix(mx2_n) ||
      (is.matrix(mx1_n) && !identical(dim(mx1_n), dim(mx2_n)))) {
    stop("After normalization, 'mx1' and 'mx2' have incompatible shapes.", call. = FALSE)
  }

  # --- derived metadata ------------------------------------------------------
  nages         <- length(age_final)
  n_causes_out  <- if (is.matrix(mx1_n)) { if (ncol(mx1_n) > 1L) ncol(mx1_n) else NULL } else NULL
  deez_dims     <- if (is.matrix(mx1_n)) dim(mx1_n) else length(mx1_n)

  # detect if caller passed stacked vectors (vector in  matrix after norm with >1 causes)
  was_stacked_vec <- orig_is_vector1 && is.matrix(mx1_n) && !is.null(n_causes_out) && n_causes_out > 1L

  # decide how outputs should be returned to echo caller's shape
  return_as <- if (orig_is_matrix1 || orig_is_matrix2) {
    "matrix"
  } else if (was_stacked_vec) {
    "stacked_vector"
  } else {
    "vector"
  }

  list(
    mx1 = mx1_n,
    mx2 = mx2_n,
    age = as.numeric(age_final),
    nages = nages,
    n_causes = n_causes_out,
    deez_dims = deez_dims,
    deez_dims_orig = deez_dims_orig,
    was_stacked_vec = was_stacked_vec,
    return_as = return_as
  )
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

# If 'age' is explicitly provided and looks like repeated 0:100,0:100,...,
# derive n_causes from repeat counts and reduce age to unique sorted values.
harmonize_repeated_age <- function(age, x_len = NULL, allow_gaps = FALSE) {
  if (is.null(age)) return(list(age = NULL, n_causes = NULL))

  a <- as.numeric(age)
  u <- sort(unique(a))
  # counts of each unique age value, in order
  counts <- as.integer(tabulate(factor(a, levels = u)))

  # If every age appears exactly once, nothing to do
  if (all(counts == 1L)) {
    return(list(age = u, n_causes = NULL))
  }

  # If repeats exist, they must all be equal (no gaps / irregular merges)
  if (length(unique(counts)) != 1L) {
    msg <- paste0(
      "Explicit 'age' appears to be repeated per cause, but counts differ by age.\n",
      "This suggests missing rows (e.g., an incomplete cross of age x cause).\n",
      "Tip: regularize your data (e.g., tidyr::complete(age, cause)) so each age appears the same number of times."
    )
    stop(msg, call. = FALSE)
  }

  k <- counts[1L]
  if (!is.null(x_len)) {
    # Optional coherence check with the length of the stacked mx vector(s)
    # When age was repeated per cause, its length should equal nages * k.
    nages <- length(u)
    if (length(a) != nages * k) {
      stop(sprintf("Length(age) = %d does not equal length(unique(age)) x n_causes.",
                   length(a), nages, k), call. = FALSE)
    }
    # If x_len is known, check that the stacked vector length matches too
    if (!is.null(x_len) && (x_len %% nages != 0L || (x_len / nages) != k)) {
      warning(sprintf(
        "Detected n_causes = from repeated ages, but length(mx) is not a multiple thereof; please verify.",
        k, x_len, nages, k
      ), call. = FALSE)
    }
  }

  list(age = u, n_causes = k)
}
