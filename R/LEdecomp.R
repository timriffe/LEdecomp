
#' Function for applying different Life-Expectancy decomposition and sensitivity methods
#' @description A variety of exact or asymptotically exact life expectancy decomposition methods are implemented. Also, several life-expectancy descomposition sensitivity methods are implemented to answer how each age will change with an increase/decrease in life expectancy.
#' These include the lifetable response experiment, using three potential sensitivity functions (lifetable-based `ltre_lt`, an instantaneous Arriaga-derived sensitivity `ltre_arriaga_instantaneous`, or numerical sensitivity `ltre_numerical`), the directional Arriaga method `arriaga`, a symmetrical Arriaga method `arriaga_symmetrical`.
#' These all give similar results, except for directional `arriaga`, which is the most different.
#'
#' @param mx1 numeric. age-structured mortality rates for population 1.
#' @param mx2 numeric. age-structured mortality rates for population 2, must be same length as `mx1`
#' @param age integer. lower bound of each age group, must be same length as `mx1`
#' @param nx integer vector of age intervals, default 1.
#' @param sex1 character. `"m"`,`"f"`, or `"t"`, affects a0 treatment.
#' @param sex2 character. `"m"`,`"f"`, or `"t"`, affects a0 treatment.
#' @param method character. `"lifetable"`, `"arriaga"`, `"arriaga_sym"`, `"sen_arriaga"'`, `"sen_arriaga_sym"`, `"sen_arriaga_inst"'`, `"sen_arriaga_inst2"'`,`"sen_arriaga_sym_inst"`,`"sen_arriaga_sym_inst2"`, `"chandrasekaran_ii"`, `"sen_chandrasekaran_ii"`, `"sen_chandrasekaran_ii_inst"`, `"sen_chandrasekaran_ii_inst2"`, `"chandrasekaran_iii"`, `"sen_chandrasekaran_iii"`, `"sen_chandrasekaran_iii_inst"`, `"sen_chandrasekaran_iii_inst2"`, `"lopez_ruzicka"`, `"lopez_ruzicka_sym"`, `"sen_lopez_ruzicka"`, `"sen_lopez_ruzicka_sym"`, `"sen_lopez_ruzicka_inst"`, `"sen_lopez_ruzicka_inst2"`, `"sen_lopez_ruzicka_sym_inst"`, `"sen_lopez_ruzicka_sym_inst2"`, "`horiuchi"`, `"stepwise"`, `"numerical"`
#' @param closeout logical. Do we handle closeout, or truncate at top age.
#' @param opt logical, default `TRUE`. For lifetable, numerical, and instantaneous sensitivity-based decomposition, shall we optimize the rate averaging to eliminate the decomposition residual?
#' @param tol numeric, default `1e-10`. Tolerance parameter for rate averaging optimization.
#' @param Num_Intervals numeric, default `FALSE`. The number of intervals to integrate over, which is used only for `method` = `"arriaga"` or the sensitvity methods of Arriaga.
#' @param symmetrical logical value, default `FALSE`. This parameter is used for `method` = `"stepwise"` and shall we want to average the results of replacing 1 with 2 and 2 with 1?
#' @param direction character, default `'up'`. This parameter is used for `method` = `"stepwise"`, and corresponds to one of \code{"up"}, \code{"down"}, or \code{"both"}.
#' @param perturb numeric, default `1e-6`, perturbation value that is a numeric constant, and a very small number.
#' @param \dots optional arguments passed to `numDeriv::grad()`
#'
#' @return A list with class \code{"LEdecomp"} including different component in the Life Expectancy decomposition method:
#' * `mx1` numerical. Age-structured mortality rates for population 1.
#' * `mx2` numerical. Age-structured mortality rates for population 2.
#' * `age` integer. Lower bound of each age group, provided by the user. Single ages only.
#' * `sex1` character. One of `"m"` (male) , `"f"` (female),`"t"` (total). Sex of first population. This affects only infant mortality.
#' * `sex2` character. One of `"m"`, `"f"`,`"t"`. Sex of second population. This affects only infant mortality.
#' * `method` character. Method selected by the user to apply the life-expectancy decomposition or sensitivity analysis.
#' * `closeout` logical. Shall we close out ax at the last age using inverse of mx (default TRUE), or assume a closed final age group?
#' * `opt` logical. For certain sensitivity methods, shall we optimize the weighting of mortality rates at which to evaluate the sensitivity? Default TRUE
#' * `tol` numerical. level of tolerance when optimizing sensitivity methods. Default `1e-10`.
#' * `Num_Intervals` integer. For horiuchi method only. The number of intervals to discretize the method over.
#' * `symmetrical` logical. for stepwise replacement algorithm only, do we average results of decomposing 1 vs 2 and 2 vs 1?
#' * `direction` character. For stepwise replacement algorithm only, up age, down age, or the average of both?
#' * `perturb` numerical. Perturbation value that is a numeric constant, and a very small number.
#' * `sens` a vector/matrix of the sensitivity life expectancy decomposition effects that is organized in the same way as \code{mx1} and \code{mx2}.
#' * `LE1` numerical. LE calculated from `mx1`.
#' * `LE1` numerical. LE calculated from `mx2`.
#' * `LEdecomp` numerical. a vector/matrix of the life expectancy decomposition effects that is organized in the same way as \code{mx1} and \code{mx2}.
#'
#' @seealso [LEdecomp::sen_e0_mx_lt()],[LEdecomp::arriaga()],[LEdecomp::arriaga_sym()],
#' [LEdecomp::sen_arriaga()], [LEdecomp::sen_arriaga_sym()]
#'
#' @references
#' \insertRef{arriaga1984measuring}{LEdecomp}
#' \insertRef{Chandrasekaran1986}{LEdecomp}
#' \insertRef{preston2000demography}{LEdecomp}
#' \insertRef{Ponnapalli2005}{LEdecomp}
#'
#' @importFrom DemoDecomp horiuchi
#' @importFrom DemoDecomp horiuchi
#' @importFrom DemoDecomp stepwise_replacement
#' @importFrom numDeriv grad
#'
#' @examples
#' a <- .001
#' b <- .07
#' x <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#'
#' LEdecomp(mx1,mx2,age=x,sex1='t',method = "sen_arriaga_inst")
#'
#' # what methods are available?
#' available_methods()
#' @export
LEdecomp <- function(mx1,
                     mx2,
                     age = (1:length(mx1))-1,
                     nx = rep(1, length(mx1)),
                     sex1 = 't',
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
                     perturb = 1e-6, ...){

  stopifnot(is.vector(mx1) | is.matrix(mx1))
  stopifnot(length(mx1) == length(mx2))

  #check that the main information has been provided
  if(is.null(mx1) || is.null(mx2) || is.null(age)){
    warning("Arguments mx1, mx2, age, and sex1, need to be provided.")
  }
  if (!(method %in% method_registry$method)) {
    stop("Method '", method, "' not found in method_registry.")
  }
  # TR: we need a way to handle case where sex1 != sex2 i.e. a sex decomp,
  # wherein a0 is handled differently for each sex. This will be a mini-
  # recursion solution, i.e. once with male vs male, again with female vs female,
  # then take the average. This is not optimal. What we'd really want is a way to
  # include the Andreev-Kingkade parameters that are used, get their contributions,
  # then add these into the m0 component. But that change would be knarly. An idea
  # to simplify: let le function take a0 arg. Before calling le function, create a0,
  # then decompose including this as parameter, and sum its contribution to m0.
  if (sex1 != sex2){
          # temporary solution for sex1 != sex2
          d1 <- LEdecomp(mx1 = mx1,
                         mx2 = mx2,
                         age = age,
                         nx = nx,
                         sex1 = sex1,
                         sex2 = sex1,
                         method = method,
                         closeout = closeout,
                         opt = opt,
                         tol = tol,
                         Num_Intervals = Num_Intervals,
                         symmetrical = symmetrical,
                         direction = direction,
                         perturb = perturb,
                         ...)
          d2 <- LEdecomp(mx1 = mx1,
                         mx2 = mx2,
                         age = age,
                         nx = nx,
                         sex1 = sex2,
                         sex2 = sex2,
                         method = method,
                         closeout = closeout,
                         opt = opt,
                         tol = tol,
                         Num_Intervals = Num_Intervals,
                         symmetrical = symmetrical,
                         direction = direction,
                         perturb = perturb,
                         ...)$LEdecomp
          decomp <- (d1$LEdecomp + d2$LEdecomp) / 2
          sen    <- (d1$sens + d2$sens) / 2
          #
          LE2 <- mx_to_e0(rowSums(as.matrix(mx2)),
                          age = age,
                          nx = nx,
                          sex = sex1,
                          closeout = closeout)
          LE1 <- mx_to_e0(rowSums(as.matrix(mx1)),
                          age = age,
                          nx = nx,
                          sex = sex1,
                          closeout = closeout)
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
                      "LEdecomp" = decomp)
          class(out) <- "LEdecomp"
          return(out)


  }

  # handle dimensions

  # incoming dimensions
  deez_dims_orig <- dim(mx1)
  # probably we should return decomp in same dims?

  nages     <- length(age)
  # recall if mx1 is matrix length() still works as if dimensionless vector
  nmx       <- length(mx1)
  if (nages < nmx){
    n_causes <- nmx / nages
    n_causes <- round(n_causes)
    stopifnot((n_causes * nages) == nmx)
    dim(mx1) <- c(nages,n_causes)
    dim(mx2) <- c(nages,n_causes)
  } else {
    n_causes <- NULL
  }

  deez_dims <- dim(mx1)

  # sort out method; defaults to first option: lifetable sensitivity method.
  # add optimization option.
  method <- tolower(method)
  method <- match.arg(method,
                      choices =
                        method_registry$method,
                      several.ok = FALSE)
  # assign generic decomp function name
  dec_fun <- get_dec_fun(method)
  # ----------------------------------------------------------------- #
  # method block for stepwise, numerical, horiuchi                    #
  # should handle vector and matrix cases                             #
  # ----------------------------------------------------------------- #
  gen_methods <- method_registry$method[method_registry$category == "general"]
  if (method %in% gen_methods){

    mx_to_e0_vec <- function(mx,n_causes, age, nx, sex, closeout,...){
      if (!is.null(n_causes)){
        dim(mx) <- c(length(mx) / n_causes, n_causes)
        mx <- rowSums(mx)
      }
      mx_to_e0(mx, age = age, nx = nx, sex = sex, closeout = closeout)
    }
    if (method == "stepwise"){
      if (!is.null(n_causes)){
        message("\nFor the case of multiple causes of death and the stepwise replacement algorithm, please note we don't arrange all possible cause-orderings for the decomposition. This method may therefore give results inconsistent with other methods.\n")
      }
      decomp <- DemoDecomp::stepwise_replacement(
                                       func = mx_to_e0_vec,
                                       pars1 = c(mx1),
                                       pars2 = c(mx2),
                                       symmetrical = symmetrical,
                                       direction = direction,
                                       n_causes = n_causes,
                                       age = age,
                                       nx = nx,
                                       sex = sex1, # see current sex solution above
                                       closeout = closeout,
                                       ...)
    }
    if (method == "horiuchi"){
      decomp <- DemoDecomp::horiuchi(func = mx_to_e0_vec,
                                     pars1 = c(mx1), pars2 = c(mx2),
                                     age = age,
                                     nx = nx,
                                     n_causes = n_causes,
                                     sex = sex1, # see current sex solution above
                                     N = Num_Intervals,
                                     closeout = closeout,
                                     ...)
    }
    dim(decomp) <- deez_dims
    delta <- rowSums(as.matrix(mx2)) - rowSums(as.matrix(mx1))
    sen   <- rowSums(as.matrix(decomp)) / delta
  }
  # ----------------------------------------------------------------- #
  # handle all "direct" methods together                              #
  # ----------------------------------------------------------------- #
  dir_methods <- method_registry$method[method_registry$category == "direct"]
  if (method %in% dir_methods){
    # handle causes per Preston Box 4.3 in this case
    if (!is.null(n_causes)){
      mx1_all     <- rowSums(mx1)
      mx2_all     <- rowSums(mx2)
      delta       <- mx2_all - mx1_all
      delta_split <- (mx2 - mx1) / delta
      decomp_all  <- dec_fun(mx1 = mx1_all,
                             mx2 = mx2_all,
                             age = age,
                             nx = nx,
                             sex1 = sex1,
                             sex2 = sex1, # see diff sex solution at start
                             closeout = closeout)
      decomp      <- delta_split * decomp_all
      sen         <- decomp_all / delta
      # Warn if we note unstable ratio (denom close to 0)
      if (any(abs(delta)<1e6)){
        warning("\nPlease check results: You have at least one all-cause rate difference < 1e-6, which can make cause partitioning of results unstable. You should compare using a sensitivity-based method or similar (lifetable, sen_*, horiuchi, stepwise\n")
      }
      # Warn if causes appear to be explaining > 3x the total age effect
      if (any(rowSums(abs(decomp))/abs(decomp) > 3)){
        warning("\nPlease check results: You have at least one age where the absolute contributions from causes are > 3x the total age contribution. You should compare using a sensitivity-based method or similar (lifetable, sen_*, horiuchi, stepwise\n")
      }

    } else {
      decomp  <- dec_fun(mx1 = mx1,
                         mx2 = mx2,
                         age = age,
                         nx = nx,
                         sex1 = sex1,
                         sex2 = sex1, # see diff sex solution at start
                         closeout = closeout)
      delta <- mx2 - mx1
      sen   <- decomp / delta
    }

  }
  # ----------------------------------------------------------------- #
  # sensitivity methods requiring an avg mx; these are instantaneous  #
  # ones plus the lifetable method;                                   #
  # TR: really all of these could use DemoDecomp::ltre(), thereby     #
  # using Num_intervals argument                                      #
  # ----------------------------------------------------------------- #
  opt_methods <- method_registry$method[method_registry$category == "opt_ok"]
  if (method %in% opt_methods){
    if (!is.null(n_causes)){
      # Collapse to all-cause mortality
      mx1_all <- rowSums(mx1)
      mx2_all <- rowSums(mx2)
      delta   <- mx2 - mx1

      # Compute sensitivities on all-cause mx
      if (opt) {
        sen_all <- sen_min(mx1 = mx1_all,
                           mx2 = mx2_all,
                           age = age,
                           nx = nx,
                           sex1 = sex1,
                           sex2 = sex1,
                           closeout = closeout,
                           sen_fun = dec_fun,
                           tol = tol)
      } else {
        mx_avg <- (mx1_all + mx2_all) / 2
        sen_all <- dec_fun(mx = mx_avg,
                           age = age,
                           nx = nx,
                           sex = sex1,
                           closeout = closeout)
      }

      # Expand sensitivity vector to matrix (recycle over causes)
      sen <- matrix(sen_all, nrow = nages, ncol = n_causes)

      # Final decomposition matrix
      decomp <- sen * delta

    } else {
      # All-cause case (no causes of death)
      delta <- mx2 - mx1

      if (opt) {
        sen <- sen_min(mx1 = mx1,
                       mx2 = mx2,
                       age = age,
                       nx = nx,
                       sex1 = sex1,
                       sex2 = sex1,
                       closeout = closeout,
                       sen_fun = dec_fun,
                       tol = tol)
      } else {
        mx_avg <- (mx1 + mx2) / 2
        sen <- dec_fun(mx = mx_avg,
                       age = age,
                       nx = nx,
                       sex = sex1,
                       closeout = closeout)
      }

      decomp <- sen * delta
    }

    # Residual check
    Delta <-
      mx_to_e0(rowSums(as.matrix(mx2)),
               age = age,
               nx = nx,
               sex = sex1,
               closeout = closeout) -
      mx_to_e0(rowSums(as.matrix(mx1)),
               age = age,
               nx = nx,
               sex = sex1,
               closeout = closeout)
    if (abs(sum(decomp) - Delta) > .1){
      msg <- paste0("\nYou used a sensitivity-based method (", method, ") ",
                    if (opt) "with opt = TRUE" else "evaluated at midpoint",
                    ", but the residual is ", round(sum(decomp) - Delta, 3),
                    ". Consider comparing with other methods.\n")
      warning(msg)
    }
  }

  # ----------------------------------------------------------------- #
  # sensitivity methods using both mx1 and mx2,                       #
  # basically all are conversions of direct methods                   #
  # ----------------------------------------------------------------- #
  ds_methods <- method_registry$method[method_registry$category == "direct_sen"]
  if (method %in% ds_methods){
    # this might have dims
    delta  <- mx2 - mx1
    # first block, causes of death
    if (!is.null(n_causes)){
      mx1_all     <- rowSums(mx1)
      mx2_all     <- rowSums(mx2)
      # For these ones, we don't have a choice to optimize; that only works
      # for sensitivity functions where we specify a single rate schedule
      # based on pre-averaging. Perhaps there is a way to imagine these
      # functions with averaging in mind, but not currently set up that way
      sen         <- dec_fun(mx1 = mx1_all,
                             mx2 = mx2_all,
                             age = age,
                             nx = nx,
                             sex1 = sex1,
                             sex2 = sex1, # see diff sex solution at start
                             closeout = closeout)
    } else {
      # second block, all cause mort was given
      sen         <- dec_fun(mx1 = mx1,
                             mx2 = mx2,
                             age = age,
                             nx = nx,
                             sex1 = sex1,
                             sex2 = sex1, # see diff sex solution at start
                             closeout = closeout)
    }
    # this is either a vector or a matrix
    decomp <- sen * delta
  }

  # presumably we've used all options by now
  LE2 <- mx_to_e0(rowSums(as.matrix(mx2)),
                  age = age,
                  nx = nx,
                  sex = sex1,
                  closeout = closeout)
  LE1 <- mx_to_e0(rowSums(as.matrix(mx1)),
                  age = age,
                  nx = nx,
                  sex = sex1,
                  closeout = closeout)

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
              "LEdecomp" = decomp)
  class(out) <- "LEdecomp"
  out

}

#' @export
print.LEdecomp <- function(x, ...) {
  if(!is.null(x)){
    if(!"LEdecomp" %in% class(x))
      stop("The 'x' does not have the 'decomp' structure of R LEdecomp package.")
  }
  if(!is.list(x)){
    stop("The 'x' is not a list. Use 'LEdecomp' function first.")
  }

  if(!is.matrix(x$LEdecomp)){
    if(x$method %in% c("arriaga", "arriaga_sym", "chandrasekaran_ii",
                       "chandrasekaran_iii", "lopez_ruzicka",
                       "lopez_ruzicka_sym", "horiuchi", "stepwise",
                       "numerical")){

      cat(paste("Estimated the", x$method, "Life-Expectancy decomposition method."))

    } else if(x$method %in% c( "sen_arriaga", "sen_arriaga_sym" , "sen_arriaga_inst" ,
                               "sen_arriaga_inst2" ,  "sen_arriaga_sym_inst", "sen_arriaga_sym_inst2", "sen_chandrasekaran_ii" , "sen_chandrasekaran_ii_inst" , "sen_chandrasekaran_ii_inst2" , "sen_chandrasekaran_iii" , "sen_chandrasekaran_iii_inst" , "sen_chandrasekaran_iii_inst2" , "sen_lopez_ruzicka" , "sen_lopez_ruzicka_sym" , "sen_lopez_ruzicka_inst" , "sen_lopez_ruzicka_sym_inst2")){

      values <- c("sen_arriaga", "sen_arriaga_sym", "sen_arriaga_inst",

                  "sen_arriaga_inst2", "sen_arriaga_sym_inst", "sen_arriaga_sym_inst2","sen_chandrasekaran_ii",
                  "sen_chandrasekaran_ii_inst", "sen_chandrasekaran_ii_inst2",
                  "sen_chandrasekaran_iii", "sen_chandrasekaran_iii_inst",
                  "sen_chandrasekaran_iii_inst2", "sen_lopez_ruzicka",
                  "sen_lopez_ruzicka_sym", "sen_lopez_ruzicka_inst",
                  "sen_lopez_ruzicka_inst2")

      names <- c("arriaga", "arriaga_sym", "arriaga_inst", "arriaga_inst2","arriaga_sym_inst","arriaga_sym_inst2",
                 "chandrasekaran_ii", "chandrasekaran_ii_inst", "chandrasekaran_ii_inst2",
                 "chandrasekaran_iii", "chandrasekaran_iii_inst", "chandrasekaran_iii_inst2",
                 "lopez_ruzicka", "lopez_ruzicka_sym", "lopez_ruzicka_inst", "lopez_ruzicka_sym__inst2")
      dicc <- names[which(x$method == values)[1]]

      cat(paste("Estimated the", dicc, "sensitivity Life-Expectancy decomposition method."))

  }
    } else if(is.matrix(x$LEdecomp)){
      if(x$method %in%c("arriaga", "arriaga_sym",  "chandrasekaran_ii" , "lopez_ruzicka","lopez_ruzicka_sym", "horiuchi","stepwise")){

        cat(paste("Estimated the", x$method, "cause-of-death Life-Expectancy decomposition method."))

      } else if(x$method %in%c( "lifetable","sen_arriaga", "sen_arriaga_sym" , "sen_arriaga_inst" , "sen_arriaga_inst2" , "sen_arriaga_sym_inst" , "sen_arriaga_sym_inst2",
                                "sen_chandrasekaran_ii" , "sen_chandrasekaran_ii_inst" , "sen_chandrasekaran_ii_inst2", "sen_chandrasekaran_iii" , "sen_chandrasekaran_iii_inst" , "sen_chandrasekaran_iii_inst2" , "sen_lopez_ruzicka" , "sen_lopez_ruzicka_sym" , "sen_lopez_ruzicka_inst" , "sen_lopez_ruzicka_sym_inst2","numerical")){

        values <- c("sen_arriaga", "sen_arriaga_sym", "sen_arriaga_inst",
                    "sen_arriaga_inst2", "sen_arriaga_sym_inst" ,
                    "sen_arriaga_sym_inst2", "sen_chandrasekaran_ii",
                    "sen_chandrasekaran_ii_inst","sen_chandrasekaran_ii_inst2",
                    "sen_chandrasekaran_iii", "sen_chandrasekaran_iii_inst",
                    "sen_chandrasekaran_iii_inst2", "sen_lopez_ruzicka",
                    "sen_lopez_ruzicka_sym", "sen_lopez_ruzicka_inst",
                    "sen_lopez_ruzicka_sym__inst2")

        names <- c("arriaga", "arriaga_sym", "arriaga_inst", "arriaga_inst2",
                   "sen_arriaga_sym_inst" ,"sen_arriaga_sym_inst2",
                   "chandrasekaran_ii", "chandrasekaran_ii_inst",
                   "chandrasekaran_ii_inst2", "chandrasekaran_iii",
                   "chandrasekaran_iii_inst", "chandrasekaran_iii_inst2",
                   "lopez_ruzicka", "lopez_ruzicka_sym", "lopez_ruzicka_inst",
                   "lopez_ruzicka_sym_inst2")
        dicc <- names[which(x$method == values)[1]]

        cat(paste("Estimated the", dicc, "cause-of-death sensitivity Life-Expectancy decomposition method.\n"))

      }
  }

  cat(paste("\nThe total difference explained is:", round(sum(x$LEdecomp), 4)))
}












