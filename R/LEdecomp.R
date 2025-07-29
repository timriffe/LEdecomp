
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
#' @param method character. `"lifetable"`, `"arriaga"`, `"arriaga_sym"`, `"sen_arriaga"'`, `"sen_arriaga_sym"`, `"sen_arriaga_inst"'`, `"sen_arriaga_inst2"'`, `"chandrasekaran_ii"`, `"sen_chandrasekaran_ii"`, `"sen_chandrasekaran_ii_inst"`, `"sen_chandrasekaran_ii_inst2"`, `"chandrasekaran_iii"`, `"sen_chandrasekaran_iii"`, `"sen_chandrasekaran_iii_inst"`, `"sen_chandrasekaran_iii_inst2"`, `"lopez_ruzicka"`, `"lopez_ruzicka_sym"`, `"sen_lopez_ruzicka"`, `"sen_lopez_ruzicka_sym"`, `"sen_lopez_ruzicka_inst"`, `"sen_lopez_ruzicka_inst2"`, `"horiuchi"`, `"stepwise"`, `"numerical"`
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
                        c("lifetable", "arriaga", "arriaga_sym",
                          "sen_arriaga", "sen_arriaga_sym",
                          "sen_arriaga_inst", "sen_arriaga_inst2",
                          "chandrasekaran_ii","chandrasekaran_ii_sym",
                          "sen_chandrasekaran_ii", "sen_chandrasekaran_ii_inst",
                          "sen_chandrasekaran_ii_inst2", "sen_chandrasekaran_ii_sym",
                          "chandrasekaran_iii","chandrasekaran_iii_sym",
                          "sen_chandrasekaran_iii", "sen_chandrasekaran_iii_inst",
                          "sen_chandrasekaran_iii_inst2", "sen_chandrasekaran_iii_sym",
                          "lopez_ruzicka", "lopez_ruzicka_sym",
                          "sen_lopez_ruzicka", "sen_lopez_ruzicka_sym",
                          "sen_lopez_ruzicka_inst", "sen_lopez_ruzicka_inst2",
                          "horiuchi", "stepwise", "numerical"),
                      several.ok = FALSE)
  # assign generic decomp function name
  dec_fun <- switch(method,
                    "lifetable" = sen_e0_mx_lt,
                    "arriaga" = arriaga,
                    "arriaga_sym" = arriaga_sym,
                    "sen_arriaga" = sen_arriaga,
                    "sen_arriaga_sym" = sen_arriaga_sym,
                    "sen_arriaga_inst" = sen_arriaga_instantaneous,
                    "sen_arriaga_inst2" = sen_arriaga_instantaneous2,
                    "chandrasekaran_ii" = chandrasekaran_II,
                    "sen_chandrasekaran_ii" = sen_chandrasekaran_II,
                    "sen_chandrasekaran_ii_inst" = sen_chandrasekaran_II_instantaneous,
                    "sen_chandrasekaran_ii_inst2" = sen_chandrasekaran_II_instantaneous2,
                    "chandrasekaran_iii"= chandrasekaran_III,
                    "sen_chandrasekaran_iii" = sen_chandrasekaran_III,
                    "sen_chandrasekaran_iii_inst" = sen_chandrasekaran_III_instantaneous,
                    "sen_chandrasekaran_iii_inst2" = sen_chandrasekaran_III_instantaneous2,
                    "lopez_ruzicka" = lopez_ruzicka,
                    "lopez_ruzicka_sym" = lopez_ruzicka_sym,
                    "sen_lopez_ruzicka" = sen_lopez_ruzicka,
                    "sen_lopez_ruzicka_sym" = sen_lopez_ruzicka_sym,
                    "sen_lopez_ruzicka_inst" = sen_lopez_ruzicka_instantaneous,
                    "sen_lopez_ruzicka_inst2" = sen_lopez_ruzicka_instantaneous2,
                    "numerical" = sen_num)
  # ----------------------------------------------------------------- #
  # method block for stepwise, numerical, horiuchi                    #
  # should handle vector and matrix cases                             #
  # ----------------------------------------------------------------- #
  if (method %in% c("stepwise","horiuchi")){

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
  if (method %in% c("arriaga",
                    "arriaga_sym",
                    "chandrasekaran_ii",
                    "chandrasekaran_iii",
                    "lopez_ruzicka",
                    "lopez_ruzicka_sym")){
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
  if (method %in% c("lifetable", "sen_arriaga_inst", "sen_arriaga_inst2",
                    "sen_chandrasekaran_ii_inst",
                    "sen_chandrasekaran_ii_inst2",
                    "sen_chandrasekaran_iii_inst", "sen_chandrasekaran_iii_inst2",
                    "sen_lopez_ruzicka_inst", "sen_lopez_ruzicka_inst2","numerical")){
    # First, either we optimize or we insist on the midpoint between mx1 and mx2
    # (1) yes, we optimize the averaging
    if (opt){
      # second case: we have causes of death
      if (!is.null(n_causes)){
        mx1_all <- rowSums(mx1)
        mx2_all <- rowSums(mx2)
        sen     <- sen_min(mx1 = mx1_all,
                           mx2 = mx2_all,
                           age = age,
                           nx = nx,
                           sex1 = sex1,
                           sex2 = sex1, # see diff sex solution at start
                           closeout = closeout,
                           sen_fun = dec_fun,
                           tol = tol)
        delta   <- mx2 - mx1
      } else {
        # all-cause mortality only
        sen     <- sen_min(mx1 = mx1,
                           mx2 = mx2,
                           age = age,
                           nx = nx,
                           sex1 = sex1,
                           sex2 = sex1, # see diff sex solution at start
                           closeout = closeout,
                           sen_fun = dec_fun,
                           tol = tol)
        delta   <- mx2 - mx1
      }
    } else {
      # (2) we don't optimize, i.e. we take average mx
      if (!is.null(n_causes)){
        mx_avg <- (rowSums(mx2) + rowSums(mx1)) / 2
      } else {
        mx_avg <- (mx1 + mx2) / 2
      }
        delta  <- mx2 - mx1
        sen <- dec_fun(mx = mx_avg,
                       age = age,
                       nx = nx,
                       sex = sex1,
                       closeout = closeout)


    }
    # all four of these cases back out the decomp in the same way;
    decomp <- sen * delta

    Delta <-
      mx_to_e0(mx2,
               age = age,
               nx = nx,
               sex = sex1,
               closeout = closeout) -
      mx_to_e0(mx1,
               age = age,
               nx = nx,
               sex = sex1,
               closeout = closeout)
    if (abs(sum(decomp) - Delta) > .1){
      if (opt){
        warning("\nYou used a sensitivity-based method (",method,") but still have a decomposition residual of",round(sum(decomp) - Delta,3),". Consider comparing with other methods\n")
      }
      if (!opt){
        warning("\nYou used a sensitivity-based method (",method,") evaluated at the midpoint between mx1 and mx1, giving a decomposition residual of",round(sum(decomp) - Delta,3),". Consider comparing with other methods or setting opt=TRUE\n")
      }
    }
  }
  # ----------------------------------------------------------------- #
  # sensitivity methods using both mx1 and mx2,                       #
  # basically all are conversions of direct methods                   #
  # ----------------------------------------------------------------- #
  if (method %in% c("sen_arriaga","sen_arriaga_sym",
                    "sen_chandrasekaran_ii","sen_chandrasekaran_iii",
                    "sen_lopez_ruzicka","sen_lopez_ruzicka_sym")){
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

    } else if(x$method %in% c( "sen_arriaga", "sen_arriaga_sym" , "sen_arriaga_inst" , "sen_arriaga_inst2" , "sen_chandrasekaran_ii" , "sen_chandrasekaran_ii_inst" , "sen_chandrasekaran_ii_inst2" , "sen_chandrasekaran_iii" , "sen_chandrasekaran_iii_inst" , "sen_chandrasekaran_iii_inst2" , "sen_lopez_ruzicka" , "sen_lopez_ruzicka_sym" , "sen_lopez_ruzicka_inst" , "sen_lopez_ruzicka_sym_inst2")){

      values <- c("sen_arriaga", "sen_arriaga_sym", "sen_arriaga_inst",
                  "sen_arriaga_inst2", "sen_chandrasekaran_ii",
                  "sen_chandrasekaran_ii_inst", "sen_chandrasekaran_ii_inst2",
                  "sen_chandrasekaran_iii", "sen_chandrasekaran_iii_inst",
                  "sen_chandrasekaran_iii_inst2", "sen_lopez_ruzicka",
                  "sen_lopez_ruzicka_sym", "sen_lopez_ruzicka_inst",
                  "sen_lopez_ruzicka_inst2")

      names <- c("arriaga", "arriaga_sym", "arriaga_inst", "arriaga_inst2",
                 "chandrasekaran_ii", "chandrasekaran_ii_inst", "chandrasekaran_ii_inst2",
                 "chandrasekaran_iii", "chandrasekaran_iii_inst", "chandrasekaran_iii_inst2",
                 "lopez_ruzicka", "lopez_ruzicka_sym", "lopez_ruzicka_inst", "lopez_ruzicka_sym__inst2")
      dicc <- names[which(x$method == values)[1]]

      cat(paste("Estimated the", dicc, "sensitivity Life-Expectancy decomposition method."))

  }
    } else if(is.matrix(x$LEdecomp)){
      if(x$method %in%c("arriaga", "arriaga_sym",  "chandrasekaran_ii" , "lopez_ruzicka","lopez_ruzicka_sym", "horiuchi","stepwise")){

        cat(paste("Estimated the", x$method, "cause-of-death Life-Expectancy decomposition method."))

      } else if(x$method %in%c( "lifetable","sen_arriaga", "sen_arriaga_sym" , "sen_arriaga_inst" , "sen_arriaga_inst2" , "sen_chandrasekaran_ii" , "sen_chandrasekaran_ii_inst" , "sen_chandrasekaran_ii_inst2", "sen_chandrasekaran_iii" , "sen_chandrasekaran_iii_inst" , "sen_chandrasekaran_iii_inst2" , "sen_lopez_ruzicka" , "sen_lopez_ruzicka_sym" , "sen_lopez_ruzicka_inst" , "sen_lopez_ruzicka_sym_inst2","numerical")){

        values <- c("sen_arriaga", "sen_arriaga_sym", "sen_arriaga_inst",
                    "sen_arriaga_inst2", "sen_chandrasekaran_ii",
                    "sen_chandrasekaran_ii_inst","sen_chandrasekaran_ii_inst2",
                    "sen_chandrasekaran_iii", "sen_chandrasekaran_iii_inst",
                    "sen_chandrasekaran_iii_inst2", "sen_lopez_ruzicka",
                    "sen_lopez_ruzicka_sym", "sen_lopez_ruzicka_inst",
                    "sen_lopez_ruzicka_sym__inst2")

        names <- c("arriaga", "arriaga_sym", "arriaga_inst", "arriaga_inst2",
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












