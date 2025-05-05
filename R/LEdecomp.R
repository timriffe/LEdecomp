
#' Function for applying different Life-Expectancy decomposition and sensitivity methods
#' @description A variety of exact or asymptotically exact life expectancy decomposition methods are implemented. Also, several life-expectancy descomposition sensitivity methods are implemented to answer how each age will change with an increase/decrease in life expectancy.
#' These include the lifetable response experiment, using three potential sensitivity functions (lifetable-based `ltre_lt`, an instantaneous Arriaga-derived sensitivity `ltre_arriaga_instantaneous`, or numerical sensitivity `ltre_numerical`), the directional Arriaga method `arriaga`, a symmetrical Arriaga method `arriaga_symmetrical`.
#' These all give similar results, except for directional `arriaga`, which is the most different.
#'
#' @param mx1 numeric. age-structured mortality rates for population 1.
#' @param mx2 numeric. age-structured mortality rates for population 2, must be same length as `mx1`
#' @param age integer. lower bound of each age group, must be same length as `mx1`
#' @param sex1 character. `"m"`,`"f"`, or `"t"`, affects a0 treatment.
#' @param sex2 character. `"m"`,`"f"`, or `"t"`, affects a0 treatment.
#' @param method character. `"lifetable"`, `"arriaga"'`, `"arriaga_sym"'`, `"sen_arriaga"'`, `"sen_arriaga_sym"'`, `"sen_arriaga_inst"'`, `"sen_arriaga_inst2"'`, `"chandrasekaran_ii"'`, `"sen_chandrasekaran_ii"'`, `"sen_chandrasekaran_ii_inst"'`, `"sen_chandrasekaran_ii_inst2"'`, `"chandrasekaran_iii"'`, `"sen_chandrasekaran_iii"'`, `"sen_chandrasekaran_iii_inst"'`, `"sen_chandrasekaran_iii_inst2"'`, `"lopez_ruzicka"'`, `"lopez_ruzicka_sym"'`, `"sen_lopez_ruzicka"'`, `"sen_lopez_ruzicka_sym"'`, `"sen_lopez_ruzicka_inst"'`, `"sen_lopez_ruzicka_inst2"'`, `"horiuchi"'`, `"stepwise"'`, `"numerical"`
#' @param closeout logical. Do we handle closeout, or truncate at top age.
#' @param opt logical, default `TRUE`. For sensitivity-based decomposition, shall we optimize the rate averaging to eliminate the decomposition residual?
#' @param func A function specified by the user. This must be able to take the vectors \code{mx1} or \code{mx2} as its argument, and to return the value of the function. This `func` is only applied for `method` = c(`"stepwise"`, `"horiuchi"`).
#' @param tol numeric, default `1e-10`. Tolerance parameter for rate averaging optimization.
#' @param Num_Intervals numeric, default `FALSE`. The number of intervals to integrate over, which is used only for `method` = `"arriaga"` or the sensitvity methods of Arriaga.
#' @param symmetrical logical value, default `FALSE`. This parameter is used for `method` = `"stepwise"` and shall we want to average the results of replacing 1 with 2 and 2 with 1?
#' @param direction character, default `'up'`. This parameter is used for `method` = `"stepwise"`, and corresponds to one of \code{"up"}, \code{"down"}, or \code{"both"}.
#' @param perturb numeric, default `1e-6`, perturbation value that is a numeric constant, and a very small number.
#' @param \dots optional arguments passed to `numDeriv::grad()`
#'
#' @return A list with class \code{"LEdecomp"} including different component in the Life Expectancy decomposition method:
#' * `mx1` age-structured mortality rates for population 1, provided by the user.
#' * `mx2` age-structured mortality rates for population 2, provided by the user.
#' * `age` age, integer. lower bound of each age group,  provided by the user.
#' * `sex1` gender selected for first population, provided by the user.
#' * `sex2` gender selected for second population, provided by the user.
#' * `method` method selected by the user to apply the life-expectancy decomposition or sensitivity analysis.
#' * `closeout` how the user have decide to handle closeout or truncate at top age.
#' * `opt` for sensitivity process how to optimize
#' * `func` function selected for estimate the life-expectancy decomposition and the sensitivity value.
#' * `tol` level of tolerance.
#' * `Num_Intervals` the number of intervals to integrate over.
#' * `symmetrical` shall the user wants to average the results of replacing 1 with 2 and 2 with 1.
#' * `direction` ----
#' * `perturb` perturbation value that is a numeric constant, and a very small number.
#' * `sens` a vector/matrix of the sensitivity life expectancy decomposition effects that is organized in the same way as \code{mx1} and \code{mx2}.
#' * `LEdecomp` a vector/matrix of the life expectancy decomposition effects that is organized in the same way as \code{mx1} and \code{mx2}.
#'
#' @seealso \\code{\link{lifetable}}, \code{\link{arriaga}},
#' \code{\link{arriaga_sym}}, \code{\link{sen_arriaga}}, \code{\link{sen_arriaga_sym}},
#' \code{\link{sen_arriaga_inst}}, \code{\link{sen_arriaga_inst2}},
#' \code{\link{chandrasekaran_ii}}, \code{\link{sen_chandrasekaran_ii}},
#' \code{\link{sen_chandrasekaran_ii_inst}}, \code{\link{sen_chandrasekaran_ii_inst2}},
#' \code{\link{chandrasekaran_iii}}, \code{\link{sen_chandrasekaran_iii}},
#' \code{\link{sen_chandrasekaran_iii_inst}}, \code{\link{sen_chandrasekaran_iii_inst2}},
#' \code{\link{lopez_ruzicka}}, \code{\link{lopez_ruzicka_sym}},
#' \code{\link{sen_lopez_ruzicka}}, \code{\link{sen_lopez_ruzicka_sym}},
#' \code{\link{sen_lopez_ruzicka_inst}}, \code{\link{sen_lopez_ruzicka_inst2}},
#' \code{\link{horiuchi}}, \code{\link{stepwise}}, \code{\link{numerical}}
#'
#' @references
#' \insertRef{arriaga1984measuring}{coddecomp}
#' \insertREf{Chandrasekaran1986}{coddecomp}
#' \insertRef{preston2000demography}{coddecomp}
#' \insertREf{Ponnapalli2005}{coddecomp}
#'
#' @import DemoTools
#' @import data.table
#' @import dplyr
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
#' LEdecomp(mx1,mx2,age=x,sex1='t',method = "arriaga_instantaneous")
#' @export
LEdecomp <- function(mx1, mx2,
                     age,
                     sex1 = 't', sex2 = sex1,
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
                                "horiuchi",
                                "stepwise", "numerical"),
                     closeout = TRUE, opt = TRUE, func = FALSE,
                     tol = 1e-10, Num_Intervals = FALSE,
                     symmetrical = TRUE, direction = "up",
                     perturb = 1e-6, ...){

  stopifnot(is.vector(mx1) | is.matrix(mx1))
  stopifnot(length(mx1) == length(mx2))

  #check that the main information has been provided
  if(is.null(mx1) || is.null(mx2) || is.null(age)){
    warning("Arguments mx1, mx2, age, and sex1, need to be provided.")
  }

  names_causes <- colnames(mx1)

  # handle dimensions
  deez_dims <- dim(mx1)
  nages <- length(age)
  nmx <- length(mx1)
  if (nages < nmx){
    ncauses <- nmx / nages
    ncauses <- round(ncauses)
    stopifnot((ncauses * nages) == nmx)
    dim(mx1) <- c(nages,ncauses)
    dim(mx2) <- c(nages,ncauses)
  }

  #if the user donot provide method, automatically,
  #the standard version of arriage method will be executed
  if(exists("method") == "FALSE"){
    method = "arriaga"
  }

  # add optimization option.
  method <- tolower(method)
  method <- match.arg(method,
                      choices =
                        c("lifetable",
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
                          "horiuchi",
                          "stepwise", "numerical"))
  if(is.vector(mx1) & length(mx1) == length(age)){
    if(method == "lifetable"){
      mx <- (mx2+mx1)/2
      sen <- sen_e0_mx_lt(mx = mx, age = x,
                          sex = sex1, closeout = closeout)

      decomp <- sen*(mx2 - mx1)

    }else if(method == "arriaga"){
      decomp <- arriaga(mx1 = mx1, mx2 = mx2,
                        age = age,
                        sex1 = sex1, sex2 = sex2,
                        closeout = closeout)
      sen <- NULL

    }else if(method == "arriaga_sym"){
      decomp <- arriaga_sym(mx1 = mx1, mx2 = mx2,
                            age = age,
                            sex1 = sex1, sex2 = sex2,
                            closeout = closeout)
      sen <- NULL

    }else if(method == "sen_arriaga"){
      prev <- sen_arriaga(mx1 = mx1, mx2 = mx2,
                          age = age,
                          sex1 = sex1, sex2 = sex2,
                          closeout = closeout)

      sen <- prev
      decomp <- prev*(mx2 - mx1)

    }else if(method == "sen_arriaga_sym"){
      prev <- sen_arriaga_sym(mx1 = mx1, mx2 = mx2,
                              age = age,
                              sex1 = sex1, sex2 = sex2,
                              closeout = closeout)

      sen <- prev
      decomp <- prev*(mx2 - mx1)

    }else if(method == "sen_arriaga_inst"){
      mx <- (mx1 + mx2)/2
      sen <- sen_arriaga_instantaneous(mx, age = x,
                                       perturb = perturb,
                                       closeout = closeout)

      decomp <- sen*(mx2 - mx1)

    }else if(method == "sen_arriaga_inst2"){
      mx <- (mx1 + mx2)/2
      sen <- sen_arriaga_instantaneous2(mx, age = x,
                                        perturb = perturb,
                                        closeout = closeout)

      decomp <- sen*(mx2 - mx1)


    }else if(method == "chandrasekaran_ii"){
      sen <- NULL
      decomp <- chandrasekaran_II(mx1 = mx1, mx2 = mx2,
                                  age = x,
                                  sex1 = sex1, sex2 = sex2,
                                  closeout = closeout)
    }else if(method == "sen_chandrasekaran_ii"){
      prev <- sen_chandrasekaran_II(mx1 = mx1, mx2 = mx2,
                                    age = age,
                                    sex1 = sex1, sex2 = sex2,
                                    closeout = closeout)
      sen <- prev
      decomp <- prev*(mx2 - mx1)

    }else if(method == "sen_chandrasekaran_ii_inst"){
      mx <- (mx1 + mx2)/2
      prev <- sen_chandrasekaran_II_instantaneous(mx = mx, age = age,
                                                  sex = sex1, closeout = closeout)
      sen <- prev
      decomp <- prev*(mx2 - mx1)

    }else if(method == "sen_chandrasekaran_ii_inst2"){
      mx <- (mx1 + mx2)/2
      prev <- sen_chandrasekaran_II_instantaneous2(mx = mx, age = age,
                                                   sex = sex1, closeout = closeout)
      sen <- prev
      decomp <- prev*(mx2 - mx1)

    }else if(method == "chandrasekaran_iii"){
      sen <- NULL
      decomp <- chandrasekaran_III(mx1 = mx1, mx2 = mx2,
                                   age = x,
                                   sex1 = sex1, sex2 = sex2,
                                   closeout = closeout)

    }else if(method == "sen_chandrasekaran_iii"){
      prev <- sen_chandrasekaran_III(mx1 = mx1, mx2 = mx2,
                                     age = age,
                                     sex1 = sex1, sex2 = sex2,
                                     closeout = closeout)
      sen <- prev
      decomp <- prev*(mx2 - mx1)

    }else if(method == "sen_chandrasekaran_iii_inst"){
      mx <- (mx1 + mx2)/2
      prev <- sen_chandrasekaran_III_instantaneous(mx = mx, age = age,
                                                   sex = sex1, closeout = closeout)
      sen <- prev
      decomp <- prev*(mx2 - mx1)

    }else if(method == "sen_chandrasekaran_iii_inst2"){
      mx <- (mx1 + mx2)/2
      prev <- sen_chandrasekaran_III_instantaneous2(mx = mx, age = age,
                                                    sex = sex1, closeout = closeout)
      sen <- prev
      decomp <- prev*(mx2 - mx1)

    }else if(method == "lopez_ruzicka"){
      sen <- NULL
      decomp <- lopez_ruzicka(mx1 = mx1, mx2 = mx2,
                              age = x,
                              sex1 = sex1, sex2 = sex2,
                              closeout = closeout)
    }else if(method == "lopez_ruzicka_sym"){
      sen <- NULL
      decomp <- lopez_ruzicka_sym(mx1 = mx1, mx2 = mx2,
                                  age = x,
                                  sex1 = sex1, sex2 = sex2,
                                  closeout = closeout)

    }else if(method == "sen_lopez_ruzicka"){
      prev <- sen_lopez_ruzicka(mx1 = mx1, mx2 = mx2,
                                age = age,
                                sex1 = sex1, sex2 = sex2,
                                closeout = closeout)

      sen <- prev
      decomp <- prev*(mx2 - mx1)

    }else if(method == "sen_lopez_ruzicka_sym"){
      prev <- sen_lopez_ruzicka_sym(mx1 = mx1, mx2 = mx2,
                                    age = age,
                                    sex1 = sex1, sex2 = sex2,
                                    closeout = closeout)

      sen <- prev
      decomp <- prev*(mx2 - mx1)

    }else if(method == "sen_lopez_ruzicka_inst"){
      mx <- (mx1 + mx2)/2
      sen <- sen_lopez_ruzicka_instantaneous(mx, age = x,
                                             perturb = perturb,
                                             closeout = closeout)

      decomp <- sen*(mx2 - mx1)

    }else if(method == "sen_lopez_ruzicka_inst2"){
      mx <- (mx1 + mx2)/2
      sen <- sen_lopez_ruzicka_instantaneous2(mx, age = x,
                                              perturb = perturb,
                                              closeout = closeout)

      decomp <- sen*(mx2 - mx1)

    }else if(method == "horiuchi"){

      if(Num_Intervals == "FALSE"){
        stop("The object 'Num_Intervals' must be provided and have to be a number.")
      } else if(!is.numeric(Num_Intervals)){
        stop("'Num_Intervals' has to be a numeric variable.")
      }
      if(is.function(func) == "FALSE"){
        stop("The object 'func' must be provided and have to be a function")
      }
      sen <- NULL
      decomp <- DemoDecomp::horiuchi(func = func,
                                     pars1 = mx1, pars2 = mx2,
                                     age = age,
                                     N = Num_Intervals)

    }else if(method == "stepwise"){
      if(symmetrical == "FALSE"){
        stop("The object 'symmetrical' must be 'TRUE'.")
      }
      if(!(direction %in% c("up", "both", "down"))){
        stop("The object 'direction' must be 'up', 'both' or 'down', please modify it." )
      }
      if(is.function(func) == "FALSE"){
        stop("The object 'func' must be provided and have to be a function")
      }
      if(is.vector(mx1) == TRUE){
        dims <- c(1, length(mx1))
      }else if(is.matrix(mx1) == TRUE | is.data.frame(mx1) == TRUE){
        dims <- dim(mx1)
      }

      sen <- NULL
      decomp <- DemoDecomp::stepwise_replacement(func = func,
                                                 pars1 = mx1, pars2 = mx2,
                                                 symmetrical = symmetrical,
                                                 direction = direction,
                                                 age = age, ...)

    }else if(method == "numerical"){
      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #don't work
      mx <- (mx1 + mx2)/2
      sen <- sen_num(mx, age=age, sex=sex1, closeout=TRUE, ...)
      decomp <- sen*(mx2 - mx1)

    }
  }else if(is.matrix(mx1) & all(dim(mx1) == dim(mx2))){
    colnames(mx1) <- names_causes
    rownames(mx1) <- age

    mx1_all <- rowSums(mx1)
    mx2_all <- rowSums(mx2)
    mx1_causes <- mx1
    mx2_causes <- mx2
    colnames(mx1_causes) <- colnames(mx1)
    delta_causes <- (mx2_causes - mx1_causes)
    delta_all <- (mx2_all - mx1_all)

    if(method == "lifetable"){
      #REVISAR!!!!!!!!!!!
      mx <- (mx2_all + mx1_all)/2
      sen <- sen_e0_mx_lt(mx = mx, age = x,
                          sex = sex1, closeout = closeout)

      decomp <- sen*delta_causes
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "arriaga" | method == "sen_arriaga"){
      sen <- sen_arriaga(mx1 = mx1_all, mx2 = mx2_all,
                         age = age,
                         sex1 = sex1, sex2 = sex2,
                         closeout = closeout)

      decomp <- sen*delta_causes
      if(sum(abs(delta_all) < tol) > 1){
        warning("There is at least one age with no differences in the cause-of-death between two considered populations.")
      }
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "arriaga_sym" | method == "sen_arriaga_sym"){
      sen <- sen_arriaga_sym(mx1 = mx1_all, mx2 = mx2_all,
                             age = age,
                             sex1 = sex1, sex2 = sex2,
                             closeout = closeout)
      decomp <- sen*delta_causes
      if(sum(abs(delta_all) < tol) > 1){
        warning("There is at least one age with no differences in the cause-of-death between two considered populations.")
      }
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "sen_arriaga_inst"){
      mx <- (mx1_all + mx2_all)/2
      sen <- sen_arriaga_instantaneous(mx, age = x,
                                       perturb = perturb,
                                       closeout = closeout)

      decomp <- sen*delta_causes
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "sen_arriaga_inst2"){
      mx <- (mx1_all + mx2_all)/2
      sen <- sen_arriaga_instantaneous2(mx, age = x,
                                        perturb = perturb,
                                        closeout = closeout)

      decomp <- sen*delta_causes
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "chandrasekaran_ii" | method == "sen_chandrasekaran_ii"){
      sen <- sen_chandrasekaran_II(mx1 = mx1_all, mx2 = mx2_all,
                                   age = x,
                                   sex1 = sex1, sex2 = sex2,
                                   closeout = closeout)
      decomp <- sen*delta_causes
      if(sum(abs(delta_all) < tol) > 1){
        warning("There is at least one age with no differences in the cause-of-death between two considered populations.")
      }
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "sen_chandrasekaran_ii_inst"){
      mx <- (mx1_all + mx2_all)/2
      prev <- sen_chandrasekaran_II_instantaneous(mx = mx, age = age,
                                                  sex = sex1, closeout = closeout)
      sen <- prev
      decomp <- prev*delta_causes
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "sen_chandrasekaran_ii_inst2"){
      mx <- (mx1_all + mx2_all)/2
      prev <- sen_chandrasekaran_II_instantaneous2(mx = mx, age = age,
                                                   sex = sex1, closeout = closeout)
      sen <- prev
      decomp <- prev*delta_causes
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "chandrasekaran_iii" | method == "sen_chandrasekaran_iii"){
      sen <- sen_chandrasekaran_III(mx1 = mx1_all, mx2 = mx2_all,
                                    age = x,
                                    sex1 = sex1, sex2 = sex2,
                                    closeout = closeout)
      decomp <- sen*delta_causes
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

      if(sum(abs(delta_all) < tol) > 1){
        warning("There is at least one age with no differences in the cause-of-death between two considered populations.")
      }

    }else if(method == "sen_chandrasekaran_iii_inst"){
      mx <- (mx1_all + mx2_all)/2
      prev <- sen_chandrasekaran_III_instantaneous(mx = mx, age = age,
                                                   sex = sex1, closeout = closeout)
      sen <- prev
      decomp <- prev*delta_causes
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "sen_chandrasekaran_iii_inst2"){
      mx <- (mx1_all + mx2_all)/2
      prev <- sen_chandrasekaran_III_instantaneous2(mx = mx, age = age,
                                                    sex = sex1, closeout = closeout)
      sen <- prev
      decomp <- prev*delta_causes
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "lopez_ruzicka" | method == "sen_lopez_ruzicka"){
      prev <- sen_lopez_ruzicka(mx1 = mx1_all, mx2 = mx2_all,
                                age = age,
                                sex1 = sex1, sex2 = sex2,
                                closeout = closeout)

      sen <- prev
      decomp <- prev*delta_causes
      if(sum(abs(delta_all) < tol) > 1){
        warning("There is at least one age with no differences in the cause-of-death between two considered populations.")
      }
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "lopez_ruzicka_sym" | method == "sen_lopez_ruzicka_sym"){
      prev <- sen_lopez_ruzicka_sym(mx1 = mx1_all, mx2 = mx2_all,
                                    age = age,
                                    sex1 = sex1, sex2 = sex2,
                                    closeout = closeout)

      sen <- prev
      decomp <- prev*delta_causes
      if(sum(abs(delta_all) < tol) > 1){
        warning("There is at least one age with no differences in the cause-of-death between two considered populations.")
      }
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "sen_lopez_ruzicka_inst"){
      mx <- (mx1_all + mx2_all)/2
      sen <- sen_lopez_ruzicka_instantaneous(mx, age = x,
                                             perturb = perturb,
                                             closeout = closeout)

      decomp <- sen*delta_causes
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "sen_lopez_ruzicka_inst2"){
      mx <- (mx1_all + mx2_all)/2
      sen <- sen_lopez_ruzicka_instantaneous2(mx, age = x,
                                              perturb = perturb,
                                              closeout = closeout)

      decomp <- sen*delta_causes
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "horiuchi"){
    #!!!!!!!!!!!!!!!!!!! REPASAR
      if(Num_Intervals == "FALSE"){
        stop("The object 'Num_Intervals' must be provided and have to be a number.")
      } else if(!is.numeric(Num_Intervals)){
        stop("'Num_Intervals' has to be a numeric variable.")
      }
      if(is.function(func) == "FALSE"){
        stop("The object 'func' must be provided and have to be a function")
      }
      sen <- NULL
      decomp <- DemoDecomp::horiuchi(func = func,
                                     pars1 = mx1_all, pars2 = mx2_all,
                                     age = age,
                                     N = Num_Intervals)
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "stepwise"){
      #!!!!!!!!!!!!!!!!!!! REPASAR
      if(symmetrical == "FALSE"){
        stop("The object 'symmetrical' must be 'TRUE'.")
      }
      if(!(direction %in% c("up", "both", "down"))){
        stop("The object 'direction' must be 'up', 'both' or 'down', please modify it." )
      }
      if(is.function(func) == "FALSE"){
        stop("The object 'func' must be provided and have to be a function")
      }
      if(is.vector(mx1) == TRUE){
        dims <- c(1, length(mx1))
      }else if(is.matrix(mx1) == TRUE | is.data.frame(mx1) == TRUE){
        dims <- dim(mx1)
      }

      sen <- NULL
      decomp <- DemoDecomp::stepwise_replacement(func = func,
                                                 pars1 = mx1, pars2 = mx2,
                                                 symmetrical = symmetrical,
                                                 direction = direction,
                                                 age = age, ...)
      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }else if(method == "numerical"){
      #!!!!!!!!!!!!!!!!!!! REPASAR
      #don't work
      mx <- (mx1_all + mx2_all)/2
      sen <- sen_num(mx, age=age, sex=sex1, closeout=TRUE, ...)
      decomp <- sen*delta_causes

      colnames(decomp) <- colnames(mx1_causes)
      rownames(decomp) <- age

    }
  }

  return <- list("mx1" = mx1,
                 "mx2" = mx2,
                 "age" = age,
                 "sex1" = sex1,
                 "sex2" = sex2,
                 "method" = method,
                 "closeout" = closeout,
                 "opt" = opt,
                 "func" = func,
                 "tol" = tol,
                 "Num_Intervals" = Num_Intervals,
                 "symmetrical" = symmetrical,
                 "direction" = direction,
                 "perturb" = perturb,
                 "sens" = sen,
                 "LEdecomp" = decomp)
  class(return) <- "LEdecomp"
  return

}

sen_num <- function(mx,age,sex='t',closeout=TRUE,...){
  numDeriv::grad(mx_to_e0,mx,age=age,sex=sex,closeout=closeout,...)
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
    if(x$method == "arriaga" | x$method == "arriaga_sym" |
       x$method == "chandrasekaran_ii" | x$method == "chandrasekaran_iii" |
       x$method == "lopez_ruzicka" | x$method == "lopez_ruzicka_sym" |
       x$method == "horiuchi" | x$method == "stepwise" | x$method == "numerical"){

      cat(paste("Estimated the", x$method, "Life-Expectancy decomposition method."))

    } else if(x$method == "sen_arriaga" | x$method == "sen_arriaga_sym" |
              x$method == "sen_arriaga_inst" | x$method == "sen_arriaga_inst2" |
              x$method == "sen_chandrasekaran_ii" | x$method == "sen_chandrasekaran_ii_inst" | x$method == "sen_chandrasekaran_ii_inst2" |
              x$method == "sen_chandrasekaran_iii" | x$method == "sen_chandrasekaran_iii_inst" | x$method == "sen_chandrasekaran_iii_inst2" |
              x$method == "sen_lopez_ruzicka" | x$method == "sen_lopez_ruzicka_sym" |
              x$method == "sen_lopez_ruzicka_inst" | x$method == "sen_lopez_ruzicka_sym__inst2"){

      values <- c("sen_arriaga", "sen_arriaga_sym", "sen_arriaga_inst", "sen_arriaga_inst2",
        "sen_chandrasekaran_ii", "sen_chandrasekaran_ii_inst", "sen_chandrasekaran_ii_inst2",
        "sen_chandrasekaran_iii", "sen_chandrasekaran_iii_inst", "sen_chandrasekaran_iii_inst2",
        "sen_lopez_ruzicka", "sen_lopez_ruzicka_sym", "sen_lopez_ruzicka_inst", "sen_lopez_ruzicka_sym__inst2")

      names <- c("arriaga", "arriaga_sym", "arriaga_inst", "arriaga_inst2",
                 "chandrasekaran_ii", "chandrasekaran_ii_inst", "chandrasekaran_ii_inst2",
                 "chandrasekaran_iii", "chandrasekaran_iii_inst", "chandrasekaran_iii_inst2",
                 "lopez_ruzicka", "lopez_ruzicka_sym", "lopez_ruzicka_inst", "lopez_ruzicka_sym__inst2")
      dicc <- names[which(x$method == values)[1]]

      cat(paste("Estimated the", dicc, "sensitivity Life-Expectancy decomposition method."))

  }
    } else if(is.matrix(x$LEdecomp)){
      if(x$method == "arriaga" | x$method == "arriaga_sym" |
         x$method == "chandrasekaran_ii" | x$method == "chandrasekaran_iii" |
         x$method == "lopez_ruzicka" | x$method == "lopez_ruzicka_sym" |
         x$method == "horiuchi" | x$method == "stepwise" | x$method == "numerical"){

        cat(paste("Estimated the", x$method, "cause-of-death Life-Expectancy decomposition method."))

      } else if(x$method == "sen_arriaga" | x$method == "sen_arriaga_sym" |
                x$method == "sen_arriaga_inst" | x$method == "sen_arriaga_inst2" |
                x$method == "sen_chandrasekaran_ii" | x$method == "sen_chandrasekaran_ii_inst" | x$method == "sen_chandrasekaran_ii_inst2" |
                x$method == "sen_chandrasekaran_iii" | x$method == "sen_chandrasekaran_iii_inst" | x$method == "sen_chandrasekaran_iii_inst2" |
                x$method == "sen_lopez_ruzicka" | x$method == "sen_lopez_ruzicka_sym" |
                x$method == "sen_lopez_ruzicka_inst" | x$method == "sen_lopez_ruzicka_sym__inst2"){

        values <- c("sen_arriaga", "sen_arriaga_sym", "sen_arriaga_inst", "sen_arriaga_inst2",
                    "sen_chandrasekaran_ii", "sen_chandrasekaran_ii_inst", "sen_chandrasekaran_ii_inst2",
                    "sen_chandrasekaran_iii", "sen_chandrasekaran_iii_inst", "sen_chandrasekaran_iii_inst2",
                    "sen_lopez_ruzicka", "sen_lopez_ruzicka_sym", "sen_lopez_ruzicka_inst", "sen_lopez_ruzicka_sym__inst2")

        names <- c("arriaga", "arriaga_sym", "arriaga_inst", "arriaga_inst2",
                   "chandrasekaran_ii", "chandrasekaran_ii_inst", "chandrasekaran_ii_inst2",
                   "chandrasekaran_iii", "chandrasekaran_iii_inst", "chandrasekaran_iii_inst2",
                   "lopez_ruzicka", "lopez_ruzicka_sym", "lopez_ruzicka_inst", "lopez_ruzicka_sym__inst2")
        dicc <- names[which(x$method == values)[1]]

        cat(paste("Estimated the", dicc, "cause-of-death sensitivity Life-Expectancy decomposition method.\n"))

      }

  }

  cat(paste("\nThe total age-different effect is:", round(sum(x$LEdecomp), 4)))
}












