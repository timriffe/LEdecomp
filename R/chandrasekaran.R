# This script contains functions that implement and adapt the
# Chandrasekaran decompositon approach

#' @title II approach of Chandrasekaran decomposition approach
#' @description Following the notation given in Ponnapalli (2005), and the decomposition method can written as:
#' \deqn{
#' _{n}\Delta_{x} =
#' \frac{(e_x^2 - e_x^1)(l_x^2 + l_x^1)}{2}
#' -
#' \frac{(e_{x+n}^2 - e_{x+n}^1)(l_{x+n}^2 + l_{x+n}^1)}{2}
#' -
#' \frac{_{n}L_{x}^{1}}{l_{x}^{1}}
#' +
#' \frac{T_{x+n}^{2}}{l_{0}^{1}} \left( \frac{l_{x}^{1}}{l_{x}^{2}} - \frac{l_{x+n}^{1}}{l_{x+n}^{2}} \right)
#' }

#' where \eqn{_{n}\Delta_{x}} is the contribution of rate differences in age \eqn{x} to the difference in life expectancy implied by `mx1` and `mx2`. This formula can be averaged between ‘effect interaction deferred’ and ‘effect interaction forwarded’ from the Ponnapalli (2005).

#' @param mx1 numeric vector of the mortality rates (central death rates) for population 1
#' @param mx2 numeric vector of the mortality rates (central death rates) for population 2
#' @param age integer vector of the lower bound of each age group (currently only single ages supported)
#' @param nx integer vector of age intervals, default 1.
#' @param sex1 character either the sex for population 1: Male (`"m"`), Female (`"f"`), or Total (`"t"`)
#' @param sex2 character either the sex for population 2: Male (`"m"`), Female (`"f"`), or Total (`"t"`) assumed same as `sex1` unless otherwise specified.
#' @param closeout logical. Default `TRUE`. Shall we use the HMD Method Protocol to close out the `ax` and `qx` values? See details.
#' @details setting `closeout` to `TRUE` will result in value of `1/mx` for the final age group, of `ax` and a value of 1 for the closeout of `qx`. This function gives numerically identical results to `arriaga_sym()`, `lopez_ruzicka_sym()`, and `chandrasekaran_III()`.
#'
#' @return `cc` numeric vector with one element per age group, and which sums to the total difference in life expectancy between population 1 and 2.
#'
#' @importFrom data.table shift
#'
#' @export
#'
#' @references
#' \insertRef{arriaga1984measuring}{LEdecomp}
#' \insertRef{preston2000demography}{LEdecomp}
#' \insertRef{Ponnapalli2005}{LEdecomp}
#'
#' @examples
#' a <- .001
#' b <- .07
#' x <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' cc <- chandrasekaran_II(mx1, mx2, age = x)
#' e01 <- mx_to_e0(mx1, age = x)
#' e02 <- mx_to_e0(mx2, age = x)
#' (delta <- e02 - e01)
#' sum(cc)
#'
#'\dontrun{
#'  plot(x, cc)
#'}

chandrasekaran_II <- function(mx1, mx2,
                              age = (1:length(mx1))-1,
                              nx = rep(1, length(mx1)),
                              sex1 = 't',
                              sex2 = sex1,
                              closeout = TRUE){
  ax1 <- mx_to_ax(mx = mx1,
                  age = age,
                  nx = nx,
                  sex = sex1,
                  closeout = closeout)
  ax2 <- mx_to_ax(mx = mx2,
                  age = age,
                  nx = nx,
                  sex = sex2,
                  closeout = closeout)
  qx1 <- mx_to_qx(mx = mx1,
                  ax = ax1,
                  nx = nx,
                  closeout = closeout)
  qx2 <- mx_to_qx(mx = mx2,
                  ax = ax2,
                  nx = nx,
                  closeout = closeout)
  lx1 <- qx_to_lx(qx1)
  lx2 <- qx_to_lx(qx2)
  dx1 <- lx_to_dx(lx1)
  dx2 <- lx_to_dx(lx2)
  Lx1 <- ald_to_Lx(ax = ax1,
                   lx = lx1,
                   dx = dx1,
                   nx = nx)
  Lx2 <- ald_to_Lx(ax = ax2,
                   lx = lx2,
                   dx = dx2,
                   nx = nx)
  Tx1 <- rcumsum(Lx1)
  Tx2 <- rcumsum(Lx2)
  ex1 <- Tx1 / lx1
  ex2 <- Tx2 / lx2

  # from here on, everything uses lx and ex only;
  # we need to lag each of these, for easier code reading
  ex1_next <- shift(ex1, n = -1, fill = 0)
  ex2_next <- shift(ex2, n = -1, fill = 0)
  lx1_next <- shift(lx1, n = -1, fill = 0)
  lx2_next <- shift(lx2, n = -1, fill = 0)

  # eq 1.3
  interaction_effect_deferred <-
    lx2 * (ex2 - ex1) - lx2_next * (ex2_next - ex1_next)

  # eq 1.4
  interaction_effect_forwarded <-
    lx1 * (ex2 - ex1) - lx1_next * (ex2_next - ex1_next)

  # eq 1.5 (open age group already handled, since
  # *_next values close with 0s)
  approachII <- (interaction_effect_deferred + interaction_effect_forwarded)/2

  # TR: this line was redundant
  # approachII[length(ax1)] <- ((ex2[length(ax1)] - ex1[length(ax1)]) * (lx2[length(ax1)] + lx1[length(ax1)])) / 2

  approachII
}

#' @title Sensitivity from Chandrasekaran II decomposition
#'
#' @description Computes the sensitivity of life expectancy to changes in age-specific mortality rates using the Chandrasekaran II decomposition approach described by Ponnapalli (2005). The sensitivity is obtained by dividing the age-specific contributions (from `chandrasekaran_II()`) by the differences in mortality rates (`mx2 - mx1`). This yields a pointwise estimate of the derivative of life expectancy with respect to each age-specific mortality rate evaluated at an imagined midpoint between the first a second set of mortality rates.
#' @details This give numerically identical results to `sen_arriaga_sym()`, `sen_lopez_ruzicka_sym()`, and `sen_chandrasekaran_III()`.
#'
#' @inheritParams chandrasekaran_II
#'
#' @return A numeric vector of sensitivity values by age group.
#'
#' @seealso
#' \code{\link{chandrasekaran_II}}, \code{\link{sen_arriaga}}, \code{\link{sen_arriaga_sym}}
#'
#' @references
#' \insertRef{Ponnapalli2005}{LEdecomp}
#'
#' @export
#'
#' @examples
#' a <- 0.001
#' b <- 0.07
#' x <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' s <- sen_chandrasekaran_II(mx1, mx2, age = x)
#'
#' # Check that multiplying sensitivity by rate difference approximates the decomposition
#' cc_check <- s * (mx2 - mx1)
#' cc <- chandrasekaran_II(mx1, mx2, age = x)
#' \dontrun{
#' plot(x, cc, type = "l")
#' lines(x, cc_check, col = "red", lty = 2)
#' }

sen_chandrasekaran_II <- function(mx1,
                                  mx2,
                                  age = (1:length(mx1))-1,
                                  nx = rep(1,length(mx1)),
                                  sex1 = 't',
                                  sex2 = sex1,
                                  closeout = TRUE){
  # TR: why not call a different function to do all this?
  approachII <- chandrasekaran_II(mx1 = mx1,
                                  mx2 = mx2,
                                  age = age,
                                  nx = nx,
                                  sex1 = sex1,
                                  sex2 = sex2,
                                  closeout = closeout)
  delta <- mx2 - mx1
  sen <- approachII / delta
  sen

}

#' @title Instantaneous sensitivity via Chandrasekaran II decomposition
#'
#' @description
#' Estimates the sensitivity of life expectancy to small changes in mortality rates using the Chandrasekaran II decomposition. This is done by perturbing the input mortality rates up and down by a small factor and calculating the directional sensitivity.
#'
#' Specifically, the function constructs:
#' \deqn{m_{x}^{1} = m_x \cdot \left(\frac{1}{1 - h}\right)}
#' \deqn{m_{x}^{2} = m_x \cdot (1 - h)}
#' and applies \code{sen_chandrasekaran_II(mx1, mx2, ...)} to the result.
#'
#' @inheritParams sen_chandrasekaran_II
#' @param mx Numeric vector of mortality rates (central death rates).
#' @param sex Character; "m" for male, "f" for female, or "t" for total.
#' @param perturb Numeric; a small constant determining the perturbation size (default 1e-6).
#'
#' @details This approach gives a reasonable approximation of the derivative of life expectancy with respect to each age-specific mortality rate. It gives numerically identical results to `sen_arriaga_sym_instantaneous()`, `sen_lopez_ruzicka_instantaneous()`, and `sen_chandrasekaran_III_instantaneous()`.
#'
#' @seealso
#' \code{\link{sen_chandrasekaran_II}}, \code{\link{sen_chandrasekaran_II_instantaneous2}}, \code{\link{sen_arriaga_sym_instantaneous}}
#'
#' @export
#'
#' @examples
#' a <- 0.001
#' b <- 0.07
#' x <- 0:100
#' mx <- a * exp(x * b)
#' s <- sen_chandrasekaran_II_instantaneous(mx, age = x)
#' \dontrun{
#' plot(x, s, type = "l")
#' }

sen_chandrasekaran_II_instantaneous <- function(mx,
                                                age = (1:length(mx1))-1,
                                                nx = rep(1,length(mx)),
                                                sex = 't',
                                                perturb = 1e-6,
                                                closeout = TRUE){
  mx1 <- mx * (1 / (1 - perturb))
  mx2 <- mx * (1 - perturb) / 1
  s1 <- sen_chandrasekaran_II(mx1 = mx1,
                              mx2 = mx2,
                              age = age,
                              nx = nx,
                              sex1 = sex,
                              sex2 = sex,
                              closeout = closeout)
  s1
}
#' @title Log-space instantaneous sensitivity via Chandrasekaran II decomposition
#'
#' @description Estimates the sensitivity of life expectancy to small changes in mortality rates using the Chandrasekaran II decomposition. This variant perturbs the mortality rates in **log space**, creating two versions of `mx` by adding and subtracting a small constant to `log(mx)`, then exponentiating.
#'
#' Specifically:
#' \deqn{m_{x}^{1} = \exp\left( \ln m_x + h \right)}
#' \deqn{m_{x}^{2} = \exp\left( \ln m_x - h \right)}
#' and applies \code{sen_chandrasekaran_II(mx1, mx2, ...)} to the result.
#'
#' @inheritParams sen_chandrasekaran_II_instantaneous
#'
#' @details This approach provides a log-linear perturbation of the mortality schedule and can be used to estimate the derivative of life expectancy with respect to logged mortality rates. It gives numerically identical results to `sen_arriaga_sym_instantaneous2()`, `sen_lopez_ruzicka_instantaneous2()`, and `sen_chandrasekaran_III_instantaneous2()`.
#'
#' @seealso
#' \code{\link{sen_chandrasekaran_II_instantaneous}},
#' \code{\link{sen_arriaga_sym_instantaneous2}},
#' \code{\link{sen_lopez_ruzicka_instantaneous2}}
#'
#' @export
#'
#' @examples
#' a <- 0.001
#' b <- 0.07
#' x <- 0:100
#' mx <- a * exp(x * b)
#' s <- sen_chandrasekaran_II_instantaneous2(mx, age = x)
#' \dontrun{
#' plot(x, s, type = "l")
#' }

sen_chandrasekaran_II_instantaneous2 <- function(mx,
                                                 age = (1:length(mx1))-1,
                                                 nx = rep(1,length(mx)),
                                                 sex = 't',
                                                 perturb = 1e-6,
                                                 closeout = TRUE){
  mx1 <- exp(log(mx) + perturb)
  mx2 <- exp(log(mx) - perturb)
  mx1 <- mx * (1 / (1 - perturb))
  mx2 <- mx * (1 - perturb) / 1
  s1 <- sen_chandrasekaran_II(mx1 = mx1,
                              mx2 = mx2,
                              age = age,
                              nx = nx,
                              sex1 = sex,
                              sex2 = sex,
                              closeout = closeout)
  s1
}

#' @title Chandrasekaran III decomposition
#'
#' @description Implements the Chandrasekaran III decomposition as described in Ponnapalli (2005), which combines multiple directional effects into a symmetric average. The method constructs a decomposition of the difference in life expectancy into four parts: the main effect, the operative effect, their average (exclusive effect), and a non-linear interaction term. These are calculated based on life table values.
#' Let \eqn{e_x^i} denote remaining life expectancy at age \eqn{x} for population \eqn{i}, and \eqn{l_x^i} the number of survivors to age \eqn{x}. Then:
#'
#' - **Main effect**:
#' \deqn{
#' \frac{l_x^1}{l_x^2} \left[ l_x^2 (e_x^2 - e_x^1) - l_{x+n}^2 (e_{x+n}^2 - e_{x+n}^1) \right]
#' }
#'
#' - **Operative effect**:
#' \deqn{
#' \frac{l_x^2}{l_x^1} \left[ l_x^1 (e_x^2 - e_x^1) - l_{x+n}^1 (e_{x+n}^2 - e_{x+n}^1) \right]
#' }
#'
#' - **Exclusive effect**:
#' \deqn{
#' \frac{\text{Main effect} + \text{Operative effect}}{2}
#' }
#'
#' - **Interaction effect**:
#' \deqn{
#' (e_{x+n}^2 - e_{x+n}^1) \cdot \frac{1}{2} \left[
#' \frac{l_x^1 \cdot l_{x+n}^2}{l_x^2} + \frac{l_x^2 \cdot l_{x+n}^1}{l_x^1}
#' - (l_{x+n}^1 + l_{x+n}^2)
#' \right]
#' }
#'
#' The final contribution by age group is the sum of exclusive and interaction effects.
#'
#' @inheritParams chandrasekaran_II
#'
#' @details This decomposition gives numerically identical results to `arriaga_sym()`, `lopez_ruzicka_sym()`, and `chandrasekaran_II()`, despite conceptual differences in their derivation. Included here for methodological completeness.
#'
#' @return Numeric vector of contributions by age group that sum to the total difference
#' in life expectancy between the two mortality schedules.
#'
#' @references
#' \insertRef{Ponnapalli2005}{LEdecomp}
#'
#' @seealso
#' \code{\link{chandrasekaran_II}}, \code{\link{arriaga_sym}}, \code{\link{lopez_ruzicka_sym}}
#'
#' @export
#'
#' @examples
#' a <- 0.001
#' b <- 0.07
#' x <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' cc <- chandrasekaran_III(mx1, mx2, age = x)
#' e01 <- mx_to_e0(mx1, age = x)
#' e02 <- mx_to_e0(mx2, age = x)
#' (delta <- e02 - e01)
#' sum(cc)
#' \dontrun{
#' plot(x, cc, type = "l")
#' }

chandrasekaran_III <- function(mx1, mx2,
                               age = (1:length(mx1))-1,
                               nx = rep(1,length(mx1)),
                               sex1 = 't',
                               sex2 = sex1,
                               closeout = TRUE){
  ax1 <- mx_to_ax(mx = mx1,
                  age = age,
                  nx = nx,
                  sex = sex1,
                  closeout = closeout)
  ax2 <- mx_to_ax(mx = mx2,
                  age = age,
                  nx = nx,
                  sex = sex2,
                  closeout = closeout)
  qx1 <- mx_to_qx(mx = mx1,
                  ax = ax1,
                  nx = nx,
                  closeout = closeout)
  qx2 <- mx_to_qx(mx = mx2,
                  ax = ax2,
                  nx = nx,
                  closeout = closeout)
  lx1 <- qx_to_lx(qx1)
  lx2 <- qx_to_lx(qx2)
  dx1 <- lx_to_dx(lx1)
  dx2 <- lx_to_dx(lx2)
  Lx1 <- ald_to_Lx(ax = ax1,
                   lx = lx1,
                   dx = dx1,
                   nx = nx)
  Lx2 <- ald_to_Lx(ax = ax2,
                   lx = lx2,
                   dx = dx2,
                   nx = nx)
  Tx1 <- rcumsum(Lx1)
  Tx2 <- rcumsum(Lx2)
  ex1 <- Tx1 / lx1
  ex2 <- Tx2 / lx2

  # from here on, everything uses lx and ex only;
  # we need to lag each of these, for easier code reading
  ex1_next <- shift(ex1, n = -1, fill = 0)
  ex2_next <- shift(ex2, n = -1, fill = 0)
  lx1_next <- shift(lx1, n = -1, fill = 0)
  lx2_next <- shift(lx2, n = -1, fill = 0)

  # eq 1.1
  main_effect <- (lx1 / lx2) *
    (lx2 * (ex2 - ex1) -
       lx2_next * (ex2_next - ex1_next))

  # eq 1.2
  operative_effect <- (lx2 / lx1) *
    (lx1 * (ex2 - ex1) - lx1_next * (ex2_next - ex1_next))

  # eq 1.7 (avg of operative and main effects)
  # It must be equal to (main_effect + operative_effect)/2
  exclusive_effect <- (main_effect + operative_effect)/2

  # eq 1.8 (per Ponnapalli from Murthy & Gandhi (2004))
  interaction_effect <- (ex2_next - ex1_next) * (((lx1*lx2_next)/lx2 + (lx2*lx1_next)/lx1) -
                                                   (lx2_next + lx1_next))/2

  total_effect <- exclusive_effect + interaction_effect
  total_effect
}

#' @title Sensitivity from Chandrasekaran III decomposition
#'
#' @description Computes the implied sensitivity of life expectancy to changes in age-specific mortality rates using the Chandrasekaran III decomposition approach described by Ponnapalli (2005). The sensitivity is obtained by dividing the age-specific contributions (from `chandrasekaran_III()`) by the differences in mortality rates (`mx2 - mx1`). This yields a pointwise estimate of the derivative of life expectancy with respect to each age-specific mortality rate evaluated at an imagined midpoint between the first a second set of mortality rates.
#'
#' @details This gives numerically identical results to `sen_arriaga_sym()`, `sen_lopez_ruzicka_sym()`, and `sen_chandrasekaran_II()`.
#'
#' @inheritParams chandrasekaran_II
#'
#' @return A numeric vector of sensitivity values by age group.
#'
#' @seealso
#' \code{\link{chandrasekaran_III}}, \code{\link{sen_chandrasekaran_II}},
#' \code{\link{sen_arriaga_sym}}, \code{\link{sen_lopez_ruzicka_sym}}
#'
#' @references
#' \insertRef{Ponnapalli2005}{LEdecomp}
#'
#' @export
#'
#' @examples
#' a <- 0.001
#' b <- 0.07
#' x <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' s <- sen_chandrasekaran_III(mx1, mx2, age = x)
#'
#' # Check that multiplying sensitivity by rate difference approximates the decomposition
#' cc_check <- s * (mx2 - mx1)
#' cc <- chandrasekaran_III(mx1, mx2, age = x)
#' \dontrun{
#' plot(x, cc, type = "l")
#' lines(x, cc_check, col = "red", lty = 2)
#' }

sen_chandrasekaran_III <- function(mx1, mx2,
                                   age = (1:length(mx1)) - 1,
                                   nx = rep(1,length(mx1)),
                                   sex1 = 't',
                                   sex2 = sex1,
                                   closeout = TRUE){
  exclusive_effect <- chandrasekaran_III(mx1 = mx1,
                                         mx2 = mx2,
                                         age = age,
                                         nx = nx,
                                         sex1 = sex1,
                                         sex2 = sex2,
                                         closeout = closeout)

  delta <- mx2 - mx1
  sen <- exclusive_effect / delta
  sen

}

#' @title Instantaneous sensitivity via Chandrasekaran III decomposition
#'
#' @description Estimates the sensitivity of life expectancy to small changes in mortality rates using the Chandrasekaran III decomposition. This is done by perturbing the input mortality rates up and down by a small factor and computing directional sensitivity from the result.
#'
#' Specifically:
#' \deqn{m_x^{1} = m_x \cdot \left( \frac{1}{1 - h} \right)}
#' \deqn{m_x^{2} = m_x \cdot (1 - h)}
#' and applies \code{sen_chandrasekaran_III(mx1, mx2, ...)} to the result.
#'
#' @inheritParams sen_chandrasekaran_III
#' @param mx Numeric vector of mortality rates (central death rates).
#' @param sex Character; "m" for male, "f" for female, or "t" for total.
#' @param perturb Numeric; a small constant determining the perturbation size (default: 1e-6).
#'
#' @details
#' This approach provides an approximation of the derivative of life expectancy with
#' respect to each age-specific mortality rate, evaluated near the input `mx`.
#' It gives numerically identical results to `sen_arriaga_sym_instantaneous()`,
#' `sen_lopez_ruzicka_instantaneous()`, and `sen_chandrasekaran_II_instantaneous()`.
#'
#' @seealso
#' \code{\link{sen_chandrasekaran_III}},
#' \code{\link{sen_chandrasekaran_III_instantaneous2}},
#' \code{\link{sen_arriaga_sym_instantaneous}},
#' \code{\link{sen_lopez_ruzicka_instantaneous}}
#'
#' @export
#'
#' @examples
#' a <- 0.001
#' b <- 0.07
#' x <- 0:100
#' mx <- a * exp(x * b)
#' s <- sen_chandrasekaran_III_instantaneous(mx, age = x)
#' \dontrun{
#' plot(x, s, type = "l")
#' }

sen_chandrasekaran_III_instantaneous <- function(mx,
                                                 age = (1:length(mx1))-1,
                                                 nx = rep(1,length(mx1)),
                                                 sex = 't',
                                                 perturb = 1e-6,
                                                 closeout = TRUE){
  mx1 <- mx * (1 / (1 - perturb))
  mx2 <- mx * (1 - perturb) / 1
  s1 <- sen_chandrasekaran_III(mx1 = mx1,
                               mx2 = mx2,
                               age = age,
                               nx = nx,
                               sex1 = sex,
                               sex2 = sex,
                               closeout = closeout)
  s1
}


#' @title Log-space instantaneous sensitivity via Chandrasekaran III decomposition
#'
#' @description Estimates the sensitivity of life expectancy to small changes in mortality rates using the Chandrasekaran III decomposition and log-transformed perturbations. The method perturbs `mx` up and down in log space and averages the directional sensitivities to approximate the derivative.
#'
#' Specifically:
#' \deqn{m_x^{1} = \exp(\ln m_x + h)}
#' \deqn{m_x^{2} = \exp(\ln m_x - h)}
#' and applies \code{sen_chandrasekaran_III(mx1, mx2, ...)} and \code{sen_chandrasekaran_III(mx2, mx1, ...)},
#' returning their average.
#'
#' @inheritParams sen_chandrasekaran_III_instantaneous
#'
#' @details This version uses symmetric log-space perturbations. It gives numerically identical results to `sen_arriaga_sym_instantaneous2()`, `sen_lopez_ruzicka_instantaneous2()`, and `sen_chandrasekaran_II_instantaneous2()`.
#'
#' @seealso
#' \code{\link{sen_chandrasekaran_III_instantaneous}},
#' \code{\link{sen_arriaga_sym_instantaneous2}},
#' \code{\link{sen_lopez_ruzicka_instantaneous2}},
#' \code{\link{sen_chandrasekaran_II_instantaneous2}}
#'
#' @export
#'
#' @examples
#' a <- 0.001
#' b <- 0.07
#' x <- 0:100
#' mx <- a * exp(x * b)
#' s <- sen_chandrasekaran_III_instantaneous2(mx, age = x)
#' \dontrun{
#' plot(x, s, type = "l")
#' }


sen_chandrasekaran_III_instantaneous2 <- function(mx,
                                                  age = (1:length(mx1))-1,
                                                  nx = rep(1,length(mx1)),
                                                  sex = 't',
                                                  perturb = 1e-6,
                                                  closeout = TRUE){
  mx1 <- exp(log(mx) + perturb)
  mx2 <- exp(log(mx) - perturb)
  s1 <- sen_chandrasekaran_III(mx1 = mx1,
                               mx2 = mx2,
                               age = age,
                               nx = nx,
                               sex1 = sex,
                               sex2 = sex,
                               closeout = closeout)
  s1
}


