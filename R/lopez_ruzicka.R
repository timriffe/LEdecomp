#' @title Lopez-Ruzicka decomposition
#'
#' @description
#' Implements the decomposition of life expectancy proposed by Lopez and Ruzicka, as described in Ponnapalli (2005). This method expresses the difference in life expectancy between two mortality schedules in terms of an exclusive effect and an interaction effect, using life table quantities.
#'
#' Let \eqn{e_x^i} denote remaining life expectancy at age \eqn{x} for population \eqn{i}, and \eqn{l_x^i} the number of survivors to age \eqn{x}. Then:
#'
#' - **Exclusive effect**:
#' \deqn{
#' \frac{l_x^1}{l_x^2} \left[ l_x^2 (e_x^2 - e_x^1) - l_{x+n}^2 (e_{x+n}^2 - e_{x+n}^1) \right]
#' }
#'
#' - **Interaction effect**:
#' \deqn{
#' (e_{x+n}^2 - e_{x+n}^1) \cdot \left[
#' \frac{l_x^1 \cdot l_{x+n}^2}{l_x^2} - l_{x+n}^1
#' \right]
#' }
#'
#' The total contribution to life expectancy difference in age group \eqn{x} is the
#' sum of the exclusive and interaction effects.
#'
#' @inheritParams arriaga
#'
#' @details This method produces **numerically identical results** to `arriaga()`.
#'
#' @return Numeric vector of contributions by age group that sum to the total difference
#' in life expectancy between the two mortality schedules.
#'
#' @references
#' \insertRef{Ponnapalli2005}{LEdecomp}
#'
#' @seealso
#' \code{\link{arriaga}}, \code{\link{chandrasekaran_III}}, \code{\link{lopez_ruzicka_sym}}
#'
#' @export
#'
#' @examples
#' a <- 0.001
#' b <- 0.07
#' x <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' cc <- lopez_ruzicka(mx1, mx2, age = x)
#' sum(cc)

lopez_ruzicka <- function(mx1,
                          mx2,
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

  exclusive_effect <- (lx1/lx2)*(lx2*(ex2-ex1) - lx2_next*(ex2_next - ex1_next))

  interaction_effect <- (ex2_next - ex1_next)*(((lx1*lx2_next)/lx2) - lx1_next)

  decomp <- exclusive_effect + interaction_effect
  decomp

}
#' @title Symmetric Lopez-Ruzicka decomposition
#'
#' @description Implements a symmetric version of the Lopez-Ruzicka decomposition by averaging the results from the forward and reverse directions. That is, \code{lopez_ruzicka_sym(mx1, mx2)} returns
#' \deqn{
#' \frac{1}{2} \left( \text{lopez\_ruzicka}(mx1, mx2) - \text{lopez\_ruzicka}(mx2, mx1) \right)
#' }
#' This symmetric adjustment ensures that the decomposition is directionally neutral.
#'
#' @inheritParams lopez_ruzicka
#'
#' @details This symmetric version gives **numerically identical results** to `arriaga_sym()`, `chandrasekaran_II()`, and `chandrasekaran_III()`.
#'
#' @return Numeric vector of contributions by age group that sum to the total difference
#' in life expectancy between the two mortality schedules.
#'
#' @seealso
#' \code{\link{lopez_ruzicka}}, \code{\link{arriaga_sym}},
#' \code{\link{chandrasekaran_II}}, \code{\link{chandrasekaran_III}}
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
#' d <- lopez_ruzicka_sym(mx1, mx2, age = x)
#'
#' # compare to arriaga_sym()
#' d2 <- arriaga_sym(mx1, mx2, age = x)
#' all.equal(d, d2)

lopez_ruzicka_sym <- function(mx1,
                              mx2,
                              age = (1:length(mx1))-1,
                              nx = rep(1,length(mx1)),
                              sex1 = 't',
                              sex2 = sex1,
                              closeout = TRUE){
  a1 <- lopez_ruzicka(mx1,
                      mx2,
                      age = age,
                      nx = nx,
                      sex1 = sex1,
                      sex2 = sex2,
                      closeout = closeout)
  a2 <- lopez_ruzicka(mx2,
                      mx1,
                      age = age,
                      nx = nx,
                      sex1 = sex2,
                      sex2 = sex1,
                      closeout = closeout)

  a_avg <- (a1 - a2) / 2
  a_avg

}


#' @title Sensitivity from Lopez-Ruzicka decomposition
#'
#' @description
#' Computes the sensitivity of life expectancy to changes in age-specific mortality rates using the Lopez-Ruzicka decomposition approach. The sensitivity is calculated by dividing the age-specific contributions (from `lopez_ruzicka()`) by the differences in mortality rates (`mx2 - mx1`). This gives a pointwise estimate of the derivative of life expectancy with respect to each age-specific mortality rate, evaluated at an imagined midpoint between the two input rate schedules.
#'
#' @details
#' This method gives numerically identical results to `sen_arriaga()`.
#'
#' @inheritParams lopez_ruzicka
#'
#' @return A numeric vector of sensitivity values by age group.
#'
#' @seealso
#' \code{\link{lopez_ruzicka}}, \code{\link{sen_arriaga}}, \code{\link{sen_chandrasekaran_III}}
#'
#' @references
#' \insertRef{Ponnapalli2005}{LEdecomp}
#'
#' @examples
#' a <- 0.001
#' b <- 0.07
#' x <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a / 2 * exp(x * b)
#' s <- sen_lopez_ruzicka(mx1, mx2, age = x)
#'
#' # Check that multiplying sensitivity by rate difference reproduces the decomposition
#' cc_check <- s * (mx2 - mx1)
#' cc <- lopez_ruzicka(mx1, mx2, age = x)
#' \dontrun{
#' plot(x, cc, type = "l")
#' lines(x, cc_check, col = "red", lty = 2)
#' }
#'
#' @export

sen_lopez_ruzicka <- function(mx1,
                              mx2,
                              age = (1:length(mx1))-1,
                              nx = rep(1,length(mx1)),
                              sex1 = 't',
                              sex2 = sex1,
                              closeout = TRUE){
  # TR: replaces redundant code from before
  decomp <- lopez_ruzicka(mx1 = mx1,
                          mx2 = mx2,
                          age = age,
                          nx = nx,
                          sex1 = sex1,
                          sex2 = sex2,
                          closeout = closeout)
  delta <- mx2 - mx1
  sen <- decomp/delta
  sen

}

#' @title Sensitivity from symmetric Lopez-Ruzicka decomposition
#'
#' @description Computes the sensitivity of life expectancy to changes in age-specific mortality rates using the symmetric version of the Lopez-Ruzicka decomposition, as described by Ponnapalli (2005). The sensitivity is obtained by dividing the symmetric decomposition result by the differences in mortality rates (`mx2 - mx1`). This yields a pointwise estimate of the derivative of life expectancy with respect to each age-specific mortality rate evaluated at an imagined midpoint between the first and second set of mortality rates.
#'
#' @details This method gives numerically identical results to `sen_arriaga_sym()`, `sen_chandrasekaran_II()`, and `sen_chandrasekaran_III()`.
#'
#' @inheritParams lopez_ruzicka
#'
#' @return A numeric vector of sensitivity values by age group.
#'
#' @seealso
#' \code{\link{lopez_ruzicka_sym}},
#' \code{\link{sen_arriaga_sym}},
#' \code{\link{sen_chandrasekaran_II}},
#' \code{\link{sen_chandrasekaran_III}}
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
#' s <- sen_lopez_ruzicka_sym(mx1, mx2, age = x)
#'
#' # Check equivalence with symmetric Arriaga
#' s2 <- sen_arriaga_sym(mx1, mx2, age = x)
#' all.equal(s, s2)

sen_lopez_ruzicka_sym <- function(mx1, mx2,
                                  age = (1:length(mx1))-1,
                                  nx = rep(1,length(mx1)),
                                  sex1 = 't',
                                  sex2 = sex1,
                                  closeout = TRUE){
  delta <- mx2 - mx1
  a_avg <- lopez_ruzicka_sym(mx1 = mx1,
                             mx2 = mx2,
                             age = age,
                             nx = nx,
                             sex1 = sex1,
                             sex2 = sex2,
                             closeout = closeout)
  a_avg / delta
}

#' @title Instantaneous sensitivity via Lopez-Ruzicka decomposition
#'
#' @description
#' Estimates the sensitivity of life expectancy to small changes in mortality rates
#' using the Lopez-Ruzicka decomposition. This is done by perturbing the input
#' mortality rates up and down by a small factor and averaging the resulting directional
#' sensitivities to approximate a symmetric derivative.
#'
#' Specifically:
#' \deqn{m_x^{1} = m_x \cdot \left( \frac{1}{1 - h} \right)}
#' \deqn{m_x^{2} = m_x \cdot (1 - h)}
#' and applies \code{sen_lopez_ruzicka(mx1, mx2, ...)} and \code{sen_lopez_ruzicka(mx2, mx1, ...)},
#' returning their average.
#'
#' @inheritParams sen_lopez_ruzicka_sym
#' @param mx Numeric vector of mortality rates (central death rates).
#' @param sex Character; "m" for male, "f" for female, or "t" for total.
#' @param perturb Numeric; a small constant determining the perturbation size (default: 1e-6).
#'
#' @details
#' This method gives numerically identical results to
#' `sen_arriaga_sym_instantaneous()`,
#' `sen_chandrasekaran_II_instantaneous()`, and
#' `sen_chandrasekaran_III_instantaneous()`.
#'
#' @seealso
#' \code{\link{sen_lopez_ruzicka}}, \code{\link{sen_lopez_ruzicka_instantaneous2}}, \code{\link{sen_arriaga_sym_instantaneous}}, \code{\link{sen_chandrasekaran_II_instantaneous}}
#'
#' @export
#'
#' @examples
#' a <- 0.001
#' b <- 0.07
#' x <- 0:100
#' mx <- a * exp(x * b)
#' s <- sen_lopez_ruzicka_instantaneous(mx, age = x)
#' \dontrun{
#' plot(x, s, type = "l")
#' }

sen_lopez_ruzicka_instantaneous <- function(mx,
                                            age = (1:length(mx1))-1,
                                            nx = rep(1,length(mx1)),
                                            sex = 't',
                                            perturb = 1e-6,
                                            closeout = TRUE){
  mx1 <- mx * (1 / (1 - perturb))
  mx2 <- mx * (1 - perturb) / 1
  s1 <- sen_lopez_ruzicka(mx1 = mx1,
                          mx2 = mx2,
                          age = age,
                          nx = nx,
                          sex1 = sex,
                          sex2 = sex,
                          closeout = closeout)
  s2 <- sen_lopez_ruzicka(mx1 = mx2,
                          mx2 = mx1,
                          age = age,
                          nx = nx,
                          sex1 = sex,
                          sex2 = sex,
                          closeout = closeout)

  sen <- (s1 + s2) / 2

  sen
}

#' @title Log-space instantaneous sensitivity via Lopez-Ruzicka decomposition
#'
#' @description Estimates the sensitivity of life expectancy to small changes in mortality rates using the Lopez-Ruzicka decomposition and log-space perturbation. This is done by shifting the log of the input mortality rates up and down by a small constant, then exponentiating, and computing the average directional sensitivity.
#'
#' Specifically:
#' \deqn{m_x^{1} = \exp(\ln m_x + h)}
#' \deqn{m_x^{2} = \exp(\ln m_x - h)}
#' and applies \code{sen_lopez_ruzicka(mx1, mx2, ...)} and \code{sen_lopez_ruzicka(mx2, mx1, ...)},
#' returning their average.
#'
#' @inheritParams sen_lopez_ruzicka_instantaneous
#'
#' @details
#' This approach gives numerically identical results to
#' `sen_arriaga_sym_instantaneous2()`,
#' `sen_chandrasekaran_II_instantaneous2()`, and
#' `sen_chandrasekaran_III_instantaneous2()`.
#'
#' @seealso
#' \code{\link{sen_lopez_ruzicka_instantaneous}},
#' \code{\link{sen_arriaga_sym_instantaneous2}},
#' \code{\link{sen_chandrasekaran_III_instantaneous2}}
#'
#' @export
#'
#' @examples
#' a <- 0.001
#' b <- 0.07
#' x <- 0:100
#' mx <- a * exp(x * b)
#' s <- sen_lopez_ruzicka_instantaneous2(mx, age = x)
#' \dontrun{
#' plot(x, s, type = "l")
#' }

sen_lopez_ruzicka_instantaneous2 <- function(mx,
                                             age = (1:length(mx1))-1,
                                             nx = rep(1,length(mx1)),
                                             sex = 't',
                                             perturb = 1e-6,
                                             closeout = TRUE){
  mx1 <- exp(log(mx) + perturb)
  mx2 <- exp(log(mx) - perturb)
  s1 <- sen_lopez_ruzicka(mx1 = mx1,
                          mx2 = mx2,
                          age = age,
                          nx = nx,
                          sex1 = sex,
                          sex2 = sex,
                          closeout = closeout)
  s2 <- sen_lopez_ruzicka(mx1 = mx2,
                          mx2 = mx1,
                          age = age,
                          nx = nx,
                          sex1 = sex,
                          sex2 = sex,
                          closeout = closeout)
  # TR: is this closeout actually needed here?
  # I have my doubts this was checked
  #To match the solution with the arriaga
  # if (closeout){
  #   s2[length(s2)] <- s2[length(s2)] * 2
  # }
  sen <- (s1 + s2) / 2

  sen
}
