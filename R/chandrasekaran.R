# This script contains functions that implement and adapt the
# Chandrasekaran decompositon approach

#' @title II approach of Chandrasekaran decomposition approach
#' @description Following the notation given in Ponnapalli (2005), and the decomposition method can written as:
#' \deqn{_{n}\Delta_{x} = \frac{\left(e_x^2 - e_x^1\right)\left(l_x^2 + l_x^1 \right)}{2} - \frac{\left(e^{2}_{x+n} - e^{1}_{x+n} \right) \left(l^{2}_{x+n} + l^{1}_{x+n} \right)}{2} - \frac{_{n}L_{x}^{1}}{l_{x}^{1}}\right) + \frac{T^{2}_{x+n}}{l_{0}^{1}} \cdot \left( \frac{l_{x}^{1}}{l_{x}^{2}} - \frac{l_{x+n}^{1}}{l_{x+n}^{2}}  \right) }
#' where \eqn{_{n}\Delta_{x}} is the contribution of rate differences in age \eqn{x} to the difference in life expectancy implied by `mx1` and `mx2`. This formula can be averaged between ‘effect interaction deferred’ and ‘effect interaction forwarded’ from the Ponnapalli (2005).
#'
#' @param mx1 numeric vector of the mortality rates (central death rates) for population 1
#' @param mx2 numeric vector of the mortality rates (central death rates) for population 2
#' @param age integer vector of the lower bound of each age group (currently only single ages supported)
#' @param sex1 character either the sex for population 1: Male (`"m"`), Female (`"f"`), or Total (`"t"`)
#' @param sex2 character either the sex for population 2: Male (`"m"`), Female (`"f"`), or Total (`"t"`) assumed same as `sex1` unless otherwise specified.
#' @param closeout logical. Default `TRUE`. Shall we use the HMD Method Protocol to close out the `ax` and `qx` values? See details.
#' @details setting `closeout` to `TRUE` will result in value of `1/mx` for the final age group, of `ax` and a value of 1 for the closeout of `qx`.
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
#' cc <- arriaga(mx1, mx2, age = x)
#' e01 <- mx_to_e0(mx1, age = x)
#' e02 <- mx_to_e0(mx2, age = x)
#' (delta <- e02 - e01)
#' sum(cc)
#'
#'\dontrun{
#'  plot(x, cc)
#'}
#' # asymmetrical with a decomposition in the opposite direction
#' cc2 <- -arriaga(mx1 = mx2, mx2 = mx1, age = x)
#' plot(x, cc)
#' lines(x,cc2)
#' # also we lose some precision?
#' sum(cc);sum(cc2)
#' # found it!
#' delta-sum(cc2); cc2[length(cc2)] / 2
#'
#' # But this is no problem if closeout = FALSE
#' -arriaga(mx1 = mx2, mx2 = mx1, age = x,closeout=FALSE) |> sum()
#' arriaga(mx1 = mx1, mx2 = mx2, age = x,closeout=FALSE) |> sum()
#  a <- .001
#  b <- .07
#  x <- 0:100
#  mx1 <- a * exp(x * b)
#  mx2 <- a/2 * exp(x * b)
#  sex1 = "m"
#  sex2 = "m"
#  closeout = TRUE
chandrasekaran_II <- function(mx1, mx2,
                              age,
                              sex1 = 't',
                              sex2 = sex1,
                              closeout = TRUE){
  ax1 <- mx_to_ax(mx = mx1,
                  age = age,
                  sex = sex1,
                  closeout = closeout)
  ax2 <- mx_to_ax(mx = mx2,
                  age = age,
                  sex = sex2,
                  closeout = closeout)
  qx1 <- mx_to_qx(mx = mx1,
                  ax = ax1,
                  closeout = closeout)
  qx2 <- mx_to_qx(mx = mx2,
                  ax = ax2,
                  closeout = closeout)
  lx1 <- qx_to_lx(qx1)
  lx2 <- qx_to_lx(qx2)
  dx1 <- lx_to_dx(lx1)
  dx2 <- lx_to_dx(lx2)
  Lx1 <- ald_to_Lx(ax = ax1,
                   lx = lx1,
                   dx = dx1)
  Lx2 <- ald_to_Lx(ax = ax2,
                   lx = lx2,
                   dx = dx2)
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

  approachII[length(ax1)] <- ((ex2[length(ax1)] - ex1[length(ax1)]) * (lx2[length(ax1)] + lx1[length(ax1)])) / 2

  approachII
}

sen_chandrasekaran_II_sym <- function(mx1,
                                      mx2,
                                      age = 0:(length(mx1) - 1),
                                      sex1 = 't',
                                      sex2 = sex1,
                                      closeout = TRUE){
  delta <- mx2 - mx1
  a_avg <- chandrasekaran_II_sym(mx1 = mx1,
                                 mx2 = mx2,
                                 age = age,
                                 sex1 = sex1,
                                 sex2 = sex2,
                                 closeout = closeout)
  a_avg / delta
}
sen_chandrasekaran_II <- function(mx1,
                                  mx2,
                                  age,
                                  sex1 = 't',
                                  sex2 = sex1,
                                  closeout = TRUE){
  # TR: why not call a different function to do all this?
  approachII <- chandrasekaran_II(mx1 = mx1,
                                  mx2 = mx2,
                                  age = age,
                                  sex1 = sex1,
                                  sex2 = sex2,
                                  closeout = closeout)
  delta <- mx2 - mx1
  sen <- approachII / delta
  sen

}
sen_chandrasekaran_II_instantaneous <- function(mx,
                                                age = 0:(length(mx1)-1),
                                                sex = 't',
                                                perturb = 1e-6,
                                                closeout = TRUE){
  mx1 <- mx * (1 / (1 - perturb))
  mx2 <- mx * (1 - perturb) / 1
  s1 <- sen_chandrasekaran_II(mx1 = mx1, mx2 = mx2,
                              age = age,
                              sex1 = sex, sex2 = sex, closeout = closeout)
  s1
}
chandrasekaran_II_sym <- function(mx1,
                                  mx2,
                                  age,
                                  sex1 = 't',
                                  sex2 = sex1,
                                  closeout = TRUE){
  a1 <- chandrasekaran_II(mx1,
                          mx2,
                          age = age,
                          sex1 = sex1,
                          sex2 = sex2,
                          closeout = closeout)
  a2 <- chandrasekaran_II(mx2,
                          mx1,
                          age = age,
                          sex1 = sex2,
                          sex2 = sex1,
                          closeout = closeout)

  a_avg <- (a1 - a2) / 2
  a_avg
}

sen_chandrasekaran_II_instantaneous2 <- function(mx,
                                                 age = 0:(length(mx1)-1),
                                                 sex = 't',
                                                 perturb = 1e-6,
                                                 closeout = TRUE){
  mx1 <- exp(log(mx) + perturb)
  mx2 <- exp(log(mx) - perturb)
  s1 <- sen_chandrasekaran_II(mx1 = mx1, mx2 = mx2,
                              age = age,
                              sex1 = sex, sex2 = sex, closeout = closeout)
  s1
}
chandrasekaran_III <- function(mx1, mx2,
                               age,
                               sex1 = 't',
                               sex2 = sex1,
                               closeout = TRUE){
  ax1 <- mx_to_ax(mx = mx1,
                  age = age,
                  sex = sex1,
                  closeout = closeout)
  ax2 <- mx_to_ax(mx = mx2,
                  age = age,
                  sex = sex2,
                  closeout = closeout)
  qx1 <- mx_to_qx(mx = mx1,
                  ax = ax1,
                  closeout = closeout)
  qx2 <- mx_to_qx(mx = mx2,
                  ax = ax2,
                  closeout = closeout)
  lx1 <- qx_to_lx(qx1)
  lx2 <- qx_to_lx(qx2)
  dx1 <- lx_to_dx(lx1)
  dx2 <- lx_to_dx(lx2)
  Lx1 <- ald_to_Lx(ax = ax1,
                   lx = lx1,
                   dx = dx1)
  Lx2 <- ald_to_Lx(ax = ax2,
                   lx = lx2,
                   dx = dx2)
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

  exclusive_effect

}

chandrasekaran_III_sym <- function(mx1,
                                  mx2,
                                  age,
                                  sex1 = 't',
                                  sex2 = sex1,
                                  closeout = TRUE){
  a1 <- chandrasekaran_III(mx1,
                           mx2,
                           age = age,
                           sex1 = sex1,
                           sex2 = sex2,
                           closeout = closeout)
  a2 <- chandrasekaran_III(mx2,
                           mx1,
                           age = age,
                           sex1 = sex2,
                           sex2 = sex1,
                           closeout = closeout)

  a_avg <- (a1 - a2) / 2
  a_avg
}
sen_chandrasekaran_III <- function(mx1, mx2,
                                   age,
                                   sex1 = 't',
                                   sex2 = sex1,
                                   closeout = TRUE){
  # TR: redundant code replaced with function call...
  exclusive_effect <- chandrasekaran_III(mx1 = mx1,
                                         mx2 = mx2,
                                         age = age,
                                         sex1 = sex1,
                                         sex2 = sex2,
                                         closeout = closeout)

  delta <- mx2 - mx1
  sen <- exclusive_effect / delta
  sen

}
sen_chandrasekaran_III_sym <- function(mx1,
                                       mx2,
                                       age = 0:(length(mx1) - 1),
                                       sex1 = 't',
                                       sex2 = sex1,
                                       closeout = TRUE){
  delta <- mx2 - mx1
  a_avg <- chandrasekaran_III_sym(mx1 = mx1,
                                  mx2 = mx2,
                                  age = age,
                                  sex1 = sex1,
                                  sex2 = sex2,
                                  closeout = closeout)
  a_avg / delta
}
sen_chandrasekaran_III_instantaneous <- function(mx,
                                                 age = 0:(length(mx1)-1),
                                                 sex = 't',
                                                 perturb = 1e-6,
                                                 closeout = TRUE){
  mx1 <- mx * (1 / (1 - perturb))
  mx2 <- mx * (1 - perturb) / 1
  s1 <- sen_chandrasekaran_III(mx1 = mx1, mx2 = mx2,
                               age = age,
                               sex1 = sex, sex2 = sex, closeout = closeout)
  s1
}
sen_chandrasekaran_III_instantaneous2 <- function(mx,
                                                  age = 0:(length(mx1)-1),
                                                  sex = 't',
                                                  perturb = 1e-6,
                                                  closeout = TRUE){
  mx1 <- exp(log(mx) + perturb)
  mx2 <- exp(log(mx) - perturb)
  s1 <- sen_chandrasekaran_III(mx1 = mx1, mx2 = mx2,
                               age = age,
                               sex1 = sex, sex2 = sex, closeout = closeout)
  s1
}


