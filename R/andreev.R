#' Andreev (1982) component decomposition of e0 difference by age
#'
#' Implements the age-interval component formula (telescoping form)
#' from Andreev (1982): eps[x1,x2) = l_x1*(e2_x1 - e1_x1) - l_x2*(e2_x2 - e1_x2),
#' where l_ and e_ without prime are from a chosen baseline life table.
#'
#' @param mx1 numeric vector of the mortality rates (central death rates) for population 1
#' @param mx2 numeric vector of the mortality rates (central death rates) for population 2
#' @param age integer vector of the lower bound of each age group (currently only single ages supported)
#' @param nx integer vector of age intervals, default 1.
#' @param sex1 character either the sex for population 1: Male (`"m"`), Female (`"f"`), or Total (`"t"`)
#' @param sex2 character either the sex for population 2: Male (`"m"`), Female (`"f"`), or Total (`"t"`) assumed same as `sex1` unless otherwise specified.
#' @param closeout logical. Default `TRUE`. Shall we use the HMD Method Protocol to close out the `ax` and `qx` values? See details.
#' @details setting `closeout` to `TRUE` will result in value of `1/mx` for the final age group, of `ax` and a value of 1 for the closeout of `qx`.
#' @return `cc` numeric vector with one element per age group, and which sums to the total difference in life expectancy between population 1 and 2.
#' @return vector of age-specific contributions to e0 gap
#' @references
#'   \insertRef{andreev1982}{LEdecomp}
#' @export
#' @examples
#' a   <- .001
#' b   <- .07
#' x   <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' cc <- andreev(mx1, mx2, age = x)
#'
#' e01 <- mx_to_e0(mx1, age = x)
#' e02 <- mx_to_e0(mx2, age = x)
#' (delta <- e02 - e01)
#'sum(cc)
andreev <- function(mx1,
                    mx2,
                    age = 0:(length(mx1) - 1),
                    nx = rep(1,length(mx1)),
                    sex1 = "t",
                    sex2 = sex1,
                    closeout = TRUE,
                    radix = 1) {

  stopifnot(is.numeric(mx1), is.numeric(mx2), is.numeric(age))
  stopifnot(length(mx1) == length(mx2), length(mx1) == length(age))
  if (any(mx1 < 0, na.rm = TRUE) || any(mx2 < 0, na.rm = TRUE)) {
    stop("mx must be nonnegative.")
  }

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

  # delta ex at age boundaries
  dex_21 <- ex2 - ex1

  lx1_next  <- c(lx1[-1], 0)
  dex_next  <- c(dex_21[-1], 0)
  cc <- lx1 * dex_21 - lx1_next * dex_next

  cc
}

#' @title the sensitivity implied by an Andreev decomposition
#' @description The sensitivity of life expectancy to a perturbation in mortality rates can be derived by dividing the Andreev decomposition result \eqn{\Delta} by the difference `mx2-mx1`.
#' \deqn{s_{x} = \frac{\Delta}{_{n}M^{2}_x - _{n}M^{1}_x}}

#' @seealso \code{\link{andreev}}
#' @inheritParams andreev
#' @return `s` numeric vector with one element per age group, and which gives the sensitivity values for each age.
#' @importFrom data.table shift
#' @export
#' @references
#'   \insertRef{andreev1982}{LEdecomp}
#' @examples
#' a <- .001
#' b <- .07
#' x <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' cc <- andreev(mx1, mx2, age = x)
#' # examples can come from above too
#' s <- sen_andreev(mx1, mx2, age = x)
#' \donttest{
#' plot(x, s)
#' }
#' cc_check <- s * (mx2 - mx1)
#' \donttest{
#' plot(x,cc)
#' lines(x,cc_check)
#' }

sen_andreev <- function(mx1,
                        mx2,
                        age = 0:(length(mx1)-1),
                        nx = rep(1,length(mx1)),
                        sex1 = 't',
                        sex2 = sex1,
                        closeout = TRUE){

  cc <- andreev(mx1 = mx1,
                mx2 = mx2,
                age = age,
                nx = nx,
                sex1 = sex1,
                sex2 = sex2,
                closeout = closeout)
  delta <- mx2 - mx1
  # watch out for 0s in delta denominator
  sen <- cc / delta
  sen
}

#' @title Estimate sensitivity of life expectancy for a set of mortality rates
#' @description This implementation gives a good approximation of the sensitivity of life expectancy to perturbations in mortality rates (central death rates). Since the Andreev approach requires two versions of mortality rates `mx1`, `mx2`, we create these by slightly perturbing `mx` up and down. Then we calculate the Andreev sensitivity in each direction and take the average. Specifically, we create `mx1` and `mx2` as
#' \deqn{m_{x}^{1} = m_x \cdot \left(\frac{1}{1 - h}\right)}
#' \deqn{m_{x}^{2} = m_x \cdot \left(1 - h\right)}
#' where `h` is small value given by the argument `perturb`.
#' @details A minor correction might be needed for the final age group for the case of the reverse-direction Andreev sensitivity. Note also for values of `perturb` (h) that are less than `1e-7` we might lose stability in results.
#' @inheritParams andreev
#' @param mx numeric vector of mortality rates (central death rates)
#' @param sex character Male (`"m"`), Female (`"f"`), or Total (`"t"`)
#' @param perturb numeric constant, a very small number
#' @importFrom data.table shift
#' @return numeric vector of sensitivity of life expectancy to perturbations in `mx`.
#' @export
#' @examples
#' a   <- .001
#' b   <- .07
#' x   <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' mx  <- (mx1 + mx2) / 2
#' s     <- sen_andreev_instantaneous(mx, age = x)
#' s1    <- sen_andreev_instantaneous(mx1, age = x)
#' s2    <- sen_andreev_instantaneous(mx2, age = x)
#' s1_d  <- sen_andreev(mx1, mx2, age = x)
#' s2_d  <- sen_andreev(mx2, mx1, age = x)
#' delta <- mx2 - mx1
#' # dots give our point estimate of sensitivity at the average of the rates,
#' # which is different from the
#'
#' plot(x,s*delta, ylim = c(0,.3))
#' lines(x,s1*delta,col = "red")
#' lines(x,s2*delta,col = "blue")
#' # the sensitivity of the average is different
#' # from the average of the sensitivities!
#' lines(x, ((s1+s2)) / 2 * delta)
#' # and these are different from the directional sensitivities
#' # covering the whole space from mx1 to mx2:
#' lines(x, s1_d*delta, col = "red", lty =2)
#' lines(x, s2_d*delta, col = "blue", lty =2)
sen_andreev_instantaneous <- function(mx,
                                      age = 0:(length(mx1)-1),
                                      sex = 't',
                                      nx = rep(1,length(mx)),
                                      perturb = 1e-6,
                                      closeout = TRUE){
  mx1 <- mx * (1 / (1 - perturb))
  mx2 <- mx * (1 - perturb) / 1
  s1 <- sen_andreev(mx1 = mx1,
                    mx2 = mx2,
                    age = age,
                    nx = nx,
                    sex1 = sex,
                    sex2 = sex,
                    closeout = closeout)
  s2 <- sen_andreev(mx1 = mx2,
                    mx2 = mx1,
                    age = age,
                    nx = nx,
                    sex1 = sex,
                    sex2 = sex,
                    closeout = closeout)

  (s1 + s2) / 2
}

#' @title Estimate sensitivity of life expectancy for a set of mortality rates by perturbing in the log space.
#' @description This is a second approach for estimating the sensitivity for a single set of rates. Here, rather than directly expanding and contracting rates to convert `mx` into `mx1` and `mx2` we instead shift the logged mortality rates up and down by the factor `perturb = h`. Specifically:
#' \deqn{m_{x}^{1} = e^{\ln\left(m_x\right) + h}}
#' \deqn{m_{x}^{2} = e^{\ln\left(m_x\right) - h}}
#' @export
#' @return numeric vector of sensitivity of life expectancy to perturbations in `mx`
#' @inheritParams sen_arriaga_instantaneous
#' @seealso \code{\link{sen_andreev_instantaneous}}
#' @examples
#' a   <- .001
#' b   <- .07
#' x   <- 0:100
#' mx <- a * exp(x * b)
#' # the multiplicative perturbation:
#' s1 <- sen_andreev_instantaneous(mx)
#' s2 <- sen_andreev_instantaneous2(mx)
#' plot(x,
#'      s1 - s2,
#'      pch = 16,
#'      cex=.5,
#'      main = "very similar")
#'
sen_andreev_instantaneous2 <- function(mx,
                                       age = 0:(length(mx1)-1),
                                       sex = 't',
                                       nx = rep(1,length(mx)),
                                       perturb = 1e-6,
                                       closeout = TRUE){
  mx1 <- exp(log(mx) + perturb)
  mx2 <- exp(log(mx) - perturb)
  s1 <- sen_andreev(mx1 = mx1,
                    mx2 = mx2,
                    age = age,
                    nx = nx,
                    sex1 = sex,
                    sex2 = sex,
                    closeout = closeout)
  s2 <- sen_andreev(mx1 = mx2,
                    mx2 = mx1,
                    age = age,
                    nx = nx,
                    sex1 = sex,
                    sex2 = sex,
                    closeout = closeout)

  (s1 + s2) / 2
}

#' @title Decompose using a symmetrical variant of the Andreev components approach.
#' @description This approach conducts a classic Andreev (1982) decomposition in both directions, averaging the (sign-adjusted) result, i.e. `a_avg = (andreev(mx1,mx2, ...) - andreev(mx2, mx1, ...)) / 2`.
#' #@note The final age group's contribution from the reversed decomposition is halved before averaging. This empirical adjustment ensures symmetry and numeric stability, though the theoretical basis requires further exploration.
#' @return numeric vector of contributions summing to the gap in life expectancy implied by `mx1` and `mx2`.
#' @export
#' @inheritParams andreev
#' @seealso \code{\link{andreev}}
#' @examples
#' a <- .001
#' b <- .07
#' x <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' d <- andreev_sym(mx1, mx2, age = x)
#'
#' e01 <- mx_to_e0(mx1,age=x)
#' e02 <- mx_to_e0(mx2,age=x)
#' (Delta <- e02 - e01)
#' sum(d)
#'
#' d12 <- andreev(mx1, mx2, age = x)
#' d21 <- andreev(mx2, mx1, age = x) # direction opposite
#' \donttest{
#' plot(x, d, type= 'l')
#'   lines(x, d12, col = "red")
#'   lines(x, -d21, col = "blue")
#' }
andreev_sym <- function(mx1,
                        mx2,
                        age = 0:(length(mx1) - 1),
                        nx = rep(1,length(mx1)),
                        sex1 = 't',
                        sex2 = sex1,
                        closeout = TRUE){
  a1 <- andreev(mx1,
                mx2,
                age = age,
                nx = nx,
                sex1 = sex1,
                sex2 = sex2,
                closeout = closeout)
  a2 <- andreev(mx2,
                mx1,
                age = age,
                nx = nx,
                sex1 = sex2,
                sex2 = sex1,
                closeout = closeout)

  a_avg <- (a1 - a2) / 2
  a_avg
}

#' @title Estimate sensitivity of life expectancy using a symmetrical Andreev approach.
#' @description This approach conducts a classic Andreev (1982) decomposition in both directions, averaging the (sign-adjusted) result, i.e. `a_avg = (andreev(mx1,mx2, ...) - andreev(mx2, mx1, ...)) / 2`, then approximates the sensitivity by dividing out the rate differences, i.e. `s = a_avg / (mx2 - mx1)`. A resulting decomposition will be exact because the two andreev directions are exact, but this method might be vulnerable to 0s in the denominator.
#' @return numeric vector of life expectancy sensitivity to perturbations in mx evaluated at the average of `mx1` and `mx2`.
#' @export
#' @inheritParams arriaga
#' @seealso \code{\link{arriaga}}
#' @examples
#' a <- .001
#' b <- .07
#' x <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' s <- sen_andreev_sym(mx1, mx2, age = x)
#'
#' e01 <- mx_to_e0(mx1,age=x)
#' e02 <- mx_to_e0(mx2,age=x)
#' (Delta <- e02 - e01)
#' deltas <- mx2- mx1
#' sum(deltas * s)
#' mx_avg <- (mx1+mx2) / 2
#' \donttest{
#' mx_avg <- (mx1 + mx2) / 2
#' plot(x, s, type = 'l')
#' lines(x, sen_andreev_instantaneous(mx_avg, age=x),col = "blue")
#' }

sen_andreev_sym <- function(mx1,
                            mx2,
                            age = 0:(length(mx1) - 1),
                            nx = rep(1,length(mx1)),
                            sex1 = 't',
                            sex2 = sex1,
                            closeout = TRUE){
  delta <- mx2 - mx1
  a_avg <- andreev_sym(mx1 = mx1,
                       mx2 = mx2,
                       age = age,
                       nx = nx,
                       sex1 = sex1,
                       sex2 = sex2,
                       closeout = closeout)
  a_avg / delta
}

#' @title Instantaneous sensitivity via symmetrical Andreev decomposition
#'
#' @description
#' Estimates the sensitivity of life expectancy to small changes in age-specific mortality rates using the symmetrical Andreev decomposition. This is done by applying a small multiplicative perturbation to the input mortality rates and using the symmetrical sensitivity function `sen_andreev_sym()`.
#'
#' Specifically, the function constructs:
#' \deqn{m_{x}^{1} = m_x \cdot \left(\frac{1}{1 - h}\right)}
#' \deqn{m_{x}^{2} = m_x \cdot (1 - h)}
#' and applies \code{sen_andreev_sym(mx1, mx2, ...)} to the result.
#'
#' @inheritParams sen_andreev_sym
#' @param mx Numeric vector of mortality rates (central death rates).
#' @param sex Character; "m" for male, "f" for female, or "t" for total.
#' @param perturb Numeric; a small constant determining the perturbation size (default 1e-6).
#'
#' @details This function yields an instantaneous approximation to the derivative of life expectancy with respect to mortality, evaluated at the input schedule. Because `sen_andreev_sym()` is itself symmetrical, only the "forward" perturbation is required.
#' @return numeric vector of life expectancy sensitivity to perturbations in `mx.`
#' @seealso \code{\link{sen_andreev_sym}}, \code{\link{sen_andreev_sym_instantaneous2}}, \code{\link{sen_lopez_ruzicka_instantaneous}}
#'
#' @export
#' @examples
#' a <- 0.001
#' b <- 0.07
#' x <- 0:100
#' mx <- a * exp(x * b)
#' s <- sen_andreev_sym_instantaneous(mx, age = x)
#' \donttest{
#' plot(x, s, type = "l")
#' }


sen_andreev_sym_instantaneous <- function(mx,
                                          age = 0:(length(mx1)-1),
                                          sex = 't',
                                          nx = rep(1,length(mx)),
                                          perturb = 1e-6,
                                          closeout = TRUE){
  mx1 <- mx * (1 / (1 - perturb))
  mx2 <- mx * (1 - perturb) / 1
  s1 <- sen_andreev_sym(mx1 = mx1,
                        mx2 = mx2,
                        age = age,
                        nx = nx,
                        sex1 = sex,
                        sex2 = sex,
                        closeout = closeout)

  s1
}

#' @title Estimate sensitivity of life expectancy for a set of mortality rates by perturbing in the log space.
#' @description This is a second approach for estimating the sensitivity for a single set of rates. Here, rather than directly expanding and contracting rates to convert `mx` into `mx1` and `mx2` we instead shift the logged mortality rates up and down by the factor `perturb = h`. Specifically:
#' \deqn{m_{x}^{1} = e^{\ln\left(m_x\right) + h}}
#' \deqn{m_{x}^{2} = e^{\ln\left(m_x\right) - h}}
#' @export
#' @inheritParams sen_andreev_instantaneous
#' @seealso \code{\link{sen_andreev_instantaneous}}
#' @return numeric vector of life expectancy sensitivity to perturbations in `mx.`
#' @examples
#' a <- 0.001
#' b <- 0.07
#' x <- 0:100
#' mx <- a * exp(x * b)
#' s <- sen_andreev_sym_instantaneous2(mx, age = x)
#' \donttest{
#' plot(x, s, type = "l")
#' }

sen_andreev_sym_instantaneous2 <- function(mx,
                                           age = 0:(length(mx1)-1),
                                           sex = 't',
                                           nx = rep(1,length(mx)),
                                           perturb = 1e-6,
                                           closeout = TRUE){
  mx1 <- exp(log(mx) + perturb)
  mx2 <- exp(log(mx) - perturb)
  s1 <- sen_andreev_sym(mx1 = mx1,
                        mx2 = mx2,
                        age = age,
                        nx = nx,
                        sex1 = sex,
                        sex2 = sex,
                        closeout = closeout)
  s1
}
