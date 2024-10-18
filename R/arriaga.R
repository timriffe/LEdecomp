# This script contains functions that implement and adapt the 
# Arriaga decompositon approach

#' @title classic Arriaga decomposition
#' @description Following the notation given in Preston et al (2000), Arriaga's decomposition method can written as:
#' \deqn{_{n}\Delta_{x} = \frac{l_x^1}{l_0^1}\cdot \left( \frac{_{n}L_{x}^{2}}{l_{x}^{2}} - \frac{_{n}L_{x}^{1}}{l_{x}^{1}}\right) + \frac{T^{2}_{x+n}}{l_{0}^{1}} \cdot \left( \frac{l_{x}^{1}}{l_{x}^{2}} - \frac{l_{x+n}^{1}}{l_{x+n}^{2}}  \right) }
#' where \eqn{_{n}\Delta_{x}} is the contribution of rate differences in age \eqn{x} to the difference in life expectancy implied by `mx1` and `mx2`. The first part of the sum refers to the so-called direct effect (in age group \eqn{x}), whereas the second part refers to indirect effects in ages above age \eqn{x}. Usually the indirect effects are far larger than the direct effects.
#' @details A little-known property of this decomposition method is that it is directional, in the sense that we are comparing a movement of `mx1` to `mx2`, and this is not exactly symmetrical with a comparison of `mx2` with `mx1`. Note also, if decomposing in reverse from the usual, you may need a slight adjustment to the closeout value in order to match sums properly (see examples for a demonstration).
#' 
#' @param mx1 numeric vector of the mortality rates (central death rates) for population 1
#' @param mx2 numeric vector of the mortality rates (central death rates) for population 2
#' @param age integer vector of the lower bound of each age group (currently only single ages supported)
#' @param sex1 character either the sex for population 1: Male (`"m"`), Female (`"f"`), or Total (`"t"`)
#' @param sex2 character either the sex for population 2: Male (`"m"`), Female (`"f"`), or Total (`"t"`) assumed same as `sex1` unless otherwise specified.
#' @param closeout logical. Default `TRUE`. Shall we use the HMD Method Protocol to close out the `ax` and `qx` values? See details.
#' @details setting `closeout` to `TRUE` will result in value of `1/mx` for the final age group, of `ax` and a value of 1 for the closeout of `qx`.
#' @return `cc` numeric vector with one element per age group, and which sums to the total difference in life expectancy between population 1 and 2.
#' @importFrom data.table shift
#' @export
#' @references
#' \insertRef{arriaga1984measuring}{coddecomp}
#' \insertRef{preston2000demography}{coddecomp}
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

arriaga <- function(mx1, 
                    mx2, 
                    age = 0:(length(mx1) - 1), 
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
  
  direct = lx1 * (Lx2 / lx2 - Lx1 / lx1)
  indirect = shift(Tx2, n = -1, fill = 0) * (lx1 / lx2 - shift(lx1, n = -1, fill = 0) / shift(lx2, n = -1, fill = 0))
  N <- length(mx1)
  indirect[N] = lx1[N] * (Tx2[N] / lx2[N] - Tx1[N] / lx1[N])

  cc = direct + indirect
  return(cc)
}

#' @title the sensitivity implied by a classic Arriaga decomposition
#' @description The sensitivity of life expectancy to a perturbation in mortality rates can be derived by dividing the Arriaga decomposition result \eqn{\Delta} by the difference `mx2-mx1`.
#' \deqn{s_{x} = \frac{\Delta}{_{n}M^{2}_x - _{n}M^{1}_x}}

#' @seealso \code{\link{arriaga}}
#' @inheritParams arriaga
#' @return `s` numeric vector with one element per age group, and which gives the sensitivity values for each age.
#' @importFrom data.table shift
#' @export
#' @references
#' \insertRef{arriaga1984measuring}{coddecomp}
#' \insertRef{preston2000demography}{coddecomp}
#' @examples
#' a <- .001
#' b <- .07
#' x <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' cc <- arriaga(mx1, mx2, age = x)
#' # examples can come from above too
#' s <- sen_arriaga(mx1, mx2, age = x)
#' \dontrun{
#' plot(x, s)
#' }
#' cc_check <- s * (mx2 - mx1)
#' \dontrun{
#' plot(x,cc)
#' lines(x,cc_check)
#' }

sen_arriaga <- function(mx1,
                        mx2, 
                        age = 0:(length(mx1)-1), 
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
  
  direct = lx1 * (Lx2 / lx2 - Lx1 / lx1)
  N = length(mx1)
  indirect = shift(Tx2, n = -1, fill = 0) * (lx1 / lx2 - shift(lx1, n = -1, fill = 0) / shift(lx2, n = -1, fill = 0))
  indirect[N] = lx1[N] * (Tx2[N] / lx2[N] - Tx1[N] / lx1[N])
  #indirect = ifelse(is.na(indirect),0,indirect),
  cc = direct + indirect
  delta = mx2 - mx1
  # watch out for 0s in delta denominator
  sen = cc / delta 
  sen
}

# well, pseudo-instantaneous anyway, and symmetrical 
# This is designed to work with ltre()

#' @title Estimate sensitivity of life expectancy for a set of mortality rates
#' @description This implementation gives a good approximation of the sensitivity of life expectancy to perturbations in mortality rates (central death rates). Since the Arriaga approach requires two versions of mortality rates `mx1`, `mx2`, we create these by slightly perturbing `mx` up and down. Then we calculate the Arriaga sensitivity in each direction and take the average. Specifically, we create `mx1` and `mx2` as
#' \deqn{m_{x}^{1} = m_x \cdot \left(\frac{1}{1 - h}\right)}
#' \deqn{m_{x}^{2} = m_x \cdot \left(1 - h\right)}
#' where `h` is small value given by the argument `perturb`.
#' @details A minor correction might be needed for the final age group for the case of the reverse-direction Arriaga sensitivity. Note also for values of `perturb` (h) that are less than `1e-7` we might lose stability in results.
#' @inheritParams arriaga
#' @param mx numeric vector of mortality rates (central death rates)
#' @param sex character Male (`"m"`), Female (`"f"`), or Total (`"t"`)
#' @param perturb numeric constant, a very small number
#' @importFrom data.table shift
#' @export
#' @examples
#' a   <- .001
#' b   <- .07
#' x   <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' mx  <- (mx1 + mx2) / 2
#' s     <- sen_arriaga_instantaneous(mx, age = x)
#' s1    <- sen_arriaga_instantaneous(mx1, age = x)
#' s2    <- sen_arriaga_instantaneous(mx2, age = x)
#' s1_d  <- sen_arriaga(mx1, mx2, age = x)
#' s2_d  <- sen_arriaga(mx2, mx1, age = x)
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
sen_arriaga_instantaneous <- function(mx, 
                                      age = 0:(length(mx1)-1), 
                                      sex = 't', 
                                      perturb = 1e-6, 
                                      closeout = TRUE){
  mx1 <- mx * (1 / (1 - perturb))
  mx2 <- mx * (1 - perturb) / 1
  s1 <- sen_arriaga(mx1 = mx1, 
                    mx2 = mx2, 
                    age = age, 
                    sex1 = sex, 
                    sex2 = sex,
                    closeout = closeout)
  s2 <- sen_arriaga(mx1 = mx2, 
                    mx2 = mx1, 
                    age = age, 
                    sex1 = sex, 
                    sex2 = sex,
                    closeout = closeout)
  # TR: this might need revision, 
  # due to a discovery in the examples of arriaga()
  if (closeout){
    s2[length(s2)] <- s2[length(s2)] * 2
  }
  (s1 + s2) / 2
}

#' @title Estimate sensitivity of life expectancy for a set of mortality rates by perturbing in the log space.
#' @description This is a second approach for estimating the sensitivity for a single set of rates. Here, rather than directly expanding and contracting rates to convert `mx` into `mx1` and `mx2` we instead shift the logged mortality rates up and down by the factor `perturb = h`. Specifically:
#' \deqn{m_{x}^{1} = e^{\ln\left(m_x\right) + h}}
#' \deqn{m_{x}^{2} = e^{\ln\left(m_x\right) - h}}
#' @export
#' @inheritParams sen_arriaga_instantaneous
#' @seealso \code{\link{sen_arriaga_instantaneous}}
#' @examples
#' a   <- .001
#' b   <- .07
#' x   <- 0:100
#' mx <- a * exp(x * b)
#' # the multiplicative perturbation:
#' s1 <- sen_arriaga_instantaneous(mx)
#' s2 <- sen_arriaga_instantaneous2(mx)
#' plot(x, 
#'      s1 - s2, 
#'      pch = 16, 
#'      cex=.5,
#'      main = "very similar.\nMaybe one is more stable than the other?")
sen_arriaga_instantaneous2 <- function(mx, 
                                       age = 0:(length(mx1)-1), 
                                       sex = 't', 
                                       perturb = 1e-6, 
                                       closeout = TRUE){
  mx1 <- exp(log(mx) + perturb)
  mx2 <- exp(log(mx) - perturb)
  s1 <- sen_arriaga(mx1 = mx1, 
                    mx2 = mx2, 
                    age = age, 
                    sex1 = sex, 
                    sex2 = sex,
                    closeout = closeout)
  s2 <- sen_arriaga(mx1 = mx2, 
                    mx2 = mx1, 
                    age = age, 
                    sex1 = sex, 
                    sex2 = sex,
                    closeout = closeout)
  # TR: this might need revision, 
  # due to a discovery in the examples of arriaga()
  if (closeout){
    s2[length(s2)] <- s2[length(s2)] * 2
  }
  (s1 + s2) / 2
}

#' @title Estimate sensitivity of life expectancy using a symmetrical Arriaga approach.
#' @description This approach conducts a classic Arriaga decomposition in both directions, averaging the (sign-adjusted) result, i.e. `a_avg = (arriaga(mx1,mx2, ...) - arriaga(mx2, mx1, ...)) / 2`. 
#' @export
#' @inheritParams arriaga
#' @seealso \code{\link{arriaga}}
#' @examples
#' a <- .001
#' b <- .07
#' x <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' d <- arriaga_sym(mx1, mx2, age = x) 
#' 
#' e01 <- mx_to_e0(mx1,age=x)
#' e02 <- mx_to_e0(mx2,age=x)
#' (Delta <- e02 - e01)
#' sum(d)
#' 
#' d12 <- arriaga(mx1, mx2, age = x)
#' d21 <- arriaga(mx2, mx1, age = x) # direction opposite
#' \dontrun{
#' plot(x, d, type= 'l')
#'   lines(x, d12, col = "red")
#'   lines(x, -d21, col = "blue")
#' }
arriaga_sym <- function(mx1, 
                        mx2, 
                        age = 0:(length(mx1) - 1), 
                        sex1 = 't', 
                        sex2 = sex1, 
                        closeout = TRUE){
  a1 <- arriaga(mx1, 
                mx2, 
                age = age, 
                sex1 = sex1, 
                sex2 = sex2, 
                closeout = closeout)
  a2 <- arriaga(mx2, 
                mx1, 
                age = age, 
                sex1 = sex2, 
                sex2 = sex1, 
                closeout = closeout)
  
  # This closeout adjustment is necessary, but I can't
  # say I fully understand why this works out.
  a2[length(a2)] <- a2[length(a2)] /2
  a_avg <- (a1 - a2) / 2
  a_avg
}



#' @title Estimate sensitivity of life expectancy using a symmetrical Arriaga approach.
#' @description This approach conducts a classic Arriaga decomposition in both directions, averaging the (sign-adjusted) result, i.e. `a_avg = (arriaga(mx1,mx2, ...) - arriaga(mx2, mx1, ...)) / 2`, then approximates the sensitivity by dividing out the rate differences, i.e. `s = a_avg / (mx2 - mx1)`. A resulting decomposition will be exact because the two arriaga directions are exact, but this method might be vulnerable to 0s in the denominator.
#' @export
#' @inheritParams arriaga
#' @seealso \code{\link{arriaga}}
#' @examples
#' a <- .001
#' b <- .07
#' x <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' s <- sen_arriaga_sym(mx1, mx2, age = x) 
#' 
#' e01 <- mx_to_e0(mx1,age=x)
#' e02 <- mx_to_e0(mx2,age=x)
#' (Delta <- e02 - e01)
#' deltas <- mx2- mx1
#' sum(deltas * s)
#' 
#' \dontrun{
#' mx_avg <- (mx1 + mx2) / 2
#' plot(x, s, type= 'l')
#' lines(x, sen_arriaga_instantaneous(mx_avg, age=x),col = "blue")
#' }

sen_arriaga_sym <- function(mx1, 
                            mx2, 
                            age = 0:(length(mx1) - 1), 
                            sex1 = 't', 
                            sex2 = sex1, 
                            closeout = TRUE){
 delta <- mx2 - mx1
 a_avg <- arriaga_sym(mx1 = mx1,
                      mx2 = mx2,
                      age = age,
                      sex1 = sex1,
                      sex2 = sex2,
                      closeout = closeout)
 a_avg / delta
}



