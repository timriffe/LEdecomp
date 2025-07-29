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
  # TR: is this closeout actually needed here?
  # I have my doubts this was checked
  #To match the solution with the arriaga
  # if (closeout){
  #   s2[length(s2)] <- s2[length(s2)] * 2
  # }
  sen <- (s1 + s2) / 2

  sen
}
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
