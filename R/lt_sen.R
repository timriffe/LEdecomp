#

#' @title A direct approximation of the sensitivity of life expectancy at birth to changes in mortality.
#' @description This function tries to get the direct discrete life expectancy sensitivity to \eqn{m(x)}, in continous math it's \eqn{-l(x)e(x)}, we just need to find the best approx with a discrete lifetable.
#' This direct lifetable-based calculation requires a few approximations to get a usable value whenever we're working with discrete data.
#' In continuous notation, we know that the sensitivity \eqn{s(x)}
#' \deqn{s(x) = -l(x)e(x)}
#' but it is not obvious what to use from a discrete lifetable. In this implementation, we use \eqn{L(x)} and an \eqn{a(x)}-weighted average of successive \eqn{e(x)} values, specifically, we calculate:
#' \deqn{
#' s_x = -L_x \cdot \left e_x \cdot \left( 1 - a_x \right) + e_{x+1} \cdot a_x\right)
#' }
#' This seems to be a very good approximation for ages >0, but we still have a small, but unaccounted-for discrepancy in age 0, at least when comparing with also-imperfect numerical derivatives.
#' @inheritParams mx_to_e0
#' @importFrom data.table shift
#' @importFrom numDeriv grad
#' @importFrom Rdpack reprompt
#' @export
#' @examples
#' x <- 0:100
#' mx <- 0.001 * exp(x * 0.07)
#' sl <-  sen_e0_mx_lt(mx,age=x,sex='t',closeout=TRUE)
#' sn <- numDeriv::grad(mx_to_e0, mx, age=x, sex = 't', closeout=TRUE)
#' \dontrun{
#' plot(x,sl)
#' lines(x,sn)
#' }
#' # examine residuals:
#' sl - sn
#' # Note discrepancies in ages >0 are due to numerical precision only
#' \dontrun{
#' plot(x, sl - sn, main = "still uncertain what accounts for the age 0 discrepancy")
#' }
sen_e0_mx_lt <- function(mx,
                      age = 0:(length(mx)-1),
                      sex = 't',
                      closeout = TRUE){
  ax <- mx_to_ax(mx = mx,
                 age = age,
                 sex = sex,
                 closeout = closeout)
  qx <- mx_to_qx(mx, ax)
  lx <- qx_to_lx(qx)
  dx <- lx_to_dx(lx)
  Lx <- ald_to_Lx(ax = ax,
                  lx = lx,
                  dx = dx)
  ex <- lL_to_ex(lx = lx, Lx = Lx)
  # ex <- mx_to_ex(mx)
  N <- length(mx)

  # This is the current-best approximation,
  # but still requires an unknown adjustment for age 0
  exs <- shift(ex, n = -1, fill = ex[N])
  ex2 <- ex * (1 - ax) + exs * ax

  ex2[N] <- ex[N]

  sen <- -Lx * ex2

  sen
}

#' @title A numerical approximation of the sensitivity of life expectancy at birth to changes in mortality.
#' @description Here we produce a numerical derivative based on the methods implemented in the `numDeriv::grad()` function. Tweaking the optional arguments of `numDeriv::grad()`, passed in via `...` might lead to greater precision, but this method actually performs usably well with its defaults.
#' @inheritParams mx_to_e0
#' @param ... optional args to pass to `mx_to_e0()`
#' @importFrom numDeriv grad
#' @importFrom Rdpack reprompt
#' @export
#' @examples
#' x <- 0:100
#' mx <- 0.001 * exp(x * 0.07)
#' sn <- sen_num(mx,age=x,sex='t',closeout=TRUE)
#' sa <- sen_arriaga_instantaneous2(mx, age=age,sex='t',perturb = 1e-4)
#' \dontrun{
#' plot(x,sa)
#' lines(x,sn)
#' }
#' # examine residuals:
#' sn - sa
#' # Note discrepancies in ages >0 are due to numerical precision only
#' \dontrun{
#' plot(x, sn - sa, main = "still uncertain what accounts for the age 0 discrepancy")
#' }
sen_num <- function(mx,age,sex='t',closeout=TRUE,...){
  numDeriv::grad(mx_to_e0,mx,age=age,sex=sex,closeout=closeout,...)
}

 #  x <- 0:100
 # mx <- 0.001 * exp(x * 0.07)
 # richardson_hack <- function(mx, age, sex = 't', closeout = TRUE, r = 5, d1 = 1e-5 ...){
 #
 #   pert <- matrix(0, nrow = length(mx), ncol = r)
 #   d <- d1
 #   for (rr in 1:r){
 #     d = d / 2
 #     for (i in 1:length(mx)){
 #       mxi1 <- mxi2 <- mx
 #       mxi1[i] <- mxi1[i] + d
 #       mxi2[i] <- mxi2[i] - d
 #       pert[i, rr] <-
 #       (mx_to_e0(mx = mxi1, age = age, sex = sex, closeout = closeout) -
 #         mx_to_e0(mx = mxi2, age = age, sex = sex, closeout = closeout)) /
 #         (2 * d)
 #     }
 #   }
 #
 #
 #
 #  }
