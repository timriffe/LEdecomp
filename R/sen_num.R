
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
#' sa <- sen_arriaga_instantaneous2(mx, age=x,sex='t',perturb = 1e-4)
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
sen_num <- function(mx,
                    age = (1:length(mx))-1,
                    nx = rep(1,length(mx)),
                    sex = 't',
                    closeout = TRUE,...){
  numDeriv::grad(mx_to_e0,
                 mx,
                 age = age,
                 nx = nx,
                 sex = sex,
                 closeout = closeout,...)
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
