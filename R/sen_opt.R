
#' @title sen_resid 
#' @description Most sensitivity methods in this packages (`sen_arriaga_sym()` excepted) are approximations; when used in decompositions they will tend to imply residuals. To acheive near-exact additivity for a decomposition using these sensitivity approaches, one can try to find a different weighting of rates from populations 1 and 2, rather than simply taking their arithmetic average. Here we turn this into an optimization problem, where we find the weighting `w` that implies an exactly additive decomposition to an arbitrary degree of tolerance. This function gives said residual, for purposes of optimizing using `sen_min()`. We export this auxiliary function because one might wish to know the value w that balances rates such that the decomposition is exact. 
#' \deqn{m_{x} = m_{x}^{1} * w + m_{x}^{2} * (1-w)}
#' @inheritParams arriaga
#' @param tol double. tolerance level for residual, passed to `optimise()`
#' @param sen_fun function name, current options include `sen_arriaga_instantaneous`, `sen_arriaga_instantaneous2`, `sen_arriaga_sym`, `sen_e0_mx_lt`, `sen_num`
#' @return age-specific sensitivity of life expectancy to changes in mortality rates.
#' @export
#' @examples
#' a   <- .001
#' b   <- .07
#' x   <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' w <- optimize(sen_resid,
#'               mx1 = mx1,
#'               mx2 = mx2,
#'               age = x,
#'               sen_fun = sen_arriaga_instantaneous,
#'               sex1 = 't',
#'               sex2 = 't',
#'               closeout = TRUE, 
#'               interval = c(.4,.6),
#'               tol = 1e-10)$minimum
#' w
#' 
sen_resid <- function(w=.5,
                    mx1,
                    mx2,
                    age,
                    sex1,
                    sex2 = sex1,
                    closeout = TRUE,
                    sen_fun = sen_arriaga_instantaneous, ...){
  
  mx    <- w * mx1 + (1-w) * mx2
  sex   <- ifelse(sex1 != sex2,"t",sex1)
  s     <- sen_fun(mx, age = age, sex = sex, closeout = closeout, ...)
  delta <- mx2 - mx1
  e02   <- mx_to_e0(mx=mx2, age=age,sex=sex2,closeout=closeout)
  e01   <- mx_to_e0(mx=mx1, age=age,sex=sex1,closeout=closeout)
  gap   <- e02 - e01
  gap_hat <- sum(s * delta)
  abs(gap - gap_hat)
}

#' @title sen_min 
#' @description Most sensitivity methods in this packages (`sen_arriaga_sym()` excepted) are approximations; when used in decompositions they will tend to imply residuals. To acheive near-exact additivity for a decomposition using these sensitivity approaches, one can try to find a different weighting of rates from populations 1 and 2, rather than simply taking their arithmetic average. Here we turn this into an optimization problem, where we find the weighting `w` that implies an exactly additive decomposition to an arbitrary degree of tolerance.
#' \deqn{m_{x} = m_{x}^{1} * w + m_{x}^{2} * (1-w)}
#' @details We expect the value `w` to be close to .5, and only search the interval `[.4,.6]`. This may need to be revisited in case that proves too narrow.
#' @inheritParams arriaga
#' @param tol double. tolerance level for residual, passed to `optimise()`
#' @param sen_fun function name, current options include `sen_arriaga_instantaneous`, `sen_arriaga_instantaneous1`, `sen_arriaga_sym`, `sen_e0_mx_lt`,  `sen_num`
#' @return age-specific sensitivity of life expectancy to changes in mortality rates.
#' @export
#' @examples
#' a   <- .001
#' b   <- .07
#' x   <- 0:100
#' mx1 <- a * exp(x * b)
#' mx2 <- a/2 * exp(x * b)
#' mx  <- (mx1 + mx2) / 2
#' s1 <- sen_min(mx1, mx2,
#'               age = x, sex1 = 't',
#'               closeout = TRUE, 
#'               sen_fun = sen_arriaga_instantaneous)
#' s2 <- sen_min(mx1, mx2,
#'               age = x, sex1 = 't',
#'               closeout = TRUE, 
#'               sen_fun = sen_e0_mx_lt,
#'               tol = 1e-12)
#' 
#' # check sums
#' e01 <- mx_to_e0(mx1,x,'t',TRUE)
#' e02 <- mx_to_e0(mx2,x,'t',TRUE)
#' (gap <- e02 - e01)
#' delta <- mx2 - mx1
#' (gap1 <- sum(s1 * delta))
#' (gap2 <- sum(s2 * delta))
#' gap2-gap
#' 
#' plot(x, s1, type= 'l')
#' lines(x, s2, col = 'red', lty = 2, lwd = 2)
#' 
#' plot(x, s2-s1, main = "age 0 difference is due to imprecision\nin lifetable approach #' for this age")

sen_min <- function(mx1,
                    mx2,
                    age,
                    sex1,
                    sex2 = sex1,
                    closeout = TRUE,
                    sen_fun = sen_arriaga_instantaneous, 
                    tol = 1e-10,
                    ...){
  w <- optimize(sen_resid,
           mx1 = mx1,
           mx2 = mx2,
           age = x,
           sen_fun = sen_fun,
           sex1 = sex1,
           sex2 = sex2,
           closeout = TRUE, 
           interval = c(.4,.6),
           tol = tol,
           ...)$minimum
  mx    <- w * mx1 + (1-w) * mx2
  sex   <- ifelse(sex1 != sex2,"t",sex1)
  s     <- sen_fun(mx, age = age, sex = sex, closeout = closeout, ...)
  s
}
