# library(DemoDecomp)
#  a <- 0.001
#  b <- 0.07
#  age <- 0:100
#  mx <- a * exp(age * b)
#  sex = 't'
#  perturb = 1e-5
sen_horiuchi <- function(mx, age = 0:(length(mx)-1), sex = 't', perturb = 1e-5){
   mx1 <- mx * (1 / (1 - perturb))
   mx2 <- mx * ((1 - perturb) / 1)

   cc <- horiuchi(mx_to_e0, mx1, mx2, age = age, sex = sex, N = 20)
   sen <- cc / (mx2 - mx1)
   sen
}
sex = 't'
horiuchi(mx_to_e0, mx1, mx2, age = age, sex = sex, N = 20)

# plot(age, sen_horiuchi(mx, age) - sen_arriaga_instantaneous(mx))
# plot(age, sen_horiuchi(mx, age) - numDeriv::grad(mx_to_e0,mx,age = age, sex = sex))
# plot(age[-1], sen_horiuchi(mx, age)[-1] - sen_e0_mx_lt(mx,age,sex)[-1])
horiuchi <- function(func, pars1, pars2, N, ...){

    d 			    <- pars2 - pars1
    n 			    <- length(pars1)
    delta 		  <- d / N
    grad        <- matrix(rep(.5:(N - .5) / N, n),
                          byrow = TRUE, ncol = N)
    x           <- pars1 + d * grad
    cc          <- matrix(0, nrow = n, ncol = N)
    # TR: added 16-9-2024 so that names can be used to reconstruct
    # data inside func()
    rownames(x) <- names(pars1)
    for (j in 1:N){
      DD <- diag(delta / 2)
      for (i in 1:n){
        cc[i,j] <- func((x[, j] + DD[, i]), ...) - func((x[, j] - DD[, i]), ...)
      }
    }

    out <- rowSums(cc)
    names(out) <- names(pars1)
    out
}
horiuchi_sym <- function(func, pars1, pars2, N, ...){
  a1 <- horiuchi(func, pars1, pars2, N = N, ...)
  a2 <- horiuchi(func, pars2, pars1, N = N, ...)

  # This closeout adjustment is necessary, but I can't
  # say I fully understand why this works out.
  a2[length(a2)] <- a2[length(a2)] /2
  a_avg <- (a1 - a2) / 2
  a_avg
}


#' @return the value of R0 for the given set of rates and proportion female of births.
#'
#' @export
#' @examples
#' data(rates1)
#' # take vec:
#' x <- c(rates1)
#' R0vec(x)
R0vec <- function(x, pfem = .4886){

    dim(x) <- c(length(x) / 2, 2)
    sum(x[, 1] * x[, 2] * pfem)
}
