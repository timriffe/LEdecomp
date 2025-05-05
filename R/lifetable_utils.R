
#' @title top-down cumulative sums, as in the lifetable Tx
#' @description Why write `x |> rev() |> cumsum() |> rev()` when you can just write `rcumsum(x)`?
#' @param x numeric vector
#' @return numeric vector the same length as `x`
#' @export
rcumsum <- function(x){
  rev(cumsum(rev(x)))
}

#' @title produce single-age `ax` values
#' @description We assume mid-interval `ax` except for age 0 and potentially the open age group. `ax` is defined as the average years lived in each age interval by those that die within the interval, and it is used to increase the precision of lifetable estimates. We allow ourselves the midpoint rule for single ages because it has little leverage. If we were working with abridged ages then we would need to use a more sophisticated method.
#' @details For the case of Total sex, we estimate the male and female \eqn{a(0)} using the Andreev-Kingkade rule of thumb, and then average them. We assume a value of 1/2 for all other ages, unless `closeout = TRUE`, in which case we close with `1/mx` for the final value.
#' @param mx numeric vector of the mortality rates (central death rates)
#' @param age integer vector of the lower bound of each age group (currently only single ages supported)
#' @param sex character: Male (`"m"`), Female (`"f"`), or Total (`"t"`)
#' @param closeout logical. Default `TRUE`.
#' @importFrom DemoTools lt_rule_ak_m0_a0
#' @export
#' @seealso \link[DemoTools]{lt_rule_ak_m0_a0}
mx_to_ax <- function(mx, age = 0:(length(mx)-1), sex = "t", closeout = TRUE){
  stopifnot(all(diff(age) == 1))
  ax <- rep(.5, length(mx))
  sex <- sex |> tolower() |> substr(1,1)
  sex <- ifelse(sex == "b","t",sex)
  stopifnot(sex %in% c("m","f","t"))
  if (min(age) == 0){
    if (sex != "t"){
      a0 <- lt_rule_ak_m0_a0(M0 = mx[1], Sex = sex)
    } else {
      a0m <- lt_rule_ak_m0_a0(M0 = mx[1], Sex = "m")
      a0f <- lt_rule_ak_m0_a0(M0 = mx[1], Sex = "f")
      a0  <- (a0m + a0f) / 2
    }
    ax[1] <- a0
  }
  if (closeout){
    ax[length(ax)] <- 1 / mx[length(mx)]
  }
  ax
}

#' @title produce single-age `qx` values
#' @description `qx` gives conditional death probabilities, in this case forced to be consistent with a set of `mx` and `ax` values per HMD Method Protocol eq 71.
#' @inheritParams mx_to_ax
#' @param ax numeric vector of `ax` values
#' @export
#' @seealso \link[DemoTools]{lt_id_ma_q}
#' @importFrom DemoTools lt_id_ma_q
#' @references
#' \insertRef{wilmoth2021methods}{LEdecomp}
mx_to_qx <- function(mx, ax, closeout = TRUE){
  lt_id_ma_q(nMx = mx,
             nax = ax,
             AgeInt = rep(1, length(mx)),
             closeout = closeout)
}

#' @title Calculate the survival curve
#' @description The survival curve is calculated as the cumulative product of the conditional survival probabilities, which are the complement of conditional death probabilities, `qx`, except we take care to start with a clean 1. This function no radix option. `lx` with a radix of 1 can be interpreted as the probability of surviving from birth to age `x`.
#' @param qx numeric vector of conditional death probabilities
#' @return numeric vector of `lx` values
#' @importFrom DemoTools lt_id_q_l
#' @seealso \link[DemoTools]{lt_id_q_l}
#' @export
qx_to_lx <- function(qx){
  lt_id_q_l(nqx = qx, radix = 1)
}

#' @title Calculate the lifetable death distribution
#' @description Minus the decumulation of the survival curve gives the death distribution. Or the element-wise product of `lx` and the conditional death probabilities `qx` gives the same thing.
#' @param lx numeric vector of lifetable survivors at exact age `x`
#' @return numeric vector of `dx` values
#' @importFrom DemoTools lt_id_l_d
#' @seealso \link[DemoTools]{lt_id_l_d}
#' @export

lx_to_dx <- function(lx){
  lt_id_l_d(lx)
}

#' @title Calculate the lifetable exposure
#' @description `Lx` is defined as the integration of `lx` in the interval `[x,x+n)`, where `n` is the width of the interval. There are many approximations for this. Here we use HMD Method Protocol equation 78. You can think of `Lx` as lifetable exposure, or person-years lived in each age interval.
#' @param ax numeric vector of ax, average time spent in the age interval by those that die in the interval
#' @param lx numeric vector of lx, lifetable survivorship at exact ages.
#' @param dx  numeric vector of dx, the lifetable deaths distribution.
#' @return numeric vector of `Lx` values
#' @importFrom DemoTools lt_id_lda_L
#' @seealso \link[DemoTools]{lt_id_lda_L}
#' @export
#' @references
#' \insertRef{wilmoth2021methods}{LEdecomp}
ald_to_Lx <- function(ax,lx,dx){
  lt_id_lda_L(nax = ax,
              lx = lx,
              ndx = dx,
              AgeInt = rep(1,length(ax)))
}

#' @title calculate remaining life expectancy `ex` for each age
#' @description Here we combine HMD Method Protocol equations 79 and 80. We calculate all the remaining years left to live at each age, then condition this on survival to each age.
#' @param lx numeric vector of lifetable survivors at exact age `x`
#' @param Lx numeric vector of lifetable exposure `Lx`
#' @return numeric vector of remaining life expectancy `ex`
#' @export
#' @references
#' \insertRef{wilmoth2021methods}{LEdecomp}
lL_to_ex <- function(lx, Lx){
  Tx <- rcumsum(Lx)
  ex <- Tx / lx
  ex
}


#' @title calculate remaining life expectancy from mortality rates
#' @description We follow the full chain of standard lifetable column calculations to translate `mx` to `ex`.
#' @inheritParams mx_to_ax
#' @return numeric vector of `ex`, the same length as `mx`
#' @export

mx_to_ex <- function(mx, age, sex = 't', closeout = TRUE){
  ax <- mx_to_ax(mx = mx,
                 age = age,
                 sex = sex,
                 closeout = closeout)
  qx <- mx_to_qx(mx = mx,
                 ax = ax,
                 closeout = closeout)
  lx <- qx_to_lx(qx)
  dx <- lx_to_dx(lx)
  Lx <- ald_to_Lx(ax = ax,
                  lx = lx,
                  dx = dx)
  ex <- lL_to_ex(lx = lx, Lx = Lx)
  ex
}

#' @title calculate life expectancy at birth from mortality rates
#' @description We follow the full chain of standard lifetable column calculations to translate `mx` to `ex`, then select the first element of `ex`. If `min(age) > 0`, then we return remaining life expectancy at the lowest given age.
#' @inheritParams mx_to_ax
#' @return numeric scalar of `e0`
#' @export

mx_to_e0 <- function(mx,
                     age,
                     sex = 't',
                     closeout = TRUE){
  ex <- mx_to_ex(mx = mx,
                 age = age,
                 sex = sex,
                 closeout = closeout)
  ex[1]
}
