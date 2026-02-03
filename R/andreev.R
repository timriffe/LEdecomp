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
                    age,
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
