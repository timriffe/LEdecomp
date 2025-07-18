% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arriaga.R
\name{arriaga_sym}
\alias{arriaga_sym}
\title{Estimate sensitivity of life expectancy using a symmetrical Arriaga approach.}
\usage{
arriaga_sym(
  mx1,
  mx2,
  age = 0:(length(mx1) - 1),
  nx = rep(1, length(mx1)),
  sex1 = "t",
  sex2 = sex1,
  closeout = TRUE
)
}
\arguments{
\item{mx1}{numeric vector of the mortality rates (central death rates) for population 1}

\item{mx2}{numeric vector of the mortality rates (central death rates) for population 2}

\item{age}{integer vector of the lower bound of each age group (currently only single ages supported)}

\item{nx}{integer vector of age intervals, default 1.}

\item{sex1}{character either the sex for population 1: Male (\code{"m"}), Female (\code{"f"}), or Total (\code{"t"})}

\item{sex2}{character either the sex for population 2: Male (\code{"m"}), Female (\code{"f"}), or Total (\code{"t"}) assumed same as \code{sex1} unless otherwise specified.}

\item{closeout}{logical. Default \code{TRUE}. Shall we use the HMD Method Protocol to close out the \code{ax} and \code{qx} values? See details.}
}
\description{
This approach conducts a classic Arriaga decomposition in both directions, averaging the (sign-adjusted) result, i.e. \code{a_avg = (arriaga(mx1,mx2, ...) - arriaga(mx2, mx1, ...)) / 2}.
}
\examples{
a <- .001
b <- .07
x <- 0:100
mx1 <- a * exp(x * b)
mx2 <- a/2 * exp(x * b)
d <- arriaga_sym(mx1, mx2, age = x)

e01 <- mx_to_e0(mx1,age=x)
e02 <- mx_to_e0(mx2,age=x)
(Delta <- e02 - e01)
sum(d)

d12 <- arriaga(mx1, mx2, age = x)
d21 <- arriaga(mx2, mx1, age = x) # direction opposite
\dontrun{
plot(x, d, type= 'l')
  lines(x, d12, col = "red")
  lines(x, -d21, col = "blue")
}
}
\seealso{
\code{\link{arriaga}}
}
