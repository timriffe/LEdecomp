% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lifetable_utils.R
\name{mx_to_ex}
\alias{mx_to_ex}
\title{calculate remaining life expectancy from mortality rates}
\usage{
mx_to_ex(mx, age, nx = rep(1, length(age)), sex = "t", closeout = TRUE)
}
\arguments{
\item{mx}{numeric vector of the mortality rates (central death rates)}

\item{age}{integer vector of the lower bound of each age group (currently only single ages supported)}

\item{nx}{age interval width, assumes 1 by default}

\item{sex}{character: Male (\code{"m"}), Female (\code{"f"}), or Total (\code{"t"})}

\item{closeout}{logical. Default \code{TRUE}.}
}
\value{
numeric vector of \code{ex}, the same length as \code{mx}
}
\description{
We follow the full chain of standard lifetable column calculations to translate \code{mx} to \code{ex}.
}
