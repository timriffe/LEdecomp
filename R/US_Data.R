#' US Mortality data
#' @description All-cause mortality rate data from the Human Mortality Database lifetables. The dataset contains information on mortality rates (mx), by year (2000-2020), sex, and (truncated) ages 0 to 100.
#' @name US_data
#' @format A data frame with 4242 rows and 4 columns with class  `"data.frame"` and  including the following columns
#' * `year` integer. Years 2000 to 2020.
#' * `sex` character. \code{"Male"} or \code{"Female"}.
#' * `age` integer. Ages 0, 1, ..., 99, and 100.
#' * `mx` numeric. Mortality rates for the corresponding sex, age and year
#' @docType data
#' @usage US_data
#' @references
#'   \insertRef{hmd2026}{LEdecomp}
#' @examples
#' #The dataset is executed with the following information
#' US_data
"US_data"

