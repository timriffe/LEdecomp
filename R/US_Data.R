#' US Mortality data
#' @description Data from the US total population from the Human Mortality Dataset and from National Center for Health Statistics (NCHS). The use of two dataset is justified because NCHS does not contain information of exposures above age 85. The dataset contains information on mortality rates (mxt), registered deaths (Dxt) and the size of the population at risk of death (Ext) by period, from 2000 to 2020, and by age, from 0 to 100 years, for both males and females.
#' @name US_data
#' @format A data frame with 4242 rows and 6 columns with class \code{"LEdecompData"} and \code{"data.frame"} including the following information
#' * `Age` a vector containing the ages considered in the dataset, 0, 1, ..., 99, and 100.
#' * `Gender` a vector containing the information regarding the gender, \code{"Male"} or \code{"Female"}.
#' * `Period` a vector containing the periods of the dataset from 2000 to 2020.
#' * `Ext` a vector containing the size of the population at risk of death by age and period.
#' * `Dxt` a vector containing the number of registered deaths by age and period.
#' * `mxt` a vector mortality rates for the corresponding age and period.
#' @docType data
#' @usage US_data
#' @examples
#' #The dataset is executed with the following information
#' US_data
"US_data"

#' @export
print.LEdecompData <- function(x, ...){
  message("Mortality Data\n")
  message(attributes(x)$label, "including", attributes(x)$series ,"\n")
  message("Periods", c(min(x$Period),":", max(x$Period)),"\n")
  message("Complete Life Table with ages from", c(min(x$Age)," to ", max(x$Age)), "\n")
}
