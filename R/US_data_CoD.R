#' US cause-of-death Mortality data
#'
#' Data from the US total population from the Human Mortality Dataset and from National Center for Health Statistics (NCHS).
#' In this case, we have information with the number of deaths by 18 different causes, for more information please review the cause_id and the National Center for Health Statistics (NCHS).
#' The use of two dataset is justified because NCHS does not contain information of exposures above age 85.
#' The dataset contains information on mortality rates (mxt), registered deaths (Dxt) and the size of the population at risk of death (Ext) by period, from 2000 to 2020, and by age, from 0 to 100 years, for both males and females.
#' In addition, we have the number of deaths by cause between 0 and 100 years of age and between 2000 and 2020.
#'
#' @name US_data_CoD
#'
#' @format A data frame with 80598 rows and 8 columns with class \code{"LEdecompData"} and \code{"data.frame"} including the following information
#' * `Age` a vector containing the ages considered in the dataset, 0, 1, ..., 99, and 100.
#' * `Gender` a vector containing the information regarding the gender, \code{"Male"} or \code{"Female"}.
#' * `Period` a vector containing the periods of the dataset from 2000 to 2020.
#' * `Ext` a vector containing the size of the population at risk of death by age and period.
#' * `Dxt` a vector containing the number of registered deaths by age and period.
#' * `mxt` a vector mortality rates for the corresponding age and period.
#' * `cause` a vector containing a brief summary of the corresponding cause of death.
#' * `cause_id` a vector containing the corresponding identification number for the cause of death.
#' @docType data
#' @usage US_data_CoD
#'
#' @examples
#' #The dataset is executed with the following information
#' US_data_CoD
#'
"US_data_CoD"
