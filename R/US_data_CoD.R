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
#' @format A data frame of of mortality rates by age, year, sex, and cause with 76356 rows and 6 columns, including the following columns. Cause of death fractions were derived from NCHS data, and constrained to HMD lifetable mx:
#' * `year` a vector containing the periods of the dataset from 2000 to 2020.
#' * `sex` a vector containing the information regarding the gender, `"Male"` or `"Female"`.
#' * `age` a vector containing the ages considered in the dataset, 0, 1, ..., 99, and 100.
#' * `cause` a vector containing a brief summary of the corresponding cause of death.
#' * `cause_id` a vector containing the corresponding identification number for the cause of death.
#' * `mxc` a vector mortality rates for the corresponding age and cause.
#' @docType data
#' @usage US_data_CoD
#' @references
#'   \insertRef{hmd2026}{LEdecomp}
#'   \insertRef{cdcwonder2024}{LEdecomp}
#' @examples
#' #The dataset is loaded by simply executing:
#' US_data_CoD
#'
"US_data_CoD"




#  library(HMDHFDplus)
#  mlt <- readHMDweb("USA","mltper_1x1", username = Sys.getenv("us"), password = Sys.getenv("pw"))|>
#    select(year = Year, age = Age, mx) |>
#    mutate(sex = "Male", .after = year)
#  flt <- readHMDweb("USA","fltper_1x1", username = Sys.getenv("us"), password = Sys.getenv("pw"))|>
#    select(year = Year, age = Age, mx) |>
#    mutate(sex = "Female", .after = year)
# #
#  hmd <- bind_rows(mlt,flt)
#  US_data_CoD <- US_data_CoD |>
#    left_join(hmd, by = join_by(year,sex,age)) |>
#    group_by(year,sex,age) |>
#    mutate(mxc = (mxc / sum(mxc)) * mx)
#
# US_data_CoD |> usethis::use_data(overwrite = TRUE)

