
library(devtools)
library(usethis)
library(rhub)
document()

devtools::check()
check_win_devel()      #
check_win_release()    #
check_win_oldrelease() #


rhub_check(platforms = c("linux", "macos", "macos-arm64", "windows"))

library(spelling)
spell_check()

library(revdepcheck)
revdepcheck::revdep_check()

devtools::release()

library(tidyverse)
data("US_data")
View(US_data)
US_data |>
  filter(Period == 2020) |>
  select(Age, Gender, mx)
  pivot_wider(names_from = Gender, )
