
library(devtools)
library(usethis)
library(rhub)
document()

devtools::check()
check_win_devel()      #
check_win_release()    #
check_win_oldrelease() #
rhub_platforms()
rhub::rhub_check(platforms = c("linux","macos","macos-arm64","windows"))

library(spelling)
spell_check()

# doesn't work until on CRAN?
#library(revdepcheck)
#revdepcheck::revdep_check(cran = FALSE)

devtools::release()
